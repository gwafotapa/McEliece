use log::info;

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

use crate::finite_field::{F2m, Field, F2};
use crate::goppa::Goppa;
use crate::matrix::{Mat, Perm, RowVec};
// use crate::polynomial::Poly;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Debug, Eq, PartialEq)]
pub struct PublicKey<'a> {
    pub sgp: Mat<'a, F2>, // disguised generator matrix S * G * P
    pub t: u32,           // degree of the goppa polynomial
}

#[derive(Eq, PartialEq)]
pub struct SecretKey<'a, 'b> {
    pub s: Mat<'a, F2>,        // singular matrix S
    pub goppa: Goppa<'b, F2m>, // goppa code of generator matrix G
    pub info_set: Vec<usize>,  // information set of matrix G
    pub p: Perm,               // permutation matrix P
}

pub fn keygen<'a, 'b>(
    f2: &'a F2,
    f2m: &'b F2m,
    // m: u32,
    n: u32,
    t: u32,
) -> (PublicKey<'a>, SecretKey<'a, 'b>) {
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f2m, n as usize, t as usize);
    // let goppa = Goppa::new(
    //     Poly::random_irreducible(&mut rng, f2m, t as usize),
    //     f2m.as_set(),
    // )
    // .unwrap();
    let h = goppa.parity_check_matrix();
    info!("Parity check matrix:{}", h);
    let h2 = h.binary_form(f2);
    info!("Parity check matrix in binary form:{}", h2);
    let (g, info_set) = Goppa::generator_matrix(&h2);
    info!("Generator matrix G:{}", g);
    let k = g.rows();
    let s = Mat::invertible_random(&mut rng, f2, k as usize);
    info!("Singular matrix S:{}", s);
    let p = Perm::random(&mut rng, n as usize);
    let sgp = &s * &g * &p;
    info!("Perturbed generator matrix ~G:{}", sgp);
    let pk = PublicKey { sgp, t };
    println!("{:?}", info_set);
    let sk = SecretKey {
        s,
        goppa,
        info_set,
        p,
    };
    (pk, sk)
}

impl<'a> PublicKey<'a> {
    pub fn encrypt(&self, m: &RowVec<'a, F2>) -> RowVec<'a, F2> {
        let mut rng = rand::thread_rng();
        let c = m * &self.sgp;
        let z = RowVec::random_with_weight(&mut rng, m.field(), self.sgp.cols(), self.t as usize);
        c + z
    }

    pub fn save_public_key(&self, file_name: &str) -> Result<()> {
        let mut f = File::create(file_name)?;
        let s = self.sgp.to_hex_string() + format!("\n{:02x}\n", self.t).as_str();
        f.write_all(s.as_bytes())?;
        Ok(())
    }

    pub fn load_public_key(file_name: &str, f2: &'a F2) -> Result<Self> {
        let f = File::open(file_name)?;
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)?;
        line.pop(); // Remove terminating newline
        let sgp = Mat::from_hex_string(&line, f2)?;
        line.clear();
        f.read_line(&mut line)?;
        line.pop(); // Remove terminating newline
        let t = u32::from_str_radix(&line, 16)?;
        Ok(PublicKey { sgp, t })
    }
}

impl<'a, 'b> SecretKey<'a, 'b> {
    pub fn decrypt(
        &self,
        c: &RowVec<'a, F2>,
    ) -> RowVec<'a, F2> {
        let m_s_g_z = c * self.p.inverse();
        let m_s_g = self.goppa.decode(&m_s_g_z);
        let m_s = m_s_g.extract_cols(&self.info_set);
        let m = m_s * self.s.inverse().unwrap();
        m
    }

    pub fn save_secret_key(&self, file_name: &str) -> Result<()> {
        let f = File::create(file_name)?;
        let mut f = BufWriter::new(f);
        f.write_all((self.s.to_hex_string() + "\n").as_bytes())?;
        f.write_all((self.goppa.to_hex_string() + "\n").as_bytes())?;

        for i in 0..self.info_set.len() - 1 {
            f.write_all(format!("{:x} ", self.info_set[i]).as_bytes())?;
        }
        f.write_all(format!("{:x}\n", self.info_set[self.info_set.len() - 1]).as_bytes())?;

        for i in 0..self.p.cols() - 1 {
            f.write_all(format!("{:x} ", self.p[i]).as_bytes())?;
        }
        f.write_all(format!("{:x}\n", self.p[self.p.cols() - 1]).as_bytes())?;
        Ok(())
    }

    pub fn load_finite_field(file_name: &str) -> Result<F2m> {
        let f = File::open(file_name)?;
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)?;
        line.clear();
        f.read_line(&mut line)?;
        let order = u32::from_str_radix(&line[..line.find('#').ok_or("Missing hashtag")?], 16)?;
        let f = Field::generate(order);
        Ok(f)
    }

    pub fn load_secret_key(file_name: &str, f2: &'a F2, f2m: &'b F2m) -> Result<Self> {
        let f = File::open(file_name)?;
        let f = BufReader::new(f);
        let mut lines_iter = f.lines().map(|l| l.unwrap());

        let line = lines_iter.next().ok_or("Cannot read line 1")?;
        let s = Mat::from_hex_string(&line, f2)?;

        let line = lines_iter.next().ok_or("Cannot read line 2")?;
        let goppa = Goppa::from_hex_string(&line, f2m)?;

        let line = lines_iter.next().ok_or("Cannot read line 3")?;
        let v: Vec<&str> = line.split(' ').collect();
        let mut info_set = Vec::new();
        for i in 0..v.len() {
            info_set.push(usize::from_str_radix(v[i], 16)?);
        }

        let line = lines_iter.next().ok_or("Cannot read line 4")?;
        let v: Vec<&str> = line.split(' ').collect();
        let mut p = Vec::new();
        for i in 0..v.len() {
            p.push(usize::from_str_radix(v[i], 16)?);
        }
        let p = Perm::new(p)?;

        Ok(SecretKey {
            s,
            goppa,
            info_set,
            p,
        })
    }
}
