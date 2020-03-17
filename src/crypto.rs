use log::info;

use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
};

use crate::finite_field::{F2m, Field, F2};
use crate::goppa::Goppa;
use crate::matrix::{Mat, RowVec};
use crate::polynomial::Poly;

#[derive(Debug, Eq, PartialEq)]
pub struct PublicKey<'a> {
    pub sgp: Mat<'a, F2>, // disguised generator matrix S * G * P
    pub t: u32,           // degree of the goppa polynomial
}

#[derive(Eq, PartialEq)]
pub struct SecretKey<'a, 'b> {
    pub s: Mat<'a, F2>,        // singular matrix S
    pub goppa: Goppa<'b, F2m>, // goppa code of generator matrix G
    pub p: Mat<'a, F2>,        // permutation matrix P
    pub info_set: Vec<usize>,  // information set of matrix G
}

pub fn keygen<'a, 'b>(
    f2: &'a F2,
    f2m: &'b F2m,
    m: u32,
    n: u32,
    t: u32,
) -> (PublicKey<'a>, SecretKey<'a, 'b>) {
    let mut rng = rand::thread_rng();
    let goppa = Goppa::new(
        // TODO: add a random goppa method
        Poly::random_irreducible(&mut rng, f2m, t as usize),
        f2m.as_set(),
    )
    .unwrap();
    let k = n - m * t;
    let h = goppa.parity_check_matrix();
    info!("Parity check matrix:{}", h);
    let h2 = h.binary_form(f2);
    info!("Parity check matrix in binary form:{}", h2);
    let (g, info_set) = Goppa::generator_matrix(&h2);
    info!("Generator matrix G:{}", g);
    let s = Mat::invertible_random(&mut rng, f2, k as usize);
    info!("Singular matrix S:{}", s);
    let p = Mat::permutation_random(&mut rng, f2, n as usize);
    info!("Permutation matrix P:{}", p);
    let sgp = &s * &g * &p;
    info!("Perturbed generator matrix ~G:{}", sgp);
    // let pk = (gg, t);
    // let sk = (s, goppa, p);
    // (pk, sk)
    let pk = PublicKey { sgp, t };
    println!("{:?}", info_set);
    let sk = SecretKey { s, goppa, p , info_set };
    // let sk = SecretKey { s, goppa, p , info_set: Vec::new() };
    (pk, sk)
}

impl<'a> PublicKey<'a> {
    pub fn encrypt(&self, m: &RowVec<'a, F2>) -> RowVec<'a, F2> {
        let mut rng = rand::thread_rng();
        // println!("sgp({}, {})", self.sgp.rows(), self.sgp.cols());
        let c = m * &self.sgp;
        let z = RowVec::random_with_weight(&mut rng, m.field(), self.sgp.cols(), self.t as usize);
        c + z
    }

    pub fn save_public_key(&self, file_name: &str) {
        // println!("{}\n", self.sgp);
        let mut f = File::create(file_name).expect("Unable to create file");
        f.write_all(self.sgp.to_hex_string().as_bytes());
        f.write_all(format!("\n{:02x}\n", self.t).as_bytes());
    }

    pub fn load_public_key(file_name: &str, f2: &'a F2) -> Self {
        let f = File::open(file_name).expect("Unable to open file");
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)
            .expect("Unable to read the first line");
        line.pop(); // Remove terminating newline
        let sgp = Mat::from_hex_string(&line, f2);
        // println!("{}\n", sgp);
        line.clear();
        f.read_line(&mut line)
            .expect("Unable to read the second line");
        line.pop(); // Remove terminating newline
        let t = u32::from_str_radix(&line, 16).unwrap();
        PublicKey { sgp, t }
    }
}

impl<'a, 'b> SecretKey<'a, 'b> {
    pub fn decrypt(
        &self,
        c: &RowVec<'a, F2>, // ciphertext
    ) -> RowVec<'a, F2> {
        let f2 = c.field();
        let m_s_g_z = c * self.p.transpose(); // msg stands for m * s * g
        let mut m_s_g = self.goppa.decode(&m_s_g_z);

        // let h = self.goppa.parity_check_matrix().binary_form(f2);
        // let (g, p) = Goppa::generator_matrix(&h);
        // let m_s_g = m_s_g * p;
        // let m_s = RowVec::new(f2, m_s_g.data()[0..self.s.rows()].to_vec());
        
        // let k = self.info_set.len();
        // let mut v = Vec::with_capacity(k);
        // for i in 0..k {
        //     v.push(m_s_g[self.info_set[i]]);
        // }
        // let m_s = RowVec::new(f2, v);

        m_s_g.permute_cols(&self.info_set);
        let m_s = RowVec::new(f2, m_s_g.data()[0..self.s.rows()].to_vec());
        
        let m = m_s * self.s.inverse().unwrap();
        m
        // self.goppa.decode(&(cpt * self.p.inverse().unwrap())) * self.s.inverse().unwrap()
    }

    pub fn save_secret_key(&self, file_name: &str) {
        let mut f = File::create(file_name).expect("Unable to create file");
        f.write_all((self.s.to_hex_string() + "\n").as_bytes());
        f.write_all((self.goppa.to_hex_string() + "\n").as_bytes());
        f.write_all((self.p.to_hex_string() + "\n").as_bytes());
    }

    pub fn load_finite_field(file_name: &str) -> F2m {
        // TODO: return result and use ? check all expect and unwrap calls
        let f = File::open(file_name).expect("Unable to open file");
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)
            .expect("Unable to read the first line");
        line.clear();
        f.read_line(&mut line)
            .expect("Unable to read the second line");
        let order = u32::from_str_radix(line.split('#').next().unwrap(), 16).unwrap();
        let f = Field::generate(order);
        f
    }

    pub fn load_secret_key(file_name: &str, f2: &'a F2, f2m: &'b F2m) -> Self {
        let f = File::open(file_name).expect("Unable to open file");
        let f = BufReader::new(f);
        let mut lines = f.lines();
        let s = Mat::from_hex_string(&lines.next().unwrap().unwrap(), f2); // TODO: double unwrap ??
        let goppa = Goppa::from_hex_string(&lines.next().unwrap().unwrap(), f2m); // TODO: double unwrap ??
        let p = Mat::from_hex_string(&lines.next().unwrap().unwrap(), f2); // TODO: double unwrap ??
        SecretKey { s, goppa, p , info_set: Vec::new() }
    }
}
