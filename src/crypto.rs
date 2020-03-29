//! The McEliece cryptosystem

use log::debug;

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

use crate::finite_field::{F2m, F2};
use crate::goppa::Goppa;
use crate::matrix::{Mat, Perm, RowVec};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// Public key of the McEliece cryptosystem
///
/// A public key is a couple (SGP, t) where:
/// - SGP is the product of three matrices S, G and P
///   (see struct [SecretKey](struct.SecretKey.html) for more)
/// - t is the correction capacity of the Goppa code (and also the degree of the Goppa polynomial)
#[derive(Debug, Eq, PartialEq)]
pub struct PublicKey<'a> {
    sgp: Mat<'a, F2>,
    t: u32,
}

/// Secret key of the McEliece cryptosystem
///
/// A secret key has three components:
/// - a random invertible matrix S
/// - a random binary irreducible Goppa code of generator matrix G
/// - a random permutation matrix P
///
/// Instead of storing P, we save the corresponding element of the symmetric group.  
/// info_set is the information set of G.  
/// See <https://en.wikipedia.org/wiki/McEliece_cryptosystem>
/// for a brief description of the McEliece cryptosystem.  
/// For more details, see Engelbert, D., Overbeck, R., & Schmidt, A. (2007),
/// A summary of McEliece-type cryptosystems and their security.
/// Journal of Mathematical Cryptology JMC, 1(2), 151-199.
#[derive(Eq, PartialEq)]
pub struct SecretKey<'a, 'b> {
    s: Mat<'a, F2>,
    goppa: Goppa<'b, F2m>,
    info_set: Vec<usize>,
    p: Perm,
}

pub fn keygen<'a, 'b>(
    f2: &'a F2,
    f2m: &'b F2m,
    n: u32,
    t: u32,
) -> (PublicKey<'a>, SecretKey<'a, 'b>) {
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f2m, n as usize, t as usize);
    debug!("{}", goppa);
    
    let xyz = goppa.parity_check_xyz();
    let (g, info_set) = Goppa::generator_from_xyz(&xyz, f2);
    debug!("Generator matrix G:{}", g);
    debug!("Information set of generator matrix G:\n{:?}\n", info_set);
    
    let k = g.rows();
    let s = Mat::invertible_random(&mut rng, f2, k as usize);
    debug!("Code dimension k = {}", k);
    debug!("Singular matrix S:{}", s);
    
    let p = Perm::random(&mut rng, n as usize);
    debug!("Permutation P:\n{:?}\n", p);
    
    let sgp = &s * &g * &p;
    debug!("Perturbed generator matrix ~G:{}", sgp);
    
    let pk = PublicKey { sgp, t };
    let sk = SecretKey {
        s,
        goppa,
        info_set,
        p,
    };
    (pk, sk)
}

impl<'a> PublicKey<'a> {
    pub fn sgp(&self) -> &Mat<'a, F2> {
        &self.sgp
    }

    pub fn t(&self) -> u32 {
        self.t
    }
    
    pub fn encrypt(&self, m: &RowVec<'a, F2>) -> RowVec<'a, F2> {
        let mut rng = rand::thread_rng();
        let c = Goppa::<F2m>::g_encode(&self.sgp, m);
        debug!("Encoded plaintext:{}", c);
        
        let z = RowVec::random_with_weight(&mut rng, m.field(), self.sgp.cols(), self.t as usize);
        debug!("Error vector:{}", z);
        
        c + z
    }

    pub fn write(&self, file_name: &str) -> Result<()> {
        let mut f = File::create(file_name)?;
        let s = self.sgp.to_hex_string() + format!("\n{:02x}\n", self.t).as_str();
        f.write_all(s.as_bytes())?;
        Ok(())
    }

    pub fn read_public_key(file_name: &str, f2: &'a F2) -> Result<Self> {
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

    pub fn read_code_dimension(file_name: &str) -> Result<u32> {
        let f = File::open(file_name)?;
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)?;
        let k = u32::from_str_radix(line.split('#').next().ok_or("Cannot read k")?, 16)?;
        Ok(k)
    }
    
}

impl<'a, 'b> SecretKey<'a, 'b> {
    pub fn s(&self) -> &Mat<'a, F2> {
        &self.s
    }

    pub fn goppa(&self) -> &Goppa<'b, F2m> {
        &self.goppa
    }

    pub fn info_set(&self) -> &Vec<usize> {
        &self.info_set
    }

    pub fn p(&self) -> &Perm {
        &self.p
    }

    pub fn decrypt(&self, c: &RowVec<'a, F2>) -> RowVec<'a, F2> {
        let m_s_g_z = c * self.p.inverse();
        let m_s_g = self.goppa.decode(&m_s_g_z);
        debug!("Decoded codeword mSG:{}", m_s_g);
        
        let m_s = m_s_g.extract_cols(&self.info_set);
        debug!(
            "Use information set {:?} to extract mS:{}",
            self.info_set, m_s
        );
        
        let m = m_s * self.s.inverse().unwrap();
        m
    }

    pub fn write(&self, file_name: &str) -> Result<()> {
        let f = File::create(file_name)?;
        let mut f = BufWriter::new(f);
        f.write_all((self.s.to_hex_string() + "\n").as_bytes())?;
        f.write_all((self.goppa.to_hex_string() + "\n").as_bytes())?;

        for i in 0..self.info_set.len() - 1 {
            f.write_all(format!("{:x} ", self.info_set[i]).as_bytes())?;
        }
        f.write_all(format!("{:x}\n", self.info_set[self.info_set.len() - 1]).as_bytes())?;

        for i in 0..self.p.len() - 1 {
            f.write_all(format!("{:x} ", self.p[i]).as_bytes())?;
        }
        f.write_all(format!("{:x}\n", self.p[self.p.len() - 1]).as_bytes())?;
        Ok(())
    }

    pub fn read_finite_field(file_name: &str) -> Result<F2m> {
        let f = File::open(file_name)?;
        let mut f = BufReader::new(f);
        let mut line = String::new();
        f.read_line(&mut line)?;
        line.clear();
        f.read_line(&mut line)?;
        let order = u32::from_str_radix(&line[..line.find('#').ok_or("Missing hashtag")?], 16)?;
        Ok(F2m::generate(order))
    }

    pub fn read_secret_key(file_name: &str, f2: &'a F2, f2m: &'b F2m) -> Result<Self> {
        let f = File::open(file_name)?;
        let f = BufReader::new(f);
        let mut lines_iter = f.lines().map(|l| l.unwrap());

        let line = lines_iter.next().ok_or("Cannot read line 1")?;
        let s = Mat::from_hex_string(&line, f2)?;
        debug!("Read matrix s:{}", s);

        let line = lines_iter.next().ok_or("Cannot read line 2")?;
        let goppa = Goppa::from_hex_string(&line, f2m)?;
        debug!("Read goppa code:{}", goppa);

        let line = lines_iter.next().ok_or("Cannot read line 3")?;
        let v: Vec<&str> = line.split(' ').collect();
        let mut info_set = Vec::new();
        for i in 0..v.len() {
            info_set.push(usize::from_str_radix(v[i], 16)?);
        }
        debug!("Read information set:{:?}", info_set);

        let line = lines_iter.next().ok_or("Cannot read line 4")?;
        let v: Vec<&str> = line.split(' ').collect();
        let mut p = Vec::new();
        for i in 0..v.len() {
            p.push(usize::from_str_radix(v[i], 16)?);
        }
        let p = Perm::new(p);

        Ok(SecretKey {
            s,
            goppa,
            info_set,
            p,
        })
    }
}
