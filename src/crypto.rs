//! The McEliece cryptosystem

use log::debug;

use std::{
    convert::TryInto,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
};

use crate::finite_field::{F2m, FiniteField, F2};
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
    t: usize,
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
    n: usize,
    t: usize,
) -> (PublicKey<'a>, SecretKey<'a, 'b>) {
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f2m, n, t);
    debug!("{}", goppa);

    let xyz = goppa.parity_check_xyz();
    let (g, info_set) = Goppa::generator_from_xyz(&xyz, f2);
    debug!("Generator matrix G:{}", g);
    debug!("Information set of generator matrix G:\n{:?}\n", info_set);

    let k = g.rows();
    let s = Mat::invertible_random(&mut rng, f2, k);
    debug!("Code dimension k = {}", k);
    debug!("Singular matrix S:{}", s);

    let p = Perm::random(&mut rng, n);
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

    pub fn t(&self) -> usize {
        self.t
    }

    pub fn encrypt(&self, m: &RowVec<'a, F2>) -> RowVec<'a, F2> {
        let mut rng = rand::thread_rng();
        let c = Goppa::<F2m>::g_encode(&self.sgp, m);
        debug!("Encoded plaintext:{}", c);

        let z = RowVec::random_with_weight(&mut rng, m.field(), self.sgp.cols(), self.t);
        debug!("Error vector:{}", z);

        c + z
    }

    pub fn write(&self, file_name: &str) -> Result<()> {
        let mut f = File::create(file_name)?;
        f.write_all(&self.sgp.to_bytes());
        f.write_all(&self.t.to_be_bytes());
        Ok(())
    }

    pub fn read_public_key(file_name: &str, f2: &'a F2) -> Result<Self> {
        let mut f = File::open(file_name)?;
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;
        let sgp = Mat::from_bytes(&vec, f2)?;
        let i = 4 + 4 + div_ceil(sgp.rows() * sgp.cols(), 8);
        let t = u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize;
        Ok(PublicKey { sgp, t })
    }

    pub fn read_code_dimension(file_name: &str) -> Result<usize> {
        let mut f = File::open(file_name)?;
        let mut buf = [0; 4];
        f.read_exact(&mut buf)?;
        let k = u32::from_be_bytes(buf) as usize;
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
        f.write_all(&self.s.to_bytes())?;
        f.write_all(&self.goppa.to_bytes())?;

        // f.write_all(self.info_set.len().to_be_bytes());
        for i in 0..self.info_set.len() {
            f.write_all(&(self.info_set[i] as u32).to_be_bytes())?;
        }

        f.write_all(&(self.p.len() as u32).to_be_bytes());
        for i in 0..self.p.len() {
            f.write_all(&(self.p[i] as u32).to_be_bytes())?;
        }
        Ok(())
    }

    // TODO: Is it possible to eliminate this function ?
    pub fn read_finite_field(file_name: &str) -> Result<F2m> {
        let mut f = File::open(file_name)?;
        // let mut f = BufReader::new(f);
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;

        // let mut buf = [0; 4];
        // f.read_exact(&mut buf)?;
        // let k = u32::from_be_bytes(buf);
        let k = u32::from_be_bytes(vec[0..4].try_into()?) as usize;
        let i = 4 + 4 + div_ceil(k * k, 8);
        let order = u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize;

        // let mut line = String::new();
        // f.read_line(&mut line)?;
        // line.clear();
        // f.read_line(&mut line)?;
        // let order = u32::from_str_radix(&line[..line.find('#').ok_or("Missing hashtag")?], 16)?;
        Ok(F2m::generate(order))
    }

    // TODO: rewrite the use of index i for all read/write functions
    pub fn read_secret_key(file_name: &str, f2: &'a F2, f2m: &'b F2m) -> Result<Self> {
        let mut f = File::open(file_name)?;
        // let f = BufReader::new(f);
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;

        let s = Mat::from_bytes(&vec, f2)?;
        debug!("Read matrix s:{}", s);

        let k = s.rows();
        let mut i = 4 + 4 + div_ceil(k * k, 8);
        let goppa = Goppa::from_bytes(&vec[i..], f2m)?;
        debug!("Read Goppa code:\n{}", goppa);

        i += 4 + 4 + 4 * (goppa.poly().degree() + 1) + div_ceil(goppa.field().order(), 8); // TODO: Can I not use order() ???
        let mut info_set = Vec::new();
        for j in 0..k {
            info_set.push(u32::from_be_bytes(vec[i + 4 * j..i + 4 * j + 4].try_into()?) as usize);
        }
        debug!("Read information set:\n{:?}", info_set);

        i += 4 * k;
        let p_len = u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize;
        i += 4;
        let mut p = Vec::with_capacity(p_len);
        for j in 0..p_len {
            p.push(u32::from_be_bytes(vec[i + 4 * j..i + 4 * j + 4].try_into()?) as usize);
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

fn div_ceil(a: usize, b: usize) -> usize {
    a / b + if a % b == 0 { 0 } else { 1 }
}
