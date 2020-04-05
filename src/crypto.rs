//! The McEliece cryptosystem

use log::debug;

use std::{
    convert::TryInto,
    error::Error,
    fs::File,
    io::{BufWriter, Read, Write},
    rc::Rc,
};

use crate::finite_field::{F2m, Field, FieldTrait, FiniteField, F2};
use crate::goppa::Goppa;
use crate::matrix::{Mat, Perm, RowVec};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// Public key of the McEliece cryptosystem
///
/// A public key is a couple (SGP, t) where:
/// - SGP is the product of three matrices S, G and P
///   (see struct [`SecretKey`] for more)
/// - t is the correction capacity of the Goppa code (and also the degree of the Goppa polynomial)
///
/// [`SecretKey`]: struct.SecretKey.html
#[derive(Debug, Eq, PartialEq)]
pub struct PublicKey {
    sgp: Mat<F2>,
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
pub struct SecretKey {
    s: Mat<F2>,
    goppa: Goppa<F2m>,
    info_set: Vec<usize>,
    p: Perm,
}

pub fn keygen(n: usize, t: usize) -> (PublicKey, SecretKey) {
    let q = if t == 1 && n.is_power_of_two() {
        2 * n
    } else {
        n.next_power_of_two()
    };
    let f2 = &Rc::new(F2 {});
    let f2m = &Rc::new(F2m::generate(q));
    let goppa = Goppa::random(Field::Some(f2m), n, t);
    debug!("{}", goppa);

    let xyz = goppa.parity_check_xyz();
    let (g, info_set) = Goppa::generator_from_xyz(&xyz, Field::Some(f2));
    debug!("Generator matrix G:{}", g);
    debug!("Information set of generator matrix G:\n{:?}\n", info_set);

    let k = g.rows();
    let s = Mat::invertible_random(Field::Some(f2), k);
    debug!("Code dimension k = {}", k);
    debug!("Singular matrix S:{}", s);

    let p = Perm::random(n);
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

impl PublicKey {
    pub fn sgp(&self) -> &Mat<F2> {
        &self.sgp
    }

    pub fn t(&self) -> usize {
        self.t
    }

    pub fn encrypt(&self, m: &RowVec<F2>) -> RowVec<F2> {
        let c = Goppa::<F2m>::g_encode(&self.sgp, m);
        debug!("Encoded plaintext:{}", c);

        let z = RowVec::random_with_weight(Field::Some(m.field()), self.sgp.cols(), self.t);
        debug!("Error vector:{}", z);

        c + z
    }

    /// Saves public key on disk
    ///
    /// The output file layout is:
    /// - bytes 0-3: number k of rows of matrix sgp
    /// - bytes 4-7: number n of columns of matrix sgp
    /// - bytes 8-x: coefficients of matrix sgp (eight per byte)
    /// - bytes x-x+4: correction capacity t
    pub fn write(&self, file_name: &str) -> Result<()> {
        let mut f = File::create(file_name)?;
        f.write_all(&self.sgp.to_bytes())?;
        f.write_all(&(self.t as u32).to_be_bytes())?;
        Ok(())
    }

    pub fn read_public_key(file_name: &str) -> Result<Self> {
        let mut f = File::open(file_name)?;
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;
        let (i, sgp) = Mat::from_bytes(&vec)?;
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

impl SecretKey {
    pub fn s(&self) -> &Mat<F2> {
        &self.s
    }

    pub fn goppa(&self) -> &Goppa<F2m> {
        &self.goppa
    }

    pub fn info_set(&self) -> &Vec<usize> {
        &self.info_set
    }

    pub fn p(&self) -> &Perm {
        &self.p
    }

    pub fn decrypt(&self, c: &RowVec<F2>) -> RowVec<F2> {
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

    /// Saves secret key on disk
    ///
    /// The output file layout is:
    /// - bytes 0-3: finite field order q
    /// - bytes 4-7: number of rows of matrix s
    /// - bytes 8-11: number of columns of matrix s
    /// - bytes 11-a: coefficients of matrix s (eight per byte)
    /// - bytes a-a+4: finite field order q
    /// - bytes a+4-a+8: correction capacity t
    /// - bytes a+8-b: coefficients of Goppa polynomial (four bytes per coefficient)
    /// - bytes b-c: Goppa set L (see [`Goppa`]::[`to_bytes()`] for more information)
    /// - bytes c-d: information set (four bytes per column index)
    /// - bytes d-d+4: length n of the code
    /// - bytes d+4-e: permutation P (four bytes per image)
    ///
    /// [`Goppa`]: ../goppa/struct.Goppa.html
    /// [`to_bytes()`]: ../goppa/struct.Goppa.html#method.to_bytes
    pub fn write(&self, file_name: &str) -> Result<()> {
        let f = File::create(file_name)?;
        let mut f = BufWriter::new(f);
        f.write_all(&(self.goppa.field().order() as u32).to_be_bytes())?;
        f.write_all(&self.s.to_bytes())?;
        f.write_all(&self.goppa.to_bytes())?;

        for i in 0..self.info_set.len() {
            f.write_all(&(self.info_set[i] as u32).to_be_bytes())?;
        }

        f.write_all(&(self.p.len() as u32).to_be_bytes())?;
        for i in 0..self.p.len() {
            f.write_all(&(self.p[i] as u32).to_be_bytes())?;
        }
        Ok(())
    }

    // // TODO: Is it possible to eliminate this function ?
    // pub fn read_finite_field(file_name: &str) -> Result<F2m> {
    //     let mut f = File::open(file_name)?;
    //     // let mut vec = Vec::new();
    //     // f.read_to_end(&mut vec)?;
    //     let mut buf = [0; 4];
    //     f.read_exact(&mut buf)?;
    //     let order = u32::from_be_bytes(buf) as usize;
    //     // let k = u32::from_be_bytes(vec[0..4].try_into()?) as usize;
    //     // let i = 4 + 4 + div_ceil(k * k, 8);
    //     // let order = u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize;
    //     Ok(F2m::generate(order))
    // }

    pub fn read_secret_key(file_name: &str) -> Result<Self> {
        let mut f = File::open(file_name)?;
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;
        let mut i = 4;

        let (bytes, s) = Mat::from_bytes(&vec[i..])?;
        i += bytes;
        debug!("Read matrix s:{}", s);

        let (bytes, goppa) = Goppa::from_bytes(&vec[i..])?;
        i += bytes;
        debug!("Read Goppa code:\n{}", goppa);

        let k = s.rows();
        let mut info_set = Vec::new();
        for _j in 0..k {
            info_set.push(u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize);
            i += 4;
        }
        debug!("Read information set:\n{:?}", info_set);

        let p_len = u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize;
        i += 4;
        let mut p = Vec::with_capacity(p_len);
        for _j in 0..p_len {
            p.push(u32::from_be_bytes(vec[i..i + 4].try_into()?) as usize);
            i += 4;
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
