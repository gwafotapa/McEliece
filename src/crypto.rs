use log::info;

use crate::finite_field::{F2m, F2};
use crate::goppa::Goppa;
use crate::matrix::{Mat, RowVec};
use crate::polynomial::Poly;

pub struct PublicKey<'a> {
    pub sgp: Mat<'a, F2>, // disguised generator matrix S * G * P
    pub t: u32,           // degree of the goppa polynomial
}

pub struct SecretKey<'a, 'b> {
    pub s: Mat<'a, F2>,        // singular matrix S
    pub goppa: Goppa<'b, F2m>, // goppa code of generator matrix G
    pub p: Mat<'a, F2>,        // permutation matrix P
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
    let g = Goppa::generator_matrix(&h2);
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
    let sk = SecretKey { s, goppa, p };
    (pk, sk)
}

// TODO: Should I implement this as a method of PublicKey and decrypt as method of SecretKey ?
pub fn encrypt<'a>(pk: &PublicKey<'a>, m: &RowVec<'a, F2>) -> RowVec<'a, F2> {
    let mut rng = rand::thread_rng();
    // let cdw = msg * &pk.0;
    println!("sgp({}, {})", pk.sgp.rows(), pk.sgp.cols());
    let c = m * &pk.sgp;
    // let z = RowVec::random_with_weight(&mut rng, f2, msg.cols(), pk.1 as usize);
    let z = RowVec::random_with_weight(&mut rng, m.field(), pk.sgp.cols(), pk.t as usize);
    c + z
}

pub fn decrypt<'a, 'b>(
    sk: &SecretKey<'a, 'b>,
    c: &RowVec<'a, F2>, // ciphertext
) -> RowVec<'a, F2> {
    // let s = &sk.0;
    // let goppa = &sk.1;
    // let p = &sk.2;
    // let m_s_g = cpt * p.inverse().unwrap();
    // let m_s = goppa.decode(&m_s_g);
    // let m = m_s * s.inverse().unwrap();
    // m

    // sk.1.decode(&(cpt * sk.2.inverse().unwrap())) * sk.0.inverse().unwrap()

    let msg = c * sk.p.inverse().unwrap(); // msg stands for m * s * g
    let ms = sk.goppa.decode(&msg);
    println!("ms({}, {}), s.inv({}, {})", ms.rows(), ms.cols(), sk.s.rows(), sk.s.cols());
    let ms = RowVec::new(c.field(), ms.data()[0..sk.s.rows()].to_vec());
    let m = ms * sk.s.inverse().unwrap();
    m
    // sk.goppa.decode(&(cpt * sk.p.inverse().unwrap())) * sk.s.inverse().unwrap()
}
