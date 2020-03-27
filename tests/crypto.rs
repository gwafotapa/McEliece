use log::warn;
// use rand::Rng;

use mceliece::{
    crypto::*,
    finite_field::{F2m, F2},
    // goppa::Goppa,
    matrix::{RowVec},
};

pub mod common;

const REPEAT: u32 = 1;

// fn setup() -> (u32, u32, u32) {
//     common::log_setup();
//     let mut rng = rand::thread_rng();
//     let n = match GOPPA_N {
//         0 => rng.gen_range(GOPPA_N_MIN, GOPPA_N_MAX),
//         value => value,
//     };
//     let q = n.next_power_of_two();
//     let m = q.trailing_zeros();
//     let t = match GOPPA_T {
//         0 => {
//             let t_min = if n == q { 2 } else { 1 };
//             let t_max = matrix::div_ceil(n, m) + if n.is_power_of_two() { 1 } else { 0 };
//             rng.gen_range(t_min, t_max)
//         }
//         value => value,
//     };
//     (m, n, t)
// }

#[test]
fn crypto_keygen() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, sk) = keygen(f2, f2m, n, t);
    let g = sk.goppa().generator_matrix(f2);
    assert!(sk.s().is_invertible());
    assert_eq!(g.rank(), g.rows());
    assert!(sk.p().is_permutation());
    assert_eq!(*pk.sgp(), sk.s() * g * sk.p());
}

#[test]
fn crypto_pk_write_read() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, _) = keygen(f2, f2m, n, t);
    let file_name = "pk_write_read_test.mce";
    pk.write(file_name).unwrap();
    let pk_read = PublicKey::read_public_key(file_name, f2).unwrap();
    assert_eq!(pk, pk_read);
}

#[test]
fn crypto_sk_write_read() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (_, sk) = keygen(f2, f2m, n, t);
    let file_name = "sk_write_read_test.mce";
    sk.write(file_name).unwrap();
    let f2m_read = &SecretKey::read_finite_field(file_name).unwrap();
    let sk_read = SecretKey::read_secret_key(file_name, f2, f2m_read).unwrap();
    assert!(f2m == f2m_read);
    assert!(sk == sk_read);
}

#[test]
fn crypto_decrypt_null_ciphertext() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, sk) = keygen(f2, f2m, n, t);
    let k = pk.sgp().rows();
    let msg = RowVec::zero(f2, k);
    let cpt = RowVec::zero(f2, n as usize);
    let dmsg = sk.decrypt(&cpt);
    assert_eq!(dmsg, msg);
}

#[test]
fn crypto_decrypt_codeword_ciphertext() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, sk) = keygen(f2, f2m, n, t);
    let k = pk.sgp().rows();
    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, k);
    let cpt = &msg * pk.sgp();
    let dmsg = sk.decrypt(&cpt);
    assert_eq!(dmsg, msg);
}

#[test]
fn crypto_encrypt_decrypt_null_message() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, sk) = keygen(f2, f2m, n, t);
    let k = pk.sgp().rows();
    let msg = RowVec::zero(f2, k);
    let cpt = pk.encrypt(&msg);
    assert_eq!(cpt.weight() as u32, t);
    let dmsg = sk.decrypt(&cpt);
    assert_eq!(dmsg, msg);
}

#[test]
fn crypto_encrypt_decrypt_random_message() {
    let (m, n, t) = common::goppa_setup();
    warn!("m={}, n={}, t={}", m, n, t);
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let (pk, sk) = keygen(f2, f2m, n, t);
    let k = pk.sgp().rows();
    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, k);
    let cpt = pk.encrypt(&msg);
    let dmsg = sk.decrypt(&cpt);
    assert_eq!(dmsg, msg);
}

#[test]
fn crypto_repeat_encrypt_decrypt() {
    let mut r = 0;
    while r < REPEAT {
        crypto_encrypt_decrypt_random_message();
        r += 1;
    }
}
