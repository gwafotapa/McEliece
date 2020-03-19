use rand::Rng;

use mceliece::{
    crypto::*,
    finite_field::{F2m, Field, F2},
    goppa::Goppa,
    matrix::RowVec,
};

#[test]
fn crypto_keygen_128() {
    let f2 = &F2::generate(2);
    let m = 7;
    let f128 = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 10;
    let (pk, sk) = keygen(f2, f128, n, t);
    let h = sk.goppa.parity_check_matrix().binary_form(&f2);
    let (g, _) = Goppa::generator_matrix(&h);
    assert!(sk.s.is_invertible());
    // assert!(sk.p.is_permutation());
    assert_eq!(g.rank(), g.rows());
    assert_eq!(pk.sgp, sk.s * g * sk.p);
}

#[test]
fn crypto_pk_save_load() {
    let f2 = &F2::generate(2);
    let m = 8;
    let f256 = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 19;
    let (pk_save, _) = keygen(f2, f256, n, t);
    let file_name = "pk_save_load_test.mce";
    pk_save.save_public_key(file_name);
    let pk_load = PublicKey::load_public_key(file_name, f2).unwrap();
    assert_eq!(pk_save, pk_load);
}

#[test]
fn crypto_sk_save_load() {
    let f2 = &F2::generate(2);
    let m = 8;
    let f256_save = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 19;
    let (_, sk_save) = keygen(f2, f256_save, n, t);
    let file_name = "sk_save_load_test.mce";
    sk_save.save_secret_key(file_name);
    let f256_load = &SecretKey::load_finite_field(file_name).unwrap();
    let sk_load = SecretKey::load_secret_key(file_name, f2, f256_load).unwrap();
    assert!(f256_save == f256_load);
    assert!(sk_save == sk_load);

    // assert!(sk_save.s == sk_load.s);
    // assert!(sk_save.goppa == sk_load.goppa);
    // assert!(sk_save.info_set == sk_load.info_set);
    // assert!(sk_save.p == sk_load.p);
}

#[test]
fn crypto_encrypt_decrypt_null_message() {
    let file_pk = "public_key_test.mce";
    let file_sk = "secret_key_test.mce";
    let file_msg = "message_test.mce";
    let file_cpt = "ciphertext_test.mce";
    let file_dcd_msg = "decoded_message_test.mce";

    let f2 = &F2::generate(2);
    let m = 7;
    let f128 = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 10;
    let (pk, sk) = keygen(f2, f128, n, t);
    pk.save_public_key(file_pk);
    sk.save_secret_key(file_sk);

    // let k = n - m * t;
    let k = pk.sgp.rows();
    assert!(k as u32 >= n - m * t);
    let msg = RowVec::zero(f2, k);
    msg.save_vector(file_msg);

    let pk = PublicKey::load_public_key(file_pk, f2).unwrap();
    let msg = RowVec::load_vector(file_msg, f2).unwrap();
    let cpt = pk.encrypt(&msg);
    cpt.save_vector(file_cpt);
    assert_eq!(cpt.weight() as u32, t);

    let f2m = &SecretKey::load_finite_field(file_sk).unwrap();
    let sk = SecretKey::load_secret_key(file_sk, f2, f2m).unwrap();
    let cpt = RowVec::load_vector(file_cpt, f2).unwrap();
    let dcd_msg = sk.decrypt(&cpt);
    dcd_msg.save_vector(file_dcd_msg);
    assert_eq!(dcd_msg, msg);
}

#[test]
fn crypto_encrypt_decrypt_random_message() {
    let file_pk = "public_key_test.mce";
    let file_sk = "secret_key_test.mce";
    let file_msg = "message_test.mce";
    let file_cpt = "ciphertext_test.mce";
    let file_dcd_msg = "decoded_message_test.mce";

    let f2 = &F2::generate(2);
    let m = 7;
    let f128 = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 10;
    let (pk, sk) = keygen(f2, f128, n, t);
    pk.save_public_key(file_pk);
    sk.save_secret_key(file_sk);

    let k = pk.sgp.rows();
    assert!(k as u32 >= n - m * t);
    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, k);
    msg.save_vector(file_msg);

    let pk = PublicKey::load_public_key(file_pk, f2).unwrap();
    // let msg = RowVec::load_vector(file_msg, f2);
    let cpt = pk.encrypt(&msg);
    cpt.save_vector(file_cpt);

    let f2m = &SecretKey::load_finite_field(file_sk).unwrap();
    let sk = SecretKey::load_secret_key(file_sk, f2, f2m).unwrap();
    // let cpt = RowVec::load_vector(file_cpt, f2);
    let dcd_msg = sk.decrypt(&cpt);
    dcd_msg.save_vector(file_dcd_msg);
    assert_eq!(dcd_msg, msg);
}

#[test]
fn crypto_encrypt_decrypt_random_message_without_error() {
    let file_pk = "public_key_test.mce";
    let file_sk = "secret_key_test.mce";
    let file_msg = "message_test.mce";
    let file_cpt = "ciphertext_test.mce";
    let file_dcd_msg = "decoded_message_test.mce";

    let f2 = &F2::generate(2);
    let m = 7;
    let f128 = &F2m::generate(1 << m);
    let n = 1 << m;
    let t = 10;
    let (pk, sk) = keygen(f2, f128, n, t);
    pk.save_public_key(file_pk);
    sk.save_secret_key(file_sk);

    let k = pk.sgp.rows();
    assert!(k as u32 >= n - m * t);
    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, k);
    msg.save_vector(file_msg);

    // let pk = PublicKey::load_public_key(file_pk, f2);
    let cpt = &msg * &pk.sgp;
    cpt.save_vector(file_cpt);
    // let h = sk.goppa.parity_check_matrix().binary_form(f2);
    // assert_eq!(&h * (&cpt * sk.p.inverse().unwrap()).transpose(), Mat::zero(f2, h.rows(), cpt.rows()));

    // let f2m = &SecretKey::load_finite_field(file_sk);
    // let sk = SecretKey::load_secret_key(file_sk, f2, f2m);
    let dcd_msg = sk.decrypt(&cpt);
    dcd_msg.save_vector(file_dcd_msg);
    assert_eq!(dcd_msg, msg);
}

// TODO: (m, n, t) = (3, 7, 1) fails
#[test]
fn crypto_encrypt_decrypt_random_message_L_not_full() {
    let file_pk = "public_key_test.mce";
    let file_sk = "secret_key_test.mce";
    let file_msg = "message_test.mce";
    let file_cpt = "ciphertext_test.mce";
    let file_dcd_msg = "decoded_message_test.mce";

    let mut rng = rand::thread_rng();
    let f2 = &F2::generate(2);
    let m = rng.gen_range(3, 10);
    let t = rng.gen_range(1, ((1 << m) - 1) / m);
    let f2m = &F2m::generate(1 << m);
    let n = rng.gen_range(m * t + 1, 1 << m);
    println!("m: {}\nn: {}\nt: {}", m, n, t);
    let (pk, sk) = keygen(f2, f2m, n, t);
    // pk.save_public_key(file_pk);
    // sk.save_secret_key(file_sk);

    let k = pk.sgp.rows();
    assert!(k as u32 >= n - m * t);
    let msg = RowVec::random(&mut rng, f2, k);
    // msg.save_vector(file_msg);

    // let pk = PublicKey::load_public_key(file_pk, f2);
    // let msg = RowVec::load_vector(file_msg, f2);
    let cpt = pk.encrypt(&msg);
    // cpt.save_vector(file_cpt);

    // let f2m = &SecretKey::load_finite_field(file_sk);
    // let sk = SecretKey::load_secret_key(file_sk, f2, f2m);
    // let cpt = RowVec::load_vector(file_cpt, f2);
    let dcd_msg = sk.decrypt(&cpt);
    // dcd_msg.save_vector(file_dcd_msg);
    assert_eq!(dcd_msg, msg);
}
