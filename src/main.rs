// use log::info;

use std::{
    env,
    // fs::File,
    // io::{BufRead, BufReader, Read, Write},
};

use crypto::{PublicKey, SecretKey};
use finite_field::{F2m, F2, Field};
// use goppa::Goppa;
use matrix::RowVec;
// use polynomial::Poly;

mod crypto;
mod finite_field;
mod goppa;
mod matrix;
mod polynomial;

const GOPPA_M: u32 = 5; // The underlying field of the code is F2^GOPPA_M
const GOPPA_N: u32 = 1 << GOPPA_M; // Code length
const GOPPA_T: u32 = 4; // Code correction capacity i.e. degree of the goppa polynomial



fn main() {
    // TODO: deal with all the unwrap() and expect() calls
    env_logger::init();
    
    if GOPPA_N > 255 {
        panic!("n should be less than 256");
    }
    if GOPPA_N % 8 != 0 {
        panic!("n should be a multiple of 8");
    }
    if GOPPA_N < GOPPA_M * GOPPA_T {
        panic!("m * t should be less than n");
    }
    let args: Vec<String> = env::args().collect();
    if args.len() == 1 {
        panic!("Command expected.");
    }
    let f2 = &F2 {};
    match args[1].as_str() {
        "keygen" => {
            let f2m = &F2m::generate(1 << GOPPA_M);
            let (pk, sk) = crypto::keygen(f2, f2m, GOPPA_N, GOPPA_T);
            pk.save_public_key("public_key.mce");
            sk.save_secret_key("secret_key.mce");
            ()
        }
        "encrypt" => {
            let pk = PublicKey::load_public_key("public_key.mce", f2);
            let m = RowVec::load_vector("message.mce", f2);
            let c = pk.encrypt(&m);
            c.save_vector("ciphertext.mce");
            ()
        }
        "decrypt" => {
            let f2m = &SecretKey::load_finite_field("secret_key.mce"); // &F2m::generate(1 << GOPPA_M);
            let sk = SecretKey::load_secret_key("secret_key.mce", f2, f2m);
            let c = RowVec::load_vector("ciphertext.mce", f2);
            let m = sk.decrypt(&c);
            m.save_vector("decoded.mce"); // TODO: pick convention for arguments order
            ()
        }
        _ => panic!("Unexpected command. Valid commands are 'keygen', 'encrypt' and 'decrypt'."),
    }
}
