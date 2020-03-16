// use log::info;

use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, Read, Write},
};

use crypto::{PublicKey, SecretKey};
use finite_field::{F2m, F2, Field};
use goppa::Goppa;
use matrix::{Mat, RowVec};
// use polynomial::Poly;

mod crypto;
mod finite_field;
mod goppa;
mod matrix;
mod polynomial;

const GOPPA_M: u32 = 5; // The underlying field of the code is F2^GOPPA_M
const GOPPA_N: u32 = 1 << GOPPA_M; // Code length
const GOPPA_T: u32 = 4; // Code correction capacity i.e. degree of the goppa polynomial

fn save_vector(file_name: &str, vec: RowVec<'_, F2>) {
    let mut f = File::create(file_name).expect("Unable to create file");
    let mut s = format!("{:x}#", vec.cols());
    let mut byte: u8 = 0;
    let mut cnt_mod_8 = 7;
    for i in 0..vec.cols() {
        byte |= (vec[i] << cnt_mod_8) as u8;
        if cnt_mod_8 == 0 {
            s.push_str(format!("{:02x}", byte).as_str());
            cnt_mod_8 = 7;
            byte = 0;
        } else {
            cnt_mod_8 -= 1;
        }
    }
    f.write_all(s.as_bytes());
}

fn load_vector<'a>(file_name: &str, f2: &'a F2) -> RowVec<'a, F2> {
    let mut f = File::open(file_name).expect("Unable to open file");
    let mut buf = String::new();
    f.read_to_string(&mut buf).expect("Unable to read data");
    let v: Vec<&str> = buf.split('#').collect();
    let cols = u32::from_str_radix(v[0], 16).unwrap() as usize;
    let mut vec = RowVec::zero(f2, cols);
    let mut cnt_mod_8 = 0;
    let mut iter = v[1].as_bytes().iter();
    let mut byte = iter.next().unwrap();
    for i in 0..cols {
        vec[i] = (*byte as u32 >> cnt_mod_8) & 1_u32;
        if cnt_mod_8 == 7 {
            byte = iter.next().unwrap();
            cnt_mod_8 = 0;
        } else {
            cnt_mod_8 += 1;
        }
    }
    vec
}

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
    match args[1].as_str() {
        "keygen" => {
            let f2 = &F2 {};
            let f2m = &F2m::generate(1 << GOPPA_M);
            let (pk, sk) = crypto::keygen(f2, f2m, GOPPA_M, GOPPA_N, GOPPA_T);
            save_public_key(&pk, "public_key.mce");
            save_secret_key(&sk, "secret_key.mce");
            ()
        }
        "encrypt" => {
            let f2 = &F2 {};
            let pk = load_public_key("public_key.mce", f2);
            let m = load_vector("message.mce", f2);
            let c = crypto::encrypt(&pk, &m);
            save_vector("ciphertext.mce", c);
            ()
        }
        "decrypt" => {
            let f2 = &F2 {};
            let f2m = &load_finite_field("secret_key.mce"); // &F2m::generate(1 << GOPPA_M);
            let sk = load_secret_key("secret_key.mce", f2, f2m);
            let c = load_vector("ciphertext.mce", f2);
            let m = crypto::decrypt(&sk, &c);
            save_vector("decoded.mce", m); // TODO: pick convention for arguments order
            ()
        }
        _ => panic!("Unexpected command. Valid commands are 'keygen', 'encrypt' and 'decrypt'."),
    }
}
