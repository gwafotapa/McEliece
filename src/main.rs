// use log::info;
use getopts::Options;
use std::{
    env,
    // fs::File,
    // io::{BufRead, BufReader, Read, Write},
    error::Error,
    result,
};

use crypto::{PublicKey, SecretKey};
use finite_field::{F2m, F2};
// use goppa::Goppa;
use matrix::RowVec;
// use polynomial::Poly;

mod crypto;
mod finite_field;
mod goppa;
mod matrix;
mod polynomial;

const GOPPA_M: u32 = 10; // F2^GOPPA_M is the underlying field of the code
const GOPPA_N: u32 = 1024; // Code length
const GOPPA_K: u32 = 524; // Code dimension
const GOPPA_T: u32 = 50; // Code correction capacity
const PLAINTEXT: &str = "plaintext.mce";
const CIPHERTEXT: &str = "ciphertext.mce";
const DECRYPTED: &str = "decrypted.mce";
const PUBLIC_KEY: &str = "public_key.mce";
const SECRET_KEY: &str = "secret_key.mce";

type Result<T> = result::Result<T, Box<dyn Error>>;

fn print_usage(program: &str, opts: Options) {
    let brief = format!(
        "Usage: {0} keygen [PK] [SK]\n\
                         {0} encrypt [PK] [PLAINTEXT] [CIPHERTEXT]\n\
                         {0} decrypt [SK] [CIPHERTEXT] [DECRYPTED]\n\
                         {0} plaintext [PLAINTEXT]\n\
                         Encrypts information using the McEliece cryptosystem.",
        program
    );
    print!("{}", opts.usage(&brief));
}

fn main() -> Result<()> {
    env_logger::init();

    if GOPPA_N < GOPPA_M * GOPPA_T {
        panic!("m * t should be less than n");
    }
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    opts.optflag("h", "help", "Display this help");
    opts.optflag("s", "silent", "Suppress display output");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!(f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&args[0], opts);
        return Ok(());
    }
    let silent = matches.opt_present("s");
    if args.len() == 1 {
        panic!("Command expected. Try --help for more information.");
    }
    let f2 = &F2 {};
    match args[1].as_str() {
        "keygen" => {
            let public_key = match args.get(2) {
                None => PUBLIC_KEY,
                Some(filename) => filename,
            };
            let secret_key = match args.get(3) {
                None => SECRET_KEY,
                Some(filename) => filename,
            };
            let f2m = &F2m::generate(1 << GOPPA_M);
            if !silent {
                print!("Generate keys...");
            }
            let (pk, sk) = crypto::keygen(f2, f2m, GOPPA_N, GOPPA_T);
            if !silent {
                print!("ok\nWrite public key to {}\n", public_key);
            }
            pk.write(public_key)?;
            if !silent {
                print!("Write secret key to {}\n", secret_key);
            }
            sk.write(secret_key)?;
            Ok(())
        }
        "encrypt" => {
            let public_key = match args.get(2) {
                None => PUBLIC_KEY,
                Some(filename) => filename,
            };
            let plaintext = match args.get(3) {
                None => PLAINTEXT,
                Some(filename) => filename,
            };
            let ciphertext = match args.get(4) {
                None => CIPHERTEXT,
                Some(filename) => filename,
            };
            if !silent {
                print!("Read public key from {}\n", public_key);
            }
            let pk = PublicKey::read_public_key(public_key, f2)?;
            if !silent {
                print!("Read plaintext from {}\n", plaintext);
            }
            let m = RowVec::read_vector(plaintext, f2)?;
            if !silent {
                print!("Encrypt plaintext...");
            }
            let c = pk.encrypt(&m);
            if !silent {
                print!("ok\nWrite ciphertext to {}\n", ciphertext);
            }
            c.write(ciphertext)?;
            Ok(())
        }
        "decrypt" => {
            let secret_key = match args.get(2) {
                None => SECRET_KEY,
                Some(filename) => filename,
            };
            let ciphertext = match args.get(3) {
                None => CIPHERTEXT,
                Some(filename) => filename,
            };
            let decrypted = match args.get(4) {
                None => DECRYPTED,
                Some(filename) => filename,
            };
            if !silent {
                print!("Read finite field from {}\n", secret_key);
            }
            let f2m = &SecretKey::read_finite_field(secret_key)?;
            if !silent {
                print!("Read secret key from {}\n", secret_key);
            }
            let sk = SecretKey::read_secret_key(secret_key, f2, f2m)?;
            if !silent {
                print!("Read ciphertext from {}\n", ciphertext);
            }
            let c = RowVec::read_vector(ciphertext, f2)?;
            if !silent {
                print!("Decrypt ciphertext...");
            }
            let m = sk.decrypt(&c);
            if !silent {
                print!("ok\nWrite decrypted plaintext to {}\n", decrypted);
            }
            m.write(decrypted)?;
            Ok(())
        }
        "plaintext" => {
            let plaintext = match args.get(2) {
                None => PLAINTEXT,
                Some(filename) => filename,
            };
            let mut rng = rand::thread_rng();
            if !silent {
                print!("Generate plaintext...");
            }
            let p = RowVec::random(&mut rng, f2, GOPPA_K as usize);
            if !silent {
                print!("ok\nWrite plaintext to {}\n", plaintext);
            }
            p.write(plaintext)?;
            Ok(())
        }
        _ => panic!(
            "Unexpected command. Valid commands are 'keygen', 'encrypt', 'decrypt' \
             and 'plaintext'. Try --help for more information."
        ),
    }
}
