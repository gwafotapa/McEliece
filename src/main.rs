// use log::info;
use getopts::{Matches, Options};
use main_error::MainError;
use std::{
    env,
    // error::Error,
    // fs::File,
    io::{self, Write},
    // result,
};

use crypto::{PublicKey, SecretKey};
use finite_field::{F2m, F2};
// use goppa::Goppa;
use matrix::RowVec;
// use polynomial::Poly;

pub mod crypto;
mod finite_field;
pub mod goppa; // TODO: Warnings about dead code without the pub. Why ?
mod matrix;
mod polynomial;

// const GOPPA_M: u32 = 10; // F2^GOPPA_M is the underlying field of the code
// const GOPPA_N: u32 = 1024; // Code length
// const GOPPA_K: u32 = 524; // Code dimension
// const GOPPA_T: u32 = 50; // Code correction capacity

// const GOPPA_M: u32 = 11; // F2^GOPPA_M is the underlying field of the code
const GOPPA_N: u32 = 2048; // Code length
const GOPPA_N_MIN: u32 = 4;
const GOPPA_N_MAX: u32 = 2048;
// const GOPPA_K: u32 = 1278; // Code dimension
const GOPPA_T: u32 = 70; // Code correction capacity

// const GOPPA_M: u32 = 12; // F2^GOPPA_M is the underlying field of the code
// const GOPPA_N: u32 = 4096; // Code length
// const GOPPA_K: u32 = 2056; // Code dimension
// const GOPPA_T: u32 = 170; // Code correction capacity

const PLAINTEXT: &str = "plaintext.mce";
const CIPHERTEXT: &str = "ciphertext.mce";
const DECRYPTED: &str = "decrypted.mce";
const PUBLIC_KEY: &str = "public_key.mce";
const SECRET_KEY: &str = "secret_key.mce";

// type Result<T> = result::Result<T, Box<dyn Error>>;

fn get_program(path: &str) -> &str {
    let i = match path.rfind('/') {
        None => 0,
        Some(index) => index + 1,
    };
    &path[i..]
}

fn print_help(program: &str, opts: Options) {
    let brief = format!(
        "Usage: {0} keygen [PK] [SK]\n\
         {0} encrypt [PK] [PLAINTEXT] [CIPHERTEXT]\n\
         {0} decrypt [SK] [CIPHERTEXT] [DECRYPTED]\n\
         {0} plaintext [PK] [PLAINTEXT]\n\
         Encrypts information using the McEliece cryptosystem.",
        program
    );
    print!("{}", opts.usage(&brief));
}

fn print_if(verbose: bool, message: &str) {
    if verbose {
        print!("{}", message);
        io::stdout().flush().unwrap();
    }
}

// TODO: move functions keygen, encrypt, decrypt and plaintext to crypto.rs ?
// Rename them ?
fn keygen(
    pk_file: &str,
    sk_file: &str,
    m: u32,
    n: u32,
    t: u32,
    verbose: bool,
) -> Result<(), MainError> {
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    print_if(verbose, "Generating keys (pk, sk).....");
    // if verbose {
    //     print!("Generating keys (pk, sk).....");
    //     io::stdout().flush().unwrap();
    // }
    let (pk, sk) = crypto::keygen(f2, f2m, n, t);
    print_if(
        verbose,
        format!("ok\nWriting pk to {}.....", pk_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nWriting pk to {}.....", pk_file);
    //     io::stdout().flush().unwrap();
    // }
    pk.write(pk_file)?;
    print_if(
        verbose,
        format!("ok\nWriting sk to {}.....", sk_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nWriting sk to {}.....", sk_file);
    //     io::stdout().flush().unwrap();
    // }
    sk.write(sk_file)?;
    print_if(verbose, "ok\n");
    // if verbose {
    //     print!("ok\n");
    // }
    Ok(())
}

fn encrypt(
    pk_file: &str,
    ptxt_file: &str,
    ctxt_file: &str,
    verbose: bool,
) -> Result<(), MainError> {
    let f2 = &F2 {};
    print_if(
        verbose,
        format!("Reading public key from {}.....", pk_file).as_str(),
    );
    // if verbose {
    //     print!("Reading public key from {}.....", public_key);
    //     io::stdout().flush().unwrap();
    // }
    let pk = PublicKey::read_public_key(pk_file, f2)?;
    print_if(
        verbose,
        format!("ok\nReading plaintext from {}.....", ptxt_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nReading plaintext from {}.....", plaintext);
    //     io::stdout().flush().unwrap();
    // }
    let m = RowVec::read_vector(ptxt_file, f2)?;
    if pk.sgp().rows() != m.cols() {
        return Err("Plaintext length does not match code dimension from public key".into());
    }
    print_if(verbose, "ok\nEncrypting plaintext.....");
    // if verbose {
    //     print!("ok\nEncrypting plaintext.....");
    //     io::stdout().flush().unwrap();
    // }
    let c = pk.encrypt(&m);
    print_if(
        verbose,
        format!("ok\nWriting ciphertext to {}.....", ctxt_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nWriting ciphertext to {}", ciphertext);
    //     io::stdout().flush().unwrap();
    // }
    c.write(ctxt_file)?;
    print_if(verbose, "ok\n");
    // if verbose {
    //     print!("ok\n");
    // }
    Ok(())
}

fn decrypt(sk_file: &str, ctxt_file: &str, dec_file: &str, verbose: bool) -> Result<(), MainError> {
    let f2 = &F2 {};
    print_if(
        verbose,
        format!("Reading finite field from {}.....", sk_file).as_str(),
    );
    // if verbose {
    //     print!("Reading finite field from {}.....", sk_file);
    //     io::stdout().flush().unwrap();
    // }
    let f2m = &SecretKey::read_finite_field(sk_file)?;
    print_if(
        verbose,
        format!("ok\nReading secret key from {}.....", sk_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nReading secret key from {}.....", sk_file);
    //     io::stdout().flush().unwrap();
    // }
    let sk = SecretKey::read_secret_key(sk_file, f2, f2m)?;
    print_if(
        verbose,
        format!("ok\nReading ctxt_file from {}.....", ctxt_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nReading ciphertext from {}.....", ctxt_file);
    //     io::stdout().flush().unwrap();
    // }
    let c = RowVec::read_vector(ctxt_file, f2)?;
    if sk.p().len() != c.cols() {
        return Err("Ciphertext length does not match code length from secret key".into());
    }
    print_if(verbose, "ok\nDecrypting ciphertext.....");
    // if verbose {
    //     print!("ok\nDecrypting ciphertext.....");
    //     io::stdout().flush().unwrap();
    // }
    let m = sk.decrypt(&c);
    print_if(
        verbose,
        format!("ok\nWriting decrypted text to {}.....", dec_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nWriting decrypted text to {}.....", dec_file);
    //     io::stdout().flush().unwrap();
    // }
    m.write(dec_file)?;
    print_if(verbose, "ok\n");
    // if verbose {
    //     print!("ok\n");
    // }
    Ok(())
}

fn plaintext(pk_file: &str, ptxt_file: &str, verbose: bool) -> Result<(), MainError> {
    print_if(
        verbose,
        format!("Reading code dimension k from {}.....", pk_file).as_str(),
    );
    // if verbose {
    //     print!("Reading code dimension k from {}.....", pk_file);
    //     io::stdout().flush().unwrap();
    // }
    let k = PublicKey::read_code_dimension(pk_file)?;
    print_if(verbose, "ok\nGenerating plaintext.....");
    // if verbose {
    //     print!("ok\nGenerating plaintext.....");
    //     io::stdout().flush().unwrap();
    // }
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let p = RowVec::random(&mut rng, f2, k as usize);
    print_if(
        verbose,
        format!("ok\nWriting plaintext to {}.....", ptxt_file).as_str(),
    );
    // if verbose {
    //     print!("ok\nWriting plaintext to {}.....", ptxt_file);
    //     io::stdout().flush().unwrap();
    // }
    p.write(ptxt_file)?;
    print_if(verbose, "ok\n");
    // if verbose {
    //     print!("ok\n");
    // }
    Ok(())
}

fn get_params(matches: &Matches) -> Result<(u32, u32, u32), MainError> {
    let n = match matches.opt_str("n") {
        None => GOPPA_N,
        Some(length) => u32::from_str_radix(&length, 10)?,
    };
    if n < GOPPA_N_MIN || GOPPA_N_MAX < n {
        return Err(format!(
            "Code length n must be greater than {} and less than {}",
            GOPPA_N_MIN, GOPPA_N_MAX,
        )
        .into());
    }
    let q = n.next_power_of_two();
    let m = q.trailing_zeros();
    let t = match matches.opt_str("t") {
        None => GOPPA_T,
        Some(correction) => u32::from_str_radix(&correction, 10)?,
    };
    if t == 1 && n == q {
        return Err("t cannot be 1 if n is a power of 2 \
                    because L cannot contain a root of the goppa polynomial"
            .into());
    } else if n <= m * t {
        return Err("Code length n must be greater than m * t".into());
    }
    Ok((m, n, t))
}

fn main() -> Result<(), MainError> {
    env_logger::init();
    let args: Vec<String> = env::args().collect();
    let program = get_program(&args[0]);
    let mut opts = Options::new();
    opts.optflag("h", "help", "Display this help");
    opts.optopt("n", "length", "Set Goppa code length", "N");
    opts.optopt("t", "correction", "Set Goppa code correction capacity", "T");
    opts.optflag("v", "verbose", "Explain what is being done");
    let matches = opts.parse(&args[1..]).map_err(|e| e.to_string())?;
    if matches.opt_present("h") {
        print_help(program, opts);
        return Ok(());
    }
    let verbose = matches.opt_present("v");
    let (m, n, t) = get_params(&matches)?;

    // println!("Code length: {}\nCode correction capacity: {}\n", n, t);
    // let f2 = &F2 {};

    let command = match matches.free.get(0) {
        Some(cmd) => cmd.as_str(),
        None => {
            return Err(format!(
                "Command expected\n\
             Try '{} --help' for more information.",
                program
            )
            .into())
        }
    };
    let files: Vec<&str> = matches.free.iter().skip(1).map(|s| s.as_str()).collect();

    match command {
        "keygen" => {
            let pk_file = files.get(0).unwrap_or(&PUBLIC_KEY);
            let sk_file = files.get(1).unwrap_or(&SECRET_KEY);
            keygen(pk_file, sk_file, m, n, t, verbose)
        }
        "encrypt" => {
            let pk_file = files.get(0).unwrap_or(&PUBLIC_KEY);
            let ptxt_file = files.get(1).unwrap_or(&PLAINTEXT);
            let ctxt_file = files.get(2).unwrap_or(&CIPHERTEXT);
            encrypt(pk_file, ptxt_file, ctxt_file, verbose)
        }
        "decrypt" => {
            let sk_file = files.get(0).unwrap_or(&SECRET_KEY);
            let ctxt_file = files.get(1).unwrap_or(&CIPHERTEXT);
            let dec_file = files.get(2).unwrap_or(&DECRYPTED);
            decrypt(sk_file, ctxt_file, dec_file, verbose)
        }
        "plaintext" => {
            let pk_file = files.get(0).unwrap_or(&PUBLIC_KEY);
            let ptxt_file = files.get(1).unwrap_or(&PLAINTEXT);
            plaintext(pk_file, ptxt_file, verbose)
        }
        _ => Err(format!(
            "Unexpected command\n\
             Try '{} --help' for more information.",
            program
        )
        .into()),
    }
}
