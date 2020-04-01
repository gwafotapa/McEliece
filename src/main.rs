use getopts::{Matches, Options};
use main_error::MainError;
use std::{
    env,
    io::{self, Write},
};

use mceliece::{
    crypto::{self, PublicKey, SecretKey},
    finite_field::{F2m, F2},
    matrix::RowVec,
};

// TODO: remove commented code

// const GOPPA_N_MIN: u32 = 3;
// const GOPPA_N_MAX: u32 = 1024;
// const GOPPA_N: u32 = 1024; // Code length
// const GOPPA_T: u32 = 50; // Code correction capacity
// const GOPPA_K: u32 = 524; // Code dimension

// const GOPPA_N_MIN: u32 = 3;
// const GOPPA_N_MAX: u32 = 2048;
// const GOPPA_N: u32 = 2048; // Code length
// const GOPPA_T: u32 = 70; // Code correction capacity
// // const GOPPA_K: u32 = 1278; // Code dimension

const GOPPA_N_MIN: usize = 3;
const GOPPA_N_MAX: usize = 4096;
const GOPPA_N: usize = 4096; // Code length
const GOPPA_T: usize = 170; // Code correction capacity
                            // const GOPPA_K: u32 = 2056; // Code dimension

const PLAINTEXT: &str = "plaintext.mce";
const CIPHERTEXT: &str = "ciphertext.mce";
const DECRYPTED: &str = "decrypted.mce";
const PUBLIC_KEY: &str = "public_key.mce";
const SECRET_KEY: &str = "secret_key.mce";

fn keygen(
    pk_file: &str,
    sk_file: &str,
    m: u32,
    n: usize,
    t: usize,
    verbose: bool,
) -> Result<(), MainError> {
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);

    if verbose {
        print!("Generating keys (pk, sk).....");
        io::stdout().flush().unwrap();
    }
    let (pk, sk) = crypto::keygen(f2, f2m, n, t);

    if verbose {
        print!("ok\nWriting pk to {}.....", pk_file);
        io::stdout().flush().unwrap();
    }
    pk.write(pk_file)?;

    if verbose {
        print!("ok\nWriting sk to {}.....", sk_file);
        io::stdout().flush().unwrap();
    }
    sk.write(sk_file)?;

    if verbose {
        print!("ok\n");
    }
    Ok(())
}

fn encrypt(
    pk_file: &str,
    ptxt_file: &str,
    ctxt_file: &str,
    verbose: bool,
) -> Result<(), MainError> {
    let f2 = &F2 {};

    if verbose {
        print!("Reading public key from {}.....", pk_file);
        io::stdout().flush().unwrap();
    }
    let pk = PublicKey::read_public_key(pk_file, f2)?;

    if verbose {
        print!("ok\nReading plaintext from {}.....", ptxt_file);
        io::stdout().flush().unwrap();
    }
    let m = RowVec::read_vector(ptxt_file, f2)?;
    if pk.sgp().rows() != m.cols() {
        return Err("Plaintext length does not match code dimension from public key".into());
    }

    if verbose {
        print!("ok\nEncrypting plaintext.....");
        io::stdout().flush().unwrap();
    }
    let c = pk.encrypt(&m);

    if verbose {
        print!("ok\nWriting ciphertext to {}", ctxt_file);
        io::stdout().flush().unwrap();
    }
    c.write(ctxt_file)?;

    if verbose {
        print!("ok\n");
    }
    Ok(())
}

fn decrypt(sk_file: &str, ctxt_file: &str, dec_file: &str, verbose: bool) -> Result<(), MainError> {
    let f2 = &F2 {};

    if verbose {
        print!("Reading finite field from {}.....", sk_file);
        io::stdout().flush().unwrap();
    }
    let f2m = &SecretKey::read_finite_field(sk_file)?;

    if verbose {
        print!("ok\nReading secret key from {}.....", sk_file);
        io::stdout().flush().unwrap();
    }
    let sk = SecretKey::read_secret_key(sk_file, f2, f2m)?;

    if verbose {
        print!("ok\nReading ciphertext from {}.....", ctxt_file);
        io::stdout().flush().unwrap();
    }
    let c = RowVec::read_vector(ctxt_file, f2)?;
    if sk.p().len() != c.cols() {
        return Err("Ciphertext length does not match code length from secret key".into());
    }

    if verbose {
        print!("ok\nDecrypting ciphertext.....");
        io::stdout().flush().unwrap();
    }
    let m = sk.decrypt(&c);

    if verbose {
        print!("ok\nWriting decrypted text to {}.....", dec_file);
        io::stdout().flush().unwrap();
    }
    m.write(dec_file)?;

    if verbose {
        print!("ok\n");
    }
    Ok(())
}

fn plaintext(pk_file: &str, ptxt_file: &str, verbose: bool) -> Result<(), MainError> {
    if verbose {
        print!("Reading code dimension k from {}.....", pk_file);
        io::stdout().flush().unwrap();
    }
    let k = PublicKey::read_code_dimension(pk_file)?;

    if verbose {
        print!("ok\nGenerating plaintext.....");
        io::stdout().flush().unwrap();
    }
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let p = RowVec::random(&mut rng, f2, k);

    if verbose {
        print!("ok\nWriting plaintext to {}.....", ptxt_file);
        io::stdout().flush().unwrap();
    }
    p.write(ptxt_file)?;

    if verbose {
        print!("ok\n");
    }
    Ok(())
}

fn get_program(path: &str) -> &str {
    let i = match path.rfind('/') {
        None => 0,
        Some(index) => index + 1,
    };
    &path[i..]
}

fn get_code_params(matches: &Matches) -> Result<(u32, usize, usize), MainError> {
    let n = match matches.opt_str("n") {
        None => GOPPA_N,
        Some(length) => u32::from_str_radix(&length, 10)? as usize,
    };
    if n < GOPPA_N_MIN || GOPPA_N_MAX < n {
        return Err(format!(
            "Code length n must be greater than {} and less than {}",
            GOPPA_N_MIN, GOPPA_N_MAX,
        )
        .into());
    }
    let t = match matches.opt_str("t") {
        None => GOPPA_T,
        Some(correction) => u32::from_str_radix(&correction, 10)? as usize,
    };
    let q = if t == 1 && n.is_power_of_two() {
        2 * n
    } else {
        n.next_power_of_two()
    };
    let m = q.trailing_zeros();
    if n <= m as usize * t {
        return Err("Code length n must be greater than m * t".into());
    }
    Ok((m, n, t))
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
    let (m, n, t) = get_code_params(&matches)?;
    if verbose {
        print!("Code length n: {}\nCode correction capacity t: {}\n", n, t);
    }
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

// TODO: add a function from_two_bytes ?
