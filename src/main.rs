use getopts::{Matches, Options};
use main_error::MainError;
use std::env;

use mceliece::{
    crypto::{self, PublicKey, SecretKey},
    finite_field::Field,
    matrix::RowVec,
};

// TODO: remove commented code
// TODO: Rework Add trait and others using self consuming
// TODO: useless line in test report
// TODO: FieldTrait called Field and Field enum called ?

// const GOPPA_N_MIN: usize = 3;
// const GOPPA_N_MAX: usize = 1024;
// const GOPPA_N: usize = 1024; // Code length
// const GOPPA_T: usize = 50; // Code correction capacity
//                            // const GOPPA_K: usize = 524; // Code dimension

// const GOPPA_N_MIN: usize = 3;
// const GOPPA_N_MAX: usize = 2048;
// const GOPPA_N: usize = 2048; // Code length
// const GOPPA_T: usize = 70; // Code correction capacity
//                            // const GOPPA_K: usize = 1278; // Code dimension

const GOPPA_N_MIN: usize = 3;
const GOPPA_N_MAX: usize = 4096;
const GOPPA_N: usize = 4096; // Code length
const GOPPA_T: usize = 170; // Code correction capacity
                            // const GOPPA_K: usize = 2056; // Code dimension

const PLAINTEXT: &str = "plaintext.mce";
const CIPHERTEXT: &str = "ciphertext.mce";
const DECRYPTED: &str = "decrypted.mce";
const PUBLIC_KEY: &str = "public_key.mce";
const SECRET_KEY: &str = "secret_key.mce";

fn get_program(path: &str) -> &str {
    let i = match path.rfind('/') {
        None => 0,
        Some(index) => index + 1,
    };
    &path[i..]
}

fn get_code_params(matches: &Matches) -> Result<(usize, usize), MainError> {
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
        return Err("The ratio n/t is not large enough. Pick a larger n or a smaller t.".into());
    }
    Ok((n, t))
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
    opts.optflag("v", "verbose", "Detail created files");
    let matches = opts.parse(&args[1..]).map_err(|e| e.to_string())?;
    if matches.opt_present("h") {
        print_help(program, opts);
        return Ok(());
    }
    let verbose = matches.opt_present("v");
    let (n, t) = get_code_params(&matches)?;
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
            let (pk, sk) = crypto::keygen(n, t);
            pk.write(pk_file)?;
            sk.write(sk_file)?;
            if verbose {
                println!(
                    "Wrote public key to '{}'.\nWrote secret key to '{}'.",
                    pk_file, sk_file
                );
            }
            Ok(())
        }
        "encrypt" => {
            let pk_file = files.get(0).unwrap_or(&PUBLIC_KEY);
            let ptxt_file = files.get(1).unwrap_or(&PLAINTEXT);
            let ctxt_file = files.get(2).unwrap_or(&CIPHERTEXT);
            let pk = PublicKey::read_public_key(pk_file)?;
            let m = RowVec::read_vector(ptxt_file)?;
            if pk.sgp().rows() != m.cols() {
                return Err(
                    "Plaintext length does not match code dimension from public key".into(),
                );
            }
            let c = pk.encrypt(&m);
            c.write(ctxt_file)?;
            if verbose {
                println!("Wrote ciphertext to '{}'.", ctxt_file);
            }
            Ok(())
        }
        "decrypt" => {
            let sk_file = files.get(0).unwrap_or(&SECRET_KEY);
            let ctxt_file = files.get(1).unwrap_or(&CIPHERTEXT);
            let dec_file = files.get(2).unwrap_or(&DECRYPTED);
            let sk = SecretKey::read_secret_key(sk_file)?;
            let c = RowVec::read_vector(ctxt_file)?;
            if sk.p().len() != c.cols() {
                return Err("Ciphertext length does not match code length from secret key".into());
            }
            let m = sk.decrypt(&c);
            m.write(dec_file)?;
            if verbose {
                println!("Wrote decrypted text to '{}'.", dec_file);
            }
            Ok(())
        }
        "plaintext" => {
            let pk_file = files.get(0).unwrap_or(&PUBLIC_KEY);
            let ptxt_file = files.get(1).unwrap_or(&PLAINTEXT);
            let k = PublicKey::read_code_dimension(pk_file)?;
            let m = RowVec::random(Field::Parameters(()), k);
            m.write(ptxt_file)?;
            if verbose {
                println!("Wrote plaintext to '{}'.", ptxt_file);
            }
            Ok(())
        }
        _ => Err(format!(
            "Unexpected command\n\
             Try '{} --help' for more information.",
            program
        )
        .into()),
    }
}
