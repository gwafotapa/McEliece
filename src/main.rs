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

// type PublicKey<'a> = (Mat<'a, F2>, u32);
// type SecretKey<'a, 'b> = (Mat<'a, F2>, Goppa<'b, F2m>, Mat<'a, F2>);

fn save_public_key(pk: &PublicKey, file_name: &str) {
    println!("{}\n", pk.sgp);
    let mut f = File::create(file_name).expect("Unable to create file");
    // f.write_all(format!("{:02x}", pk.sgp.rows()).as_bytes());
    // f.write_all(format!("{:02x}", pk.sgp.cols()).as_bytes());
    // for i in 0..pk.sgp.rows() {
    //     for j in (0..pk.sgp.cols()).step_by(8) {
    //         let mut tmp = 0;
    //         for k in 0..8 {
    //             tmp |= pk.sgp[(i, j + k)] << 7 - k;
    //         }
    //         f.write_all(format!("{:02x}", tmp).as_bytes());
    //     }
    // }
    f.write_all(pk.sgp.to_hex_string().as_bytes());
    f.write_all(format!("\n{:02x}\n", pk.t).as_bytes());
}

fn save_secret_key(sk: &SecretKey, file_name: &str) {
    let mut f = File::create(file_name).expect("Unable to create file");
    f.write_all((sk.s.to_hex_string() + "\n").as_bytes());
    f.write_all((sk.goppa.to_hex_string() + "\n").as_bytes());
    f.write_all((sk.p.to_hex_string() + "\n").as_bytes());
}

fn load_public_key<'a>(file_name: &str, f2: &'a F2) -> PublicKey<'a> {
    let mut f = File::open(file_name).expect("Unable to open file");
    let mut f = BufReader::new(f);
    // let mut buf = String::new();
    // let mut lines = f.lines();
    // f.read_to_string(&mut data).expect("Unable to read data");
    // buf.pop(); // Remove terminating newline character
    let mut line = String::new();
    f.read_line(&mut line).expect("Unable to read the first line");
    line.pop(); // Remove terminating newline
    let sgp = Mat::from_hex_string(&line, f2);
                                                                         // let mut data = hex::decode(data).expect("Hex decoding failed");
                                                                         // let t = data.pop().unwrap() as u32;
                                                                         // let rows = data[0];
                                                                         // let cols = data[1];
                                                                         // println!("rows: {}\ncols: {}\n", rows, cols);
                                                                         // let mut sgp = Mat::zero(f2, rows.into(), cols as usize); // TODO: into() or as ?
                                                                         // for i in 2..data.len() {
                                                                         //     for j in 0..8 {
                                                                         //         let bit_index = 8 * (i - 2) + j;
                                                                         //         let row = bit_index / cols as usize;
                                                                         //         let col = bit_index % cols as usize;
                                                                         //         sgp[(row, col)] = ((data[i] >> 7 - j) & 1).into();
                                                                         //     }
                                                                         // }
    println!("{}\n", sgp);
    // f.readline(&mut buf).expect("Cannot read second line");
    // buf.pop();
    line.clear();
    f.read_line(&mut line).expect("Unable to read the second line");
    line.pop(); // Remove terminating newline
    let t = u32::from_str_radix(&line, 16).unwrap();
    PublicKey { sgp, t }
}

fn load_finite_field(file_name: &str) -> F2m {
    let mut f = File::open(file_name).expect("Unable to open file"); // TODO: return result and use ?
    let mut f = BufReader::new(f);
    let mut line = String::new();
    f.read_line(&mut line).expect("Unable to read the first line");
    line.clear();
    f.read_line(&mut line).expect("Unable to read the second line");
    let order = u32::from_str_radix(line.split('#').next().unwrap(), 16).unwrap();
    let f = F2m::generate(order);
    f
}

fn load_secret_key<'a, 'b>(file_name: &str, f2: &'a F2, f2m: &'b F2m) -> SecretKey<'a, 'b> {
    let mut f = File::open(file_name).expect("Unable to open file"); // TODO: return result and use ?
    let f = BufReader::new(f);
    let mut lines = f.lines();
    let s = Mat::from_hex_string(&lines.next().unwrap().unwrap(), f2); // TODO: double unwrap ??
    let goppa = Goppa::from_hex_string(&lines.next().unwrap().unwrap(), f2m); // TODO: double unwrap ??
    let p = Mat::from_hex_string(&lines.next().unwrap().unwrap(), f2); // TODO: double unwrap ??
    SecretKey { s, goppa, p }
}

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
    let mut f = File::open(file_name).expect("Unable to open file"); // TODO: return result and use ?
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
    // TODO: deal with all the unwrap()s
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

    // // TODO: Check the length of data to avoid mul overflow ?
    // let data = fs::read(&args[1]).expect("Unable to read file");
    // println!("{:?}", data);
    // let k = 8 * data.len();
    // if k != 72 {
    //     panic!("Plaintext is {} bits and should be 72", k);
    // }
    // let f2 = &F2 {};
    // let mut msg = RowVec::zero(f2, k);
    // for i in 0..(data.len() - 1) {
    //     let nbr = char::from(data[i]).to_digit(16).unwrap();
    //     println!("{} ", nbr);
    //     for j in 0..8 {
    //         msg[8 * i + j] = ((nbr >> (7 - j)) & 1).into();
    //     }
    // }
    // println!("{:?}", msg);

    // let file_name = args[1];

    // let mut file = match OpenOptions::new()
    //     .write(true)
    //     .create_new(true)
    //     .open(file_name)
    // {
    //     Ok(f) => f,
    //     Err(e) => panic!("Cannot create file: {}", e),
    // };
}
