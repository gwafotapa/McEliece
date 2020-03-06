// Is there a way to maintain the integrity of a polynomial (i.e. to ensure that the leading coefficient is not zero) ? With a trait like for characteristic 2 ?
// Add functions for rows and columns operations on matrices
// Try rr

// mod matrix;
// mod goppa;
// mod finite_field;
// mod finite_field_2;
// mod finite_field_7;
// mod finite_field_1024;
// mod matrix;

#[macro_use]
extern crate log;

mod finite_field;
mod finite_field_2;
mod finite_field_1024;
mod polynomial;

use crate::finite_field_2::F2;
use crate::polynomial::Poly;

// use std::path::Path;
use std::char;
use std::env;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{Write, Seek, SeekFrom};

fn main() {
    let args: Vec<String> = env::args().collect();

    println!("{}", &args[0]);

    if args.len() != 2 {
        panic!("One argument expected (the field order)");
    }

    let order = match u32::from_str_radix(&args[1], 10) {
        Ok(o) => o,
        Err(e) => panic!("Cannot read argument: {}", e),
    };
    let mut prime_factors = polynomial::trial_division(order);
    let m = prime_factors.len();
    prime_factors.dedup();
    if prime_factors.len() != 1 {
        panic!("Order must be a power of a prime");
    }
    let q = prime_factors[0];

    println!("{}^{}", q, m);
    // TODO: check bounds on q and m

    let file_name = "src/finite_field_".to_owned() + &args[1] + ".rs";

    // if Path::new(&file_name).exists() {
    //     panic!("File {} already exists", file_name);
    // }

    let mut file = match OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(file_name)
    {
        Ok(f) => f,
        Err(e) => panic!("Cannot create file: {}", e),
    };

    let content = "
use crate::finite_field;
use finite_field::FiniteFieldElement;

use rand::distributions;
use rand::Rng;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

macro_rules! array_init {
    ( $( $x:expr ),+ ) => {
        [ $( F***($x) ),+ ]
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct F***(u32);

const CARD: u32 = ***;

const EXP: [F***; CARD as usize] = array_init![
    ";
    
    let content = content.replace("***", &order.to_string());
    file.write_all(content.as_bytes()).expect("Cannot write file");
    
    // https://www.partow.net/programming/polynomials/index.html
    let primitive_poly = match m {
        2 => 0x7,
        3 => 0xB,
        4 => 0x11,
        5 => 0x25,
        6 => 0x43,
        7 => 0x83,
        8 => 0x11D,
        9 => 0x211,
        10 => 0x409,
        11 => 0x805,
        12 => 0x1053,
        13 => 0x201B,
        14 => 0x4143,
        15 => 0x8003,
        16 => 0x110B,
        17 => 0x2009,
        _ => panic!("m must be at least 2 and at most 17"),
    };

    let mut log = vec![0_u32; order as usize];
    let mut elt = 1_u32;
    file.write_all(b"1, ").expect("Cannot write file");
    for i in 1..order {
        elt *= 2;
        if elt >= order {
            elt ^= primitive_poly;
        }
        file.write_all((elt.to_string() + ", ").as_bytes()).expect("Cannot write file");
            
        log[elt as usize] = i;
    }
    log[1] = 0;
    file.seek(SeekFrom::Current(-2)).expect("Cannot seek file"); // erase the last comma
    file.write_all(b"\n];\n\n").expect("Cannot write file");

    file.write_all(b"const LOG: [u32; CARD as usize] = [\n    CARD, ").expect("Cannot write file");
        
    for i in 1..order as usize {
        file.write_all((log[i].to_string() + ", ").as_bytes()).expect("Cannot write file");
            
    }
    file.seek(SeekFrom::Current(-2)).expect("Cannot seek file"); // erase the last comma
    file.write_all(b"\n];\n").expect("Cannot seek file");

    let content = "
fn modulo(a: u32) -> u32 {
    if a >= CARD {
        a - (CARD - 1)
    } else {
        a
    }
}

impl fmt::Debug for F*** {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, \"{}\", self.0)
    }
}

impl finite_field::CharacteristicTwo for F*** {}

impl finite_field::FiniteFieldElement for F*** {
    fn finite_field_q() -> u32 {
        2
    }

    fn finite_field_m() -> u32 {
        10
    }

    fn zero() -> Self {
        Self(0)
    }

    fn one() -> Self {
        Self(1)
    }

    fn exp(i: u32) -> Self {
        EXP[i as usize]
    }

    fn log(self) -> Option<u32> {
        if self == Self(0) {
            None
        } else {
            Some(LOG[self.0 as usize])
        }
    }

    fn to_u32(self) -> u32 {
        self.0
    }

    fn inv(self) -> Option<Self> {
        match self {
            Self(0) => None,
            _ => Some(Self::exp(CARD - 1 - Self::log(self).unwrap())),
        }
    }
}

impl Neg for F*** {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for F*** {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl Sub for F*** {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for F*** {
    fn add_assign(&mut self, other: Self) {
        *self = Self(self.0 ^ other.0);
    }
}

impl SubAssign for F*** {
    fn sub_assign(&mut self, other: Self) {
        *self += other;
    }
}

impl Mul for F*** {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            Self::exp(modulo(Self::log(self).unwrap() + Self::log(other).unwrap()))
        }
    }
}

impl MulAssign for F*** {
    fn mul_assign(&mut self, other: Self) {
        if *self == Self(0) || other == Self(0) {
            *self = Self(0);
        } else {
            *self = Self::exp(modulo(
                Self::log(*self).unwrap() + Self::log(other).unwrap(),
            ));
        }
    }
}

impl fmt::Display for F*** {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, \"{:>4}\", self.0,)
    }
}

impl distributions::Distribution<F***> for distributions::Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F*** {
        F***(rng.gen_range(0, CARD) as u32)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::finite_field::FiniteFieldElement;

    #[test]
    fn f***_add() {
        let mut rng = rand::thread_rng();
        let a: F*** = rng.gen();
        let b: F*** = rng.gen();
        let c: F*** = rng.gen();
        let z = F***::zero();
        assert_eq!(a + (b + c), (a + b) + c);
        assert_eq!(a + b, b + a);
        assert_eq!(a + z, a);
    }

    #[test]
    fn f***_sub() {
        let mut rng = rand::thread_rng();
        let a: F*** = rng.gen();
        let z = F***::zero();
        assert_eq!(a - z, a);
        assert_eq!(a - a, z);
    }

    #[test]
    fn f***_mul() {
        let mut rng = rand::thread_rng();
        let a: F*** = rng.gen();
        let b: F*** = rng.gen();
        let c: F*** = rng.gen();
        let i = F***::one();
        let z = F***::zero();
        assert_eq!(a * (b * c), (a * b) * c);
        assert_eq!(a * b, b * a);
        assert_eq!(a * i, a);
        assert_eq!(a * z, z);
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a * (b - c), a * b - a * c);
    }

    #[test]
    fn f***_inv() {
        let mut rng = rand::thread_rng();
        let a: F*** = rng.gen();
        let i = F***::one();
        let z = F***::zero();
        assert_eq!(z.inv(), None);
        assert_eq!(i.inv(), Some(i));
        if a != F***::zero() {
            assert_eq!(a.inv().unwrap().inv().unwrap(), a);
            assert_eq!(a * a.inv().unwrap(), i);
        }
    }

    #[test]
    fn f***_neg() {
        let mut rng = rand::thread_rng();
        let a: F*** = rng.gen();
        let b: F*** = rng.gen();
        let z = F***::zero();
        assert_eq!(-z, z);
        assert_eq!(--a, a);
        assert_eq!(a + -b, a - b);
    }
}";

    let content = content.replace("***", &order.to_string());
    file.write_all(content.as_bytes()).expect("Cannot write file");
}

// TODO: need to open lib.rs and add the line "mod finite_field_{q};"
// modify functions to return proper values (q instead of 2 for characteristic...)
// make changes for odd characteristic (neg, ...)
