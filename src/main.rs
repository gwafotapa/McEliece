// // Try rr

// use std::{
//     env,
//     fs::OpenOptions,
//     io::{Seek, SeekFrom, Write},
// };

// mod finite_field;
// // mod matrix;
// mod polynomial;
// // mod goppa;

// fn main() {
//     let args: Vec<String> = env::args().collect();

//     println!("{}", &args[0]);

//     if args.len() != 2 {
//         panic!("One argument expected (the field order)");
//     }

//     let order = match u32::from_str_radix(&args[1], 10) {
//         Ok(o) => o,
//         Err(e) => panic!("Cannot read argument: {}", e),
//     };
//     let mut prime_factors = polynomial::characteristic_two::trial_division(order);
//     let m = prime_factors.len();
//     prime_factors.dedup();
//     if prime_factors.len() != 1 {
//         panic!("Order must be a power of a prime");
//     }
//     let p = prime_factors[0];
//     if p != 2 {
//         panic!("Only characteristic 2 is supported");
//     }
//     if m < 2 || m > 17 {
//         panic!("Field order must be between 2^2 and 2^17");
//     }
//     println!("{}^{}", p, m);

//     let file_name = "src/finite_field_".to_owned() + &args[1] + ".rs";

//     let mut file = match OpenOptions::new()
//         .write(true)
//         .create_new(true)
//         .open(file_name)
//     {
//         Ok(f) => f,
//         Err(e) => panic!("Cannot create file: {}", e),
//     };

//     let content =
//         "use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, Inv};

// use rand::{
//     distributions::{Distribution, Standard},
//     Rng,
// };
// use std::{
//     fmt::{Debug, Display, Formatter, Result},
//     ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
// };

// macro_rules! array_init {
//     ( $( $x:expr ),+ ) => {
//         [ $( F{order}($x) ),+ ]
//     }
// }

// #[derive(Clone, Copy, Eq, PartialEq)]
// pub struct F{order}(pub u32);

// const CARD: u32 = {order};

// const EXP: [F{order}; CARD as usize] = array_init![
//     ";

//     let content = content.replace("{order}", &order.to_string());
//     file.write_all(content.as_bytes())
//         .expect("Cannot write file");

//     // https://www.partow.net/programming/polynomials/index.html
//     let primitive_poly = match m {
//         2 => 0x7,
//         3 => 0xB,
//         4 => 0x11,
//         5 => 0x25,
//         6 => 0x43,
//         7 => 0x83,
//         8 => 0x11D,
//         9 => 0x211,
//         10 => 0x409,
//         11 => 0x805,
//         12 => 0x1053,
//         13 => 0x201B,
//         14 => 0x4143,
//         15 => 0x8003,
//         16 => 0x110B,
//         17 => 0x2009,
//         _ => panic!("m must be at least 2 and at most 17"),
//     };

//     let mut log = vec![0_u32; order as usize];
//     let mut elt = 1_u32;
//     file.write_all(b"1, ").expect("Cannot write file");
//     for i in 1..order {
//         elt *= 2;
//         if elt >= order {
//             elt ^= primitive_poly;
//         }
//         file.write_all((elt.to_string() + ", ").as_bytes())
//             .expect("Cannot write file");

//         log[elt as usize] = i;
//     }
//     log[1] = 0;
//     file.seek(SeekFrom::Current(-2)).expect("Cannot seek file"); // erase the last comma
//     file.write_all(b"\n];\n\n").expect("Cannot write file");

//     file.write_all(b"const LOG: [u32; CARD as usize] = [\n    CARD, ")
//         .expect("Cannot write file");

//     for i in 1..order as usize {
//         file.write_all((log[i].to_string() + ", ").as_bytes())
//             .expect("Cannot write file");
//     }
//     file.seek(SeekFrom::Current(-2)).expect("Cannot seek file"); // erase the last comma
//     file.write_all(b"\n];\n").expect("Cannot seek file");

//     let content = "
// impl CharacteristicTwo for F{order} {}

// impl FiniteFieldElement for F{order} {
//     fn characteristic_exponent() -> u32 {
//         {characteristic_exponent}
//     }

//     fn exp(i: u32) -> Self {
//         EXP[i as usize]
//     }

//     fn log(self) -> Option<u32> {
//         if self == Self(0) {
//             None
//         } else {
//             Some(LOG[self.0 as usize])
//         }
//     }

//     fn to_canonical_basis(self) -> u32 {
//         self.0
//     }
// }

// impl FieldElement for F{order} {
//     fn zero() -> Self {
//         Self(0)
//     }

//     fn one() -> Self {
//         Self(1)
//     }

//     fn characteristic() -> u32 {
//         2
//     }
// }

// impl Add for F{order} {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         Self(self.0 ^ other.0)
//     }
// }

// impl AddAssign for F{order} {
//     fn add_assign(&mut self, other: Self) {
//         *self = *self + other;
//     }
// }

// impl Sub for F{order} {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         self + other
//     }
// }

// impl SubAssign for F{order} {
//     fn sub_assign(&mut self, other: Self) {
//         *self = *self - other;
//     }
// }

// impl Mul for F{order} {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         let modulo = |a| if a >= CARD { a - (CARD - 1) } else { a };

//         if self == Self(0) || other == Self(0) {
//             Self(0)
//         } else {
//             EXP[modulo(LOG[self.0 as usize] + LOG[other.0 as usize]) as usize]
//         }
//     }
// }

// impl MulAssign for F{order} {
//     fn mul_assign(&mut self, other: Self) {
//         *self = *self * other;
//     }
// }

// impl Neg for F{order} {
//     type Output = Self;

//     fn neg(self) -> Self {
//         self
//     }
// }

// impl Inv for F{order} {
//     type Output = Self;

//     fn inv(self) -> Option<Self::Output> {
//         match self {
//             Self(0) => None,
//             _ => Some(EXP[(CARD - 1 - LOG[self.0 as usize]) as usize]),
//         }
//     }
// }

// impl Distribution<F{order}> for Standard {
//     fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F{order} {
//         F{order}(rng.gen_range(0, CARD) as u32)
//     }
// }

// impl Debug for F{order} {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         write!(f, \"{:b}\", self.0)
//     }
// }

// impl Display for F{order} {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         match self {
//             Self(0) => write!(f, \"0\"),
//             Self(1) => write!(f, \"1\"),
//             Self(2) => write!(f, \"a\"),
//             _ => write!(f, \"a^{}\", self.log().unwrap()),
//         }
//     }
// }

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn f{order}_add() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let b: F{order} = rng.gen();
//         let c: F{order} = rng.gen();
//         let z = F{order}::zero();
//         assert_eq!(a + (b + c), (a + b) + c);
//         assert_eq!(a + b, b + a);
//         assert_eq!(a + z, a);
//     }

//     #[test]
//     fn f{order}_characteristic() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let z = F{order}::zero();
//         assert_eq!(a + a, z);
//     }

//     #[test]
//     fn f{order}_sub() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let b: F{order} = rng.gen();
//         assert_eq!(a + b, a - b);
//     }

//     #[test]
//     fn f{order}_mul() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let b: F{order} = rng.gen();
//         let c: F{order} = rng.gen();
//         let i = F{order}::one();
//         let z = F{order}::zero();
//         assert_eq!(a * (b * c), (a * b) * c);
//         assert_eq!(a * b, b * a);
//         assert_eq!(a * i, a);
//         assert_eq!(a * z, z);
//         assert_eq!(a * (b + c), a * b + a * c);
//         assert_eq!(a * (b - c), a * b - a * c);
//     }

//     #[test]
//     fn f{order}_neg() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let b: F{order} = rng.gen();
//         let z = F{order}::zero();
//         assert_eq!(-z, z);
//         assert_eq!(--a, a);
//         assert_eq!(a + -b, a - b);
//     }

//     #[test]
//     fn f{order}_inv() {
//         let mut rng = rand::thread_rng();
//         let a: F{order} = rng.gen();
//         let i = F{order}::one();
//         let z = F{order}::zero();
//         assert_eq!(z.inv(), None);
//         assert_eq!(i.inv(), Some(i));
//         if a != F{order}::zero() {
//             assert_eq!(a.inv().unwrap().inv().unwrap(), a);
//             assert_eq!(a * a.inv().unwrap(), i);
//         }
//     }
// }";

//     let content = content.replace("{order}", &order.to_string());
//     // let content = content.replace("{characteristic}", &p.to_string());
//     let content = content.replace("{characteristic_exponent}", &m.to_string());
//     file.write_all(content.as_bytes())
//         .expect("Cannot write file");
// }

// // TODO: need to open lib.rs and add the line "mod finite_field_{q};"

mod finite_field;

fn main() {}
