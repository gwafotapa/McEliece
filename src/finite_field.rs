// use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
// use std::fmt;

pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}

pub trait Inv {
    fn inv(self) -> Option<Self>
    where
        Self: Sized;
}

pub trait Exp {
    fn exp(n: usize) -> Self;
}

pub trait Log {
    fn log(self) -> Option<usize>;
}

pub trait AsU32 {
    fn as_u32(self) -> u32;
}

// // Finite field element
// struct FFElt(u32);

// // Finite field of order 2^m
// struct FF {
//     order: u32,
//     // p: usize, // characteristic
//     m: u32,
//     exp: Vec<FFElt>,
//     log: Vec<u32>,
// }

// impl Zero for FFElt {
//     fn zero() -> FFElt {
//         FFElt(0)
//     }
// }

// impl One for FFElt {
//     fn one() -> FFElt {
//         FFElt(1)
//     }
// }

// impl Add for FFElt {
//     type Output = FFElt;

//     fn add(self, other: FFElt) -> FFElt {
//         FFElt(self.0 ^ other.0)
//     }
// }

// impl Sub for FFElt {
//     type Output = FFElt;

//     fn sub(self, other: FFElt) -> FFElt {
//         // FFElt(self.0 ^ other.0)
//         self.add(other)
//     }
// }

// impl AddAssign for FFElt {
//     fn add_assign(&mut self, other: FFElt) {
//         // *self = FFElt(self.0 ^ other.0);
//         *self = self + other;
//     }
// }

// impl SubAssign for FFElt {
//     fn sub_assign(&mut self, other: FFElt) {
//         // *self += other;
//         *self = self - other;
//     }
// }

// // impl Mul for FFElt {
// //     type Output = FFElt;

// //     fn mul(self, other: FFElt) -> FFElt {
// //         if self == FFElt(0) || other == FFElt(0) {
// //             FFElt(0)
// //         } else {
// //             FFElt::exp(modulo(
// //                 FFElt::log(self).unwrap() + FFElt::log(other).unwrap(),
// //             ))
// //         }
// //     }
// // }

// // impl MulAssign for FFElt {
// //     fn mul_assign(&mut self, other: FFElt) {
// //         if *self == FFElt(0) || other == FFElt(0) {
// //             *self = FFElt(0);
// //         } else {
// //             *self = FFElt::exp(modulo(
// //                 FFElt::log(*self).unwrap() + FFElt::log(other).unwrap(),
// //             ));
// //         }
// //     }
// // }

// impl fmt::Display for FFElt {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         write!(f, "{:>4}", self.0)
//     }
// }

// // impl distributions::Distribution<FFElt> for distributions::Standard {
// //     fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FFElt {
// //         FFElt(rng.gen_range(0, CARD) as u16)
// //     }
// // }

// // impl finite_field::Exp for FFElt {
// //     fn exp(i: usize) -> FFElt {
// //         EXP[i]
// //     }
// // }

// // impl finite_field::Log for FFElt {
// //     fn log(self) -> Option<usize> {
// //         if self == FFElt(0) {
// //             None
// //         } else {
// //             Some(LOG[self.0 as usize])
// //         }
// //     }
// // }

// impl Neg for FFElt {
//     type Output = FFElt;

//     fn neg(self) -> FFElt {
//         self
//     }
// }

// // impl finite_field::Inv for FFElt {
// //     fn inv(self) -> Option<FFElt> {
// //         match self {
// //             FFElt(0) => None,
// //             _ => Some(FFElt::exp(CARD - 1 - FFElt::log(self).unwrap())),
// //         }
// //     }
// // }

// impl FF {
//     fn primitive_poly(m: u32) -> u32 {
//         match m {
//             3 => 0xB,
//             4 => 0x11,
//             5 => 0x25,
//             6 => 0x43,
//             7 => 0x83,
//             8 => 0x11D,
//             9 => 0x211,
//             10 => 0x409,
//             11 => 0x805,
//             12 => 0x1053,
//             13 => 0x201B,
//             14 => 0x4143,
//             15 => 0x8003,
//             _ => panic!("m is too big"),
//         }
//     }

//     fn generate(m: u32) -> Self {
//         let order = 2_u32.pow(m);
//         let mut ff = Self {
//             order: order,
//             m: m,
//             exp: Vec::with_capacity(order as usize),
//             log: Vec::with_capacity(order as usize),
//         };
//         ff.log.resize(order as usize, 0);

//         let poly = Self::primitive_poly(m);
//         let mut elt = 1;
//         let gen = 2; // primitive element

//         for i in 0..order - 1 {
//             ff.exp.push(FFElt(elt));
//             ff.log[elt as usize] = i;

//             elt *= gen;
//             if elt >= order {
//                 elt ^= poly;
//             }
//         }

//         ff.exp[order as usize - 1] = FFElt(1);
//         ff.log[0] = order; // by convention
//         ff
//     }
// }
