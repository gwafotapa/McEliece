use rand::{rngs::ThreadRng, Rng};
use std::fmt::{Formatter, Result};

use super::{Field, FiniteField};

type F2mElt = u32; // Finite Field Element

#[derive(Eq)]
pub struct F2m {
    order: u32,
    m: u32, // characteristic exponent
    exp: Vec<F2mElt>,
    log: Vec<u32>,
}

impl PartialEq for F2m {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}

impl Field for F2m {
    type FElt = F2mElt;

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        2
    }

    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a ^ b
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        self.add(a, b)
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        let modulo = |a| {
            if a >= self.order {
                a - (self.order - 1)
            } else {
                a
            }
        };

        if a == 0 || b == 0 {
            0
        } else {
            self.exp[modulo(self.log[a as usize] + self.log[b as usize]) as usize]
        }
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        a
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            Some(self.exp[(self.order - 1 - self.log[a as usize]) as usize])
        }
    }

    fn random(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, self.order)
    }

    fn to_string_debug(&self, a: Self::FElt) -> String {
        format!("{:b}", a)
    }
    
    fn to_string_display(&self, a: Self::FElt) -> String {
        match a {
            0 => "0".to_owned(),
            1 => "1".to_owned(),
            2 => "a".to_owned(),
            _ => format!("a^{}", self.log[a as usize])
        }
    }
}

impl FiniteField for F2m {
    fn characteristic_exponent(&self) -> u32 {
        self.m
    }
    // fn order(&self) -> u32 {
    //     self.characteristic().pow(self.characteristic_exponent())
    // }

    fn exp(&self, n: u32) -> Self::FElt {
        self.exp[n as usize]
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            Some(self.log[a as usize])
        }
    }
}

// impl FF {
//     pub fn order(&self) -> u32 {
//         self.order
//     }

//     pub fn characteristic(&self) -> u32 {
//         self.p
//     }

//     pub fn characteristic_exponent(&self) -> u32 {
//         self.m
//     }

//     pub fn exp(&self, i: u32) -> FFElt {
//         FFElt {
//             ff: self,
//             val: self.exp[(i % (self.order - 1)) as usize],
//         }
//     }

//     pub fn log(&self, elt: FFElt) -> Option<u32> {
//         match elt.val {
//             0 => None,
//             x => Some(self.log[x as usize]),
//         }
//     }

//     pub fn elt(&self, val: u32) -> FFElt {
//         if val > self.order {
//             panic!("Element does not live in finite field");
//         }
//         FFElt { ff: self, val }
//     }

impl F2m {
    fn primitive_poly(order: u32) -> u32 {
        match prime_power(order) {
            Ok((_, m)) => match m {
                // 1 => 0x2,
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
            },
            Err(s) => panic!(s),
        }
    }

    pub fn generate(order: u32) -> Self {
        let (_, m) = match prime_power(order) {
            Ok(r) => r,
            Err(s) => panic!(s),
        };
        let mut f = Self {
            order,
            m,
            exp: vec![0; order as usize],
            log: vec![0; order as usize],
        };
        f.exp[0] = 1;
        f.log[1] = 0;
        let mut elt = 1;
        for i in 1..order {
            elt *= 2;
            if elt >= order {
                elt ^= Self::primitive_poly(order);
            }
            f.exp[i as usize] = elt;
            f.log[elt as usize] = i;
        }
        f
    }

    // pub fn random(&self, rng: &mut ThreadRng) -> F2mElt {
    //     rng.gen_range(0, self.order)
    // }

    // pub fn fmt(&self, f: &mut Formatter<'_>) -> Result {

    // }
}

// impl<'a> CharacteristicTwo for FFElt<'a> {}

// impl<'a> FiniteFieldElement for FFElt<'a> {
//     fn characteristic_exponent() -> u32 {
//         3
//     }

//     fn exp(i: u32) -> Self {
//         unsafe { EXP[i as usize] }
//     }

//     fn log(self) -> Option<u32> {
//         if self == Self(0) {
//             None
//         } else {
//             unsafe { Some(LOG[self.0 as usize]) }
//         }
//     }

//     fn to_canonical_basis(self) -> u32 {
//         self.0
//     }
// }

// impl<'a> FieldElement for FFElt<'a> {
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

// impl<'a> Add for FFElt<'a> {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         Self {
//             ff: self.ff,
//             val: self.val ^ other.val,
//         }
//     }
// }

// impl<'a> AddAssign for FFElt<'a> {
//     fn add_assign(&mut self, other: Self) {
//         *self = *self + other;
//     }
// }

// impl<'a> Sub for FFElt<'a> {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         self + other
//     }
// }

// impl<'a> SubAssign for FFElt<'a> {
//     fn sub_assign(&mut self, other: Self) {
//         *self = *self - other;
//     }
// }

// impl<'a> Mul for FFElt<'a> {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         let ff = self.ff;
//         let modulo = |a| {
//             if a >= ff.order {
//                 a - (ff.order - 1)
//             } else {
//                 a
//             }
//         };

//         if self.val == 0 || other.val == 0 {
//             Self { ff: ff, val: 0 }
//         } else {
//             Self {
//                 ff: ff,
//                 val: ff.exp
//                     [modulo(ff.log[self.val as usize] + ff.log[other.val as usize]) as usize],
//             }
//         }
//     }
// }

// impl<'a> MulAssign for FFElt<'a> {
//     fn mul_assign(&mut self, other: Self) {
//         *self = *self * other;
//     }
// }

// impl<'a> Neg for FFElt<'a> {
//     type Output = Self;

//     fn neg(self) -> Self {
//         self
//     }
// }

// impl<'a> Inv for FFElt<'a> {
//     type Output = Self;

//     fn inv(self) -> Option<Self::Output> {
//         match self.val {
//             0 => None,
//             _ => Some(Self {
//                 ff: self.ff,
//                 val: self.ff.exp[(self.ff.order - 1 - self.ff.log[self.val as usize]) as usize],
//             }),
//         }
//     }
// }

// impl<'a> Distribution<FFElt<'a>> for Standard {
//     fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FFElt<'a> {
//         unsafe { FFElt<'a>(rng.gen_range(0, CARD) as u32) }
//     }
// }

// impl<'a> Debug for FFElt<'a> {
//     fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{:b}", self.val)
//     }
// }

// impl<'a> Display for FFElt<'a> {
//     fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
//         match self.val {
//             0 => write!(f, "0"),
//             1 => write!(f, "1"),
//             2 => write!(f, "a"),
//             _ => write!(f, "a^{}", self.ff.log(*self).unwrap()),
//         }
//     }
// }

pub fn prime_power(q: u32) -> std::result::Result<(u32, u32), &'static str> {
    let mut prime_factors = trial_division(q);
    let m = prime_factors.len() as u32;
    prime_factors.dedup();
    if prime_factors.len() != 1 {
        return Err("Number is not prime");
    }
    let p = prime_factors[0];
    Ok((p, m))
}

/// Computes the prime factors of (non zero) integer n by trial division
/// https://en.wikipedia.org/wiki/Trial_division
/// ```
/// # use mceliece::finite_field::ff2m::trial_division;
/// assert_eq!(trial_division(1), vec![]);
/// assert_eq!(trial_division(19), vec![19]);
/// assert_eq!(trial_division(77), vec![7, 11]);
/// assert_eq!(trial_division(12), vec![2, 2, 3]);
/// ```
pub fn trial_division(mut n: u32) -> Vec<u32> {
    if n == 0 {
        panic!("0 is an invalid input for trial division");
    }

    let mut prime_factors = Vec::new();
    while n % 2 == 0 {
        prime_factors.push(2);
        n /= 2;
    }
    let mut f = 3;
    while f * f <= n {
        if n % f == 0 {
            prime_factors.push(f);
            n /= f;
        } else {
            f += 2;
        }
    }
    if n != 1 {
        prime_factors.push(n);
    }
    prime_factors
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn f256_add() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);
        let c = f.random(&mut rng);
        let z = f.zero();

        assert_eq!(f.add(a, f.add(b, c)), f.add(f.add(a, b), c));
        assert_eq!(f.add(a, b), f.add(b, a));
        assert_eq!(f.add(a, z), a);
    }

    #[test]
    fn f256_characteristic() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let z = f.zero();

        assert_eq!(f.add(a, a), z);
    }

    #[test]
    fn f256_sub() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);

        assert_eq!(f.add(a, b), f.sub(a, b));
    }

    #[test]
    fn f256_mul() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);
        let c = f.random(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.mul(a, f.mul(b, c)), f.mul(f.mul(a, b), c));
        assert_eq!(f.mul(a, b), f.mul(b, a));
        assert_eq!(f.mul(a, i), a);
        assert_eq!(f.mul(a, z), z);
        assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
        assert_eq!(f.mul(a, f.sub(b, c)), f.sub(f.mul(a, b), f.mul(a, c)));
    }

    #[test]
    fn f256_neg() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);
        let z = f.zero();

        assert_eq!(f.neg(z), z);
        assert_eq!(f.neg(f.neg(a)), a);
        assert_eq!(f.add(a, f.neg(b)), f.sub(a, b));
    }

    #[test]
    fn f256_inv() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.inv(z), None);
        assert_eq!(f.inv(i), Some(i));
        if a != f.zero() {
            assert_eq!(f.inv(f.inv(a).unwrap()), Some(a));
            assert_eq!(f.mul(a, f.inv(a).unwrap()), i);
        }
    }
}
