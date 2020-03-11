use rand::{rngs::ThreadRng, Rng};
use std::fmt::{Formatter, Result};

use super::{Field, FiniteField};

const CARD: u32 = 7;

const EXP: [F7Elt; CARD as usize] = [1, 3, 2, 6, 4, 5, 1];

const LOG: [u32; CARD as usize] = [CARD as u32, 0, 2, 1, 4, 5, 3];

type F7Elt = u32;

#[derive(Eq)]
pub struct F7 {
    // exp: Vec<F7Elt>,
    // log: Vec<u32>,
}

impl PartialEq for F7 {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl Field for F7 {
    type FElt = F7Elt;

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        CARD
    }

    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        (a + b) % CARD
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        (CARD + a - b) % CARD
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        let modulo = |x| {
            if x >= CARD {
                x - (CARD - 1)
            } else {
                x
            }
        };

        if a == 0 || b == 0 {
            0
        } else {
            // self.exp[modulo(self.log[a as usize] + self.log[b as usize]) as usize]
            EXP[modulo(LOG[a as usize] + LOG[b as usize]) as usize]
        }
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        (CARD - a) % CARD
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            // Some(self.exp[(CARD - 1 - self.log[a as usize]) as usize])
            Some(EXP[(CARD - 1 - LOG[a as usize]) as usize])
        }
    }

    fn random(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, CARD)
    }

    fn to_string_debug(&self, a: Self::FElt) -> String {
        a.to_string()
    }
    
    fn to_string_display(&self, a: Self::FElt) -> String {
        match a {
            0 => "0".to_owned(),
            1 => "1".to_owned(),
            3 => "a".to_owned(),
            _ => format!("a^{}", LOG[a as usize]),
        }
    }
}

impl FiniteField for F7 {
    fn characteristic_exponent(&self) -> u32 {
        1
    }
    // fn order(&self) -> u32 {
    //     self.characteristic().pow(self.characteristic_exponent())
    // }

    fn exp(&self, n: u32) -> Self::FElt {
        // self.exp[n as usize]
        EXP[n as usize]
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            // Some(self.log[a as usize])
            Some(LOG[a as usize])
        }
    }
}







// #[derive(Clone, Copy, Eq, PartialEq)]
// pub struct F7(pub u32);

// const CARD: u32 = 7;

// const EXP: [F7; CARD as usize] = [F7(1), F7(3), F7(2), F7(6), F7(4), F7(5), F7(1)];

// const LOG: [u32; CARD as usize] = [CARD as u32, 0, 2, 1, 4, 5, 3];

// impl FiniteFieldElement for F7 {
//     fn characteristic_exponent() -> u32 {
//         1
//     }

//     fn exp(i: u32) -> Self {
//         EXP[i as usize]
//     }

//     fn log(self) -> Option<u32> {
//         if self.0 == 0 {
//             None
//         } else {
//             Some(LOG[self.0 as usize])
//         }
//     }

//     fn to_canonical_basis(self) -> u32 {
//         panic!("Method canonical_basis has not been implemented for type F7");
//     }
// }

// impl Field for F7 {
//     fn zero() -> Self {
//         Self(0)
//     }

//     fn one() -> Self {
//         Self(1)
//     }

//     fn characteristic() -> u32 {
//         7
//     }
// }

// impl Add for F7 {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         Self((self.0 + other.0) % CARD)
//     }
// }

// impl AddAssign for F7 {
//     fn add_assign(&mut self, other: Self) {
//         *self = *self + other;
//     }
// }

// impl Sub for F7 {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         Self((CARD + self.0 - other.0) % CARD)
//     }
// }

// impl SubAssign for F7 {
//     fn sub_assign(&mut self, other: Self) {
//         *self = *self - other;
//     }
// }

// impl Mul for F7 {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         let modulo = |a| if a >= CARD { a - (CARD - 1) } else { a };

//         if self == Self(0) || other == Self(0) {
//             Self(0)
//         } else {
//             // Self::exp(modulo(Self::log(self).unwrap() + Self::log(other).unwrap()))
//             EXP[modulo(LOG[self.0 as usize] + LOG[other.0 as usize]) as usize]
//         }
//     }
// }

// impl MulAssign for F7 {
//     fn mul_assign(&mut self, other: Self) {
//         *self = *self * other;
//         // if *self == Self(0) || other == Self(0) {
//         //     *self = Self(0);
//         // } else {
//         //     *self = Self::exp(modulo(
//         //         Self::log(*self).unwrap() + Self::log(other).unwrap(),
//         //     ))
//         // }
//     }
// }

// impl Neg for F7 {
//     type Output = Self;

//     fn neg(self) -> Self {
//         Self(0) - self
//         // match self {
//         //     Self(0) => Self(0),
//         //     Self(1) => Self(6),
//         //     Self(2) => Self(5),
//         //     Self(3) => Self(4),
//         //     Self(4) => Self(3),
//         //     Self(5) => Self(2),
//         //     Self(6) => Self(1),
//         //     _ => Self(0),
//         // }
//     }
// }

// impl Inv for F7 {
//     type Output = Self;

//     fn inv(self) -> Option<Self::Output> {
//         match self {
//             Self(0) => None,
//             _ => Some(EXP[(CARD - 1 - LOG[self.0 as usize]) as usize]),
//         }
//     }
// }

// impl Distribution<F7> for Standard {
//     fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F7 {
//         F7(rng.gen_range(0, CARD) as u32)
//     }
// }

// impl Debug for F7 {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         write!(f, "{:?}", self.0,)
//     }
// }

// impl Display for F7 {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         // write!(f, "{}", self.0,)
//         match self {
//             Self(0) => write!(f, "0"),
//             Self(1) => write!(f, "1"),
//             Self(2) => write!(f, "a"),
//             _ => write!(f, "a^{}", self.log().unwrap()),
//         }
//     }
// }

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn f7_add() {
        let f = &F7 {};
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
    fn f7_characteristic() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let z = f.zero();
        let mut s = f.zero();
        for _i in 0..7 {
            s = f.add(s, a);
        }
        
        assert_eq!(s, z);
    }
    
    #[test]
    fn f7_sub() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let z = f.zero();
        
        assert_eq!(f.sub(a, z), a);
        assert_eq!(f.sub(a, a), z);
    }

    #[test]
    fn f7_mul() {
        let f = &F7 {};
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
    fn f7_neg() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);
        let z = f.zero();

        assert_eq!(f.neg(z), z);
        assert_eq!(f.neg(f.neg(a)), a);
        assert_eq!(f.add(a, f.neg(b)), f.sub(a, b));
    }

    #[test]
    fn f7_inv() {
        let f = &F7 {};
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
    
