// use rand::distributions::{Distribution, Standard};
use rand::Rng;
use std::fmt::{Formatter, Result};
// use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::rngs::ThreadRng;

// use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, Inv};
// use self::F2Elt::{One, Zero};
use super::{Field, FiniteField};

// #[derive(Clone, Copy, Eq, PartialEq)]
// pub enum F2Elt {
//     Zero,
//     One,
// }
pub type F2Elt = u32;

#[derive(Eq, PartialEq)]
pub struct F2;

impl Field for F2 {
    type FElt = F2Elt;

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        2
    }

    /// Adds two finite field elements
    /// ```
    /// # use mceliece::finite_field::Field;
    /// # use mceliece::finite_field::F2;
    ///
    /// let f2 = F2 {};
    /// assert_eq!(f2.add(0, 0), 0);
    /// assert_eq!(f2.add(0, 1), 1);
    /// assert_eq!(f2.add(1, 0), 1);
    /// assert_eq!(f2.add(1, 1), 0);
    /// ```    
    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a ^ b
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        self.add(a, b)
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a * b
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        a
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            Some(1)
        }
    }

    fn random(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, 2)
    }

    fn to_string_debug(&self, a: Self::FElt) -> String {
        a.to_string()
    }
    
    fn to_string_display(&self, a: Self::FElt) -> String {
        a.to_string()
    }
}

impl FiniteField for F2 {
    fn characteristic_exponent(&self) -> u32 {
        1
    }

    fn exp(&self, _n: u32) -> Self::FElt {
        1
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            Some(0)
        }
    }
}

// impl F2 {
//     pub fn random(&self, rng: &mut ThreadRng) -> F2Elt {
//         rng.gen_range(0, 2)
//     }
// }

// impl CharacteristicTwo for F2Elt {}

// impl FiniteFieldElement for F2Elt {
//     fn characteristic_exponent() -> u32 {
//         1
//     }

//     fn exp(_i: u32) -> Self {
//         Self::One
//     }

//     fn log(self) -> Option<u32> {
//         match self {
//             Self::Zero => None,
//             Self::One => Some(0),
//         }
//     }

//     fn to_canonical_basis(self) -> u32 {
//         match self {
//             Self::Zero => 0,
//             Self::One => 1,
//         }
//     }
// }

// impl Add for F2Elt {
//     type Output = Self;

//     /// Adds two finite field elements
//     /// ```
//     /// # use mceliece::finite_field::FieldElement;
//     /// # use mceliece::finite_field::F2Elt;
//     ///
//     /// assert_eq!(F2Elt::zero() + F2Elt::zero(), F2Elt::zero());
//     /// assert_eq!(F2Elt::zero() + F2Elt::one(), F2Elt::one());
//     /// assert_eq!(F2Elt::one() + F2Elt::zero(), F2Elt::one());
//     /// assert_eq!(F2Elt::one() + F2Elt::one(), F2Elt::zero());
//     /// ```

//     fn add(self, other: Self) -> Self {
//         match self {
//             Self::Zero => other,
//             Self::One => match other {
//                 Self::Zero => Self::One,
//                 Self::One => Self::Zero,
//             },
//         }
//     }
// }

// impl AddAssign for F2Elt {
//     fn add_assign(&mut self, other: Self) {
//         *self = *self + other;
//     }
// }

// impl Sub for F2Elt {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         self + other
//     }
// }

// impl SubAssign for F2Elt {
//     fn sub_assign(&mut self, other: Self) {
//         *self = *self - other;
//     }
// }

// impl Mul for F2Elt {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         match self {
//             Self::Zero => Self::Zero,
//             Self::One => other,
//         }
//     }
// }

// impl MulAssign for F2Elt {
//     fn mul_assign(&mut self, other: Self) {
//         *self = *self * other;
//     }
// }

// impl Neg for F2Elt {
//     type Output = Self;

//     fn neg(self) -> Self {
//         self
//     }
// }

// impl Inv for F2Elt {
//     type Output = Self;

//     fn inv(self) -> Option<Self::Output> {
//         match self {
//             Self::Zero => None,
//             Self::One => Some(Self::One),
//         }
//     }
// }

// impl Distribution<F2Elt> for Standard {
//     fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F2Elt {
//         match rng.gen_bool(0.5) {
//             false => F2Elt::Zero,
//             true => F2Elt::One,
//         }
//     }
// }

// impl Debug for F2Elt {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         write!(f, "{:?}", self)
//     }
// }

// impl Display for F2Elt {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         write!(f, "{}", self)
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f2_add() {
        let f2 = F2 {};
        assert_eq!(f2.add(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.add(f2.zero(), f2.one()), f2.one());
        assert_eq!(f2.add(f2.one(), f2.zero()), f2.one());
        assert_eq!(f2.add(f2.one(), f2.one()), f2.zero());
    }

    #[test]
    fn f2_sub() {
        let f2 = F2 {};
        assert_eq!(f2.sub(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.sub(f2.zero(), f2.one()), f2.one());
        assert_eq!(f2.sub(f2.one(), f2.zero()), f2.one());
        assert_eq!(f2.sub(f2.one(), f2.one()), f2.zero());
    }

    #[test]
    fn f2_mul() {
        let f2 = F2 {};
        assert_eq!(f2.mul(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.mul(f2.zero(), f2.one()), f2.zero());
        assert_eq!(f2.mul(f2.one(), f2.zero()), f2.zero());
        assert_eq!(f2.mul(f2.one(), f2.one()), f2.one());
    }

    #[test]
    fn f2_inv() {
        let f2 = F2 {};
        assert_eq!(f2.inv(f2.zero()), None);
        assert_eq!(f2.inv(f2.one()), Some(f2.one()));
    }

    #[test]
    fn f2_neg() {
        let f2 = F2 {};
        assert_eq!(f2.neg(f2.zero()), f2.zero());
        assert_eq!(f2.neg(f2.one()), f2.one());
    }
}
