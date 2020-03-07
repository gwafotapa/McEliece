use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, Inv};

use rand::distributions::{Distribution, Standard};
use rand::Rng;
use std::fmt::{Debug, Display, Formatter, Result};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Eq, PartialEq)]
pub enum F2 {
    Zero,
    One,
}

impl CharacteristicTwo for F2 {}

impl FiniteFieldElement for F2 {
    fn characteristic_exponent() -> u32 {
        1
    }

    fn exp(_i: u32) -> Self {
        Self::One
    }

    fn log(self) -> Option<u32> {
        match self {
            Self::Zero => None,
            Self::One => Some(0),
        }
    }

    fn to_canonical_basis(self) -> u32 {
        match self {
            Self::Zero => 0,
            Self::One => 1,
        }
    }
}

impl FieldElement for F2 {
    fn zero() -> Self {
        Self::Zero
    }

    fn one() -> Self {
        Self::One
    }

    fn characteristic() -> u32 {
        2
    }
}

impl Add for F2 {
    type Output = Self;

    /// Adds two finite field elements
    /// ```
    /// # use mceliece::finite_field::FieldElement;
    /// # use mceliece::finite_field::F2;
    ///
    /// assert_eq!(F2::zero() + F2::zero(), F2::zero());
    /// assert_eq!(F2::zero() + F2::one(), F2::one());
    /// assert_eq!(F2::one() + F2::zero(), F2::one());
    /// assert_eq!(F2::one() + F2::one(), F2::zero());
    /// ```

    fn add(self, other: Self) -> Self {
        match self {
            Self::Zero => other,
            Self::One => match other {
                Self::Zero => Self::One,
                Self::One => Self::Zero,
            },
        }
    }
}

impl AddAssign for F2 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for F2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + other
    }
}

impl SubAssign for F2 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Mul for F2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::One => other,
        }
    }
}

impl MulAssign for F2 {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Neg for F2 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Inv for F2 {
    type Output = Self;

    fn inv(self) -> Option<Self::Output> {
        match self {
            Self::Zero => None,
            Self::One => Some(Self::One),
        }
    }
}

impl Distribution<F2> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F2 {
        match rng.gen_bool(0.5) {
            false => F2::Zero,
            true => F2::One,
        }
    }
}

impl Debug for F2 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", *self as u32)
    }
}

impl Display for F2 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", *self as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f2_add() {
        assert_eq!(F2::zero() + F2::zero(), F2::zero());
        assert_eq!(F2::zero() + F2::one(), F2::one());
        assert_eq!(F2::one() + F2::zero(), F2::one());
        assert_eq!(F2::one() + F2::one(), F2::zero());
    }

    #[test]
    fn f2_sub() {
        assert_eq!(F2::zero() - F2::zero(), F2::zero());
        assert_eq!(F2::zero() - F2::one(), F2::one());
        assert_eq!(F2::one() - F2::zero(), F2::one());
        assert_eq!(F2::one() - F2::one(), F2::zero());
    }

    #[test]
    fn f2_mul() {
        assert_eq!(F2::zero() * F2::zero(), F2::zero());
        assert_eq!(F2::zero() * F2::one(), F2::zero());
        assert_eq!(F2::one() * F2::zero(), F2::zero());
        assert_eq!(F2::one() * F2::one(), F2::one());
    }

    #[test]
    fn f2_inv() {
        assert_eq!(F2::zero().inv(), None);
        assert_eq!(F2::one().inv(), Some(F2::one()));
    }

    #[test]
    fn f2_neg() {
        assert_eq!(-F2::zero(), F2::zero());
        assert_eq!(-F2::one(), F2::one());
    }
}
