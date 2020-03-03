use crate::finite_field;
// use finite_field::FiniteFieldElement;

use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F2 {
    Zero,
    One,
}

impl finite_field::CharacteristicTwo for F2 {}

impl finite_field::FiniteFieldElement for F2 {
    fn finite_field_q() -> u32 {
        2
    }

    fn finite_field_m() -> u32 {
        1
    }

    fn zero() -> Self {
        Self::Zero
    }

    fn one() -> Self {
        Self::One
    }

    fn inv(self) -> Option<Self> {
        match self {
            Self::Zero => None,
            Self::One => Some(Self::One),
        }
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

    fn to_u32(self) -> u32 {
        match self {
            Self::Zero => 0,
            Self::One => 1,
        }
    }
}

impl Neg for F2 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for F2 {
    type Output = Self;

    /// Adds two finite field elements
    /// ```
    /// use mceliece::finite_field::FiniteFieldElement;
    /// use mceliece::finite_field_2::F2;
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

impl Sub for F2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
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
        *self = match self {
            Self::Zero => other,
            Self::One => match other {
                Self::Zero => Self::One,
                Self::One => Self::Zero,
            },
        };
    }
}

impl SubAssign for F2 {
    fn sub_assign(&mut self, other: Self) {
        // *self = match self {
        //     Self::Zero => other,
        //     Self::One => match other {
        //         Self::Zero => Self::One,
        //         Self::One => Self::Zero,
        //     },
        // };
        *self += other;
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
        *self = match self {
            Self::Zero => Self::Zero,
            Self::One => other,
        };
    }
}

impl fmt::Display for F2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Zero => 0,
                Self::One => 1,
            }
        )
    }
}

impl distributions::Distribution<F2> for distributions::Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F2 {
        match rng.gen_bool(0.5) {
            false => F2::Zero,
            true => F2::One,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_field::FiniteFieldElement;

    #[test]
    fn f2_add() {
        assert_eq!(F2::zero() + F2::zero(), F2::zero());
        assert_eq!(F2::zero() + F2::one(), F2::one());
        assert_eq!(F2::one() + F2::zero(), F2::one());
        assert_eq!(F2::one() + F2::one(), F2::zero());
    }

    #[test]
    fn f2_sub() {
        assert_eq!(F2::zero() + F2::zero(), F2::zero());
        assert_eq!(F2::zero() + F2::one(), F2::one());
        assert_eq!(F2::one() + F2::zero(), F2::one());
        assert_eq!(F2::one() + F2::one(), F2::zero());
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
        let mut rng = rand::thread_rng();
        let a: F2 = rng.gen();
        assert_eq!(-a, a);
    }
}
