use rand::{
    distributions::{Distribution, Standard},
    Rng,
};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::{FieldElement, FiniteFieldElement, Inv};

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct F7(pub u32);

const CARD: u32 = 7;

const EXP: [F7; CARD as usize] = [F7(1), F7(3), F7(2), F7(6), F7(4), F7(5), F7(1)];

const LOG: [u32; CARD as usize] = [CARD as u32, 0, 2, 1, 4, 5, 3];

impl FiniteFieldElement for F7 {
    fn characteristic_exponent() -> u32 {
        1
    }

    fn exp(i: u32) -> Self {
        EXP[i as usize]
    }

    fn log(self) -> Option<u32> {
        if self.0 == 0 {
            None
        } else {
            Some(LOG[self.0 as usize])
        }
    }

    fn to_canonical_basis(self) -> u32 {
        panic!("Method canonical_basis has not been implemented for type F7");
    }
}

impl FieldElement for F7 {
    fn zero() -> Self {
        Self(0)
    }

    fn one() -> Self {
        Self(1)
    }

    fn characteristic() -> u32 {
        7
    }
}

impl Add for F7 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self((self.0 + other.0) % CARD)
    }
}

impl AddAssign for F7 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for F7 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self((CARD + self.0 - other.0) % CARD)
    }
}

impl SubAssign for F7 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Mul for F7 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulo = |a| if a >= CARD { a - (CARD - 1) } else { a };

        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            // Self::exp(modulo(Self::log(self).unwrap() + Self::log(other).unwrap()))
            EXP[modulo(LOG[self.0 as usize] + LOG[other.0 as usize]) as usize]
        }
    }
}

impl MulAssign for F7 {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
        // if *self == Self(0) || other == Self(0) {
        //     *self = Self(0);
        // } else {
        //     *self = Self::exp(modulo(
        //         Self::log(*self).unwrap() + Self::log(other).unwrap(),
        //     ))
        // }
    }
}

impl Neg for F7 {
    type Output = Self;

    fn neg(self) -> Self {
        Self(0) - self
        // match self {
        //     Self(0) => Self(0),
        //     Self(1) => Self(6),
        //     Self(2) => Self(5),
        //     Self(3) => Self(4),
        //     Self(4) => Self(3),
        //     Self(5) => Self(2),
        //     Self(6) => Self(1),
        //     _ => Self(0),
        // }
    }
}

impl Inv for F7 {
    type Output = Self;

    fn inv(self) -> Option<Self::Output> {
        match self {
            Self(0) => None,
            _ => Some(EXP[(CARD - 1 - LOG[self.0 as usize]) as usize]),
        }
    }
}

impl Distribution<F7> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F7 {
        F7(rng.gen_range(0, CARD) as u32)
    }
}

impl Debug for F7 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}", self.0,)
    }
}

impl Display for F7 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        // write!(f, "{}", self.0,)
        match self {
            Self(0) => write!(f, "0"),
            Self(1) => write!(f, "1"),
            Self(2) => write!(f, "a"),
            _ => write!(f, "a^{}", self.log().unwrap()),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn f7_add() {
        let mut rng = rand::thread_rng();
        let a: F7 = rng.gen();
        let b: F7 = rng.gen();
        let c: F7 = rng.gen();
        let z = F7::zero();
        assert_eq!(a + (b + c), (a + b) + c);
        assert_eq!(a + b, b + a);
        assert_eq!(a + z, a);
        assert_eq!(a + a + a + a + a + a + a, z);
    }

    #[test]
    fn f7_sub() {
        let mut rng = rand::thread_rng();
        let a: F7 = rng.gen();
        let z = F7::zero();
        assert_eq!(a - z, a);
        assert_eq!(a - a, z);
    }

    #[test]
    fn f7_mul() {
        let mut rng = rand::thread_rng();
        let a: F7 = rng.gen();
        let b: F7 = rng.gen();
        let c: F7 = rng.gen();
        let i = F7::one();
        let z = F7::zero();
        assert_eq!(a * (b * c), (a * b) * c);
        assert_eq!(a * b, b * a);
        assert_eq!(a * i, a);
        assert_eq!(a * z, z);
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a * (b - c), a * b - a * c);
    }

    #[test]
    fn f7_inv() {
        let mut rng = rand::thread_rng();
        let a: F7 = rng.gen();
        let i = F7::one();
        let z = F7::zero();
        assert_eq!(z.inv(), None);
        assert_eq!(i.inv(), Some(i));
        if a != F7::zero() {
            assert_eq!(a.inv().unwrap().inv().unwrap(), a);
            assert_eq!(a * a.inv().unwrap(), i);
        }
    }

    #[test]
    fn f7_neg() {
        let mut rng = rand::thread_rng();
        let a: F7 = rng.gen();
        let b: F7 = rng.gen();
        let z = F7::zero();
        assert_eq!(-z, z);
        assert_eq!(--a, a);
        assert_eq!(a + -b, a - b);
    }
}
