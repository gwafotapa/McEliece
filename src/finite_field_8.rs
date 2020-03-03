use crate::finite_field;
use finite_field::FiniteFieldElement;

use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F8(pub u32);

const CARD: u32 = 8;

const EXP: [F8; CARD as usize] = [F8(1), F8(2), F8(4), F8(3), F8(6), F8(7), F8(5), F8(1)];

const LOG: [u32; CARD as usize] = [CARD, 0, 1, 3, 2, 6, 4, 5];

fn modulo(a: u32) -> u32 {
    if a >= CARD {
        a - (CARD - 1)
    } else {
        a
    }
}

// pub fn log(a: F8) -> usize {
//     LOG[a.0 as usize]
// }

// pub fn exp(i: usize) -> F8 {
//     EXP[i]
// }

impl finite_field::CharacteristicTwo for F8 {}

impl finite_field::FiniteFieldElement for F8 {
    fn finite_field_q() -> u32 {
        2
    }

    fn finite_field_m() -> u32 {
        3
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
        if self.0 == 0 {
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

impl Neg for F8 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for F8 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl Sub for F8 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for F8 {
    fn add_assign(&mut self, other: Self) {
        *self = Self(self.0 ^ other.0);
    }
}

impl SubAssign for F8 {
    fn sub_assign(&mut self, other: Self) {
        *self += other;
    }
}

impl Mul for F8 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            Self::exp(modulo(Self::log(self).unwrap() + Self::log(other).unwrap()))
        }
    }
}

impl MulAssign for F8 {
    fn mul_assign(&mut self, other: Self) {
        if *self == Self(0) || other == Self(0) {
            *self = Self(0);
        } else {
            *self = Self::exp(modulo(
                Self::log(*self).unwrap() + Self::log(other).unwrap(),
            ))
        }
    }
}

impl fmt::Display for F8 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0,)
    }
}

impl distributions::Distribution<F8> for distributions::Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F8 {
        F8(rng.gen_range(0, CARD) as u32)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::finite_field::FiniteFieldElement;

    #[test]
    fn f8_add() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let b: F8 = rng.gen();
        let c: F8 = rng.gen();
        let z = F8::zero();
        assert_eq!(a + (b + c), (a + b) + c);
        assert_eq!(a + b, b + a);
        assert_eq!(a + z, a);
        assert_eq!(a + a, z);
    }

    #[test]
    fn f8_sub() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let z = F8::zero();
        assert_eq!(a - z, a);
        assert_eq!(a - a, z);
    }

    #[test]
    fn f8_mul() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let b: F8 = rng.gen();
        let c: F8 = rng.gen();
        let i = F8::one();
        let z = F8::zero();
        assert_eq!(a * (b * c), (a * b) * c);
        assert_eq!(a * b, b * a);
        assert_eq!(a * i, a);
        assert_eq!(a * z, z);
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a * (b - c), a * b - a * c);
    }

    #[test]
    fn f8_inv() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let i = F8::one();
        let z = F8::zero();
        assert_eq!(z.inv(), None);
        assert_eq!(i.inv(), Some(i));
        if a != F8::zero() {
            assert_eq!(a.inv().unwrap().inv().unwrap(), a);
            assert_eq!(a * a.inv().unwrap(), i);
        }
    }

    #[test]
    fn f8_neg() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        assert_eq!(-a, a);
    }
}
