use crate::finite_field;
use finite_field::FiniteFieldElement;

use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F7(u32);

const CARD: u32 = 7;

const EXP: [F7; CARD as usize] = [F7(1), F7(3), F7(2), F7(6), F7(4), F7(5), F7(1)];

const LOG: [u32; CARD as usize] = [CARD as u32, 0, 2, 1, 4, 5, 3];

fn modulo(a: u32) -> u32 {
    if a >= CARD {
        a - (CARD - 1)
    } else {
        a
    }
}

// pub fn log(a: F7) -> usize {
//     LOG[a.0 as usize]
// }

// pub fn exp(i: usize) -> F7 {
//     EXP[i]
// }

impl finite_field::FiniteFieldElement for F7 {
    fn finite_field_q() -> u32 {
        7
    }

    fn finite_field_m() -> u32 {
        1
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

impl Neg for F7 {
    type Output = Self;

    fn neg(self) -> Self {
        match self {
            Self(0) => Self(0),
            Self(1) => Self(6),
            Self(2) => Self(5),
            Self(3) => Self(4),
            Self(4) => Self(3),
            Self(5) => Self(2),
            Self(6) => Self(1),
            _ => Self(0),
        }
    }
}

impl Add for F7 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self((self.0 + other.0) % CARD)
    }
}

impl Sub for F7 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(
            ((((self.0 as isize - other.0 as isize) % CARD as isize) + CARD as isize)
                % CARD as isize) as u32,
        )
    }
}

impl AddAssign for F7 {
    fn add_assign(&mut self, other: Self) {
        // *self = Self((self.0 + other.0) % CARD);
        *self = *self + other;
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
        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            Self::exp(modulo(Self::log(self).unwrap() + Self::log(other).unwrap()))
        }
    }
}

impl MulAssign for F7 {
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

impl fmt::Display for F7 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0,)
    }
}

impl distributions::Distribution<F7> for distributions::Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F7 {
        F7(rng.gen_range(0, CARD) as u32)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::finite_field::FiniteFieldElement;

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
