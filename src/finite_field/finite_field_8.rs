use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, Inv};

use rand::{
    distributions::{Distribution, Standard},
    Rng,
};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

macro_rules! array_init {
    ( $( $x:expr ),+ ) => {
        [ $( F8($x) ),+ ]
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct F8(pub u32);

const CARD: u32 = 8;

const EXP: [F8; CARD as usize] = array_init![
    1, 2, 4, 3, 6, 7, 5, 1
];

const LOG: [u32; CARD as usize] = [
    CARD, 0, 1, 3, 2, 6, 4, 5
];

impl CharacteristicTwo for F8 {}

impl FiniteFieldElement for F8 {
    fn characteristic_exponent() -> u32 {
        3
    }

    fn exp(i: u32) -> Self {
        EXP[i as usize]
    }

    fn log(self) -> Option<u32> {
        if self == Self(0) {
            None
        } else {
            Some(LOG[self.0 as usize])
        }
    }

    fn to_canonical_basis(self) -> u32 {
        self.0
    }
}

impl FieldElement for F8 {
    fn zero() -> Self {
        Self(0)
    }

    fn one() -> Self {
        Self(1)
    }

    fn characteristic() -> u32 {
        2
    }
}

impl Add for F8 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for F8 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for F8 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + other
    }
}

impl SubAssign for F8 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Mul for F8 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulo = |a| if a >= CARD { a - (CARD - 1) } else { a };

        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            EXP[modulo(LOG[self.0 as usize] + LOG[other.0 as usize]) as usize]
        }
    }
}

impl MulAssign for F8 {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Neg for F8 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Inv for F8 {
    type Output = Self;

    fn inv(self) -> Option<Self::Output> {
        match self {
            Self(0) => None,
            _ => Some(EXP[(CARD - 1 - LOG[self.0 as usize]) as usize]),
        }
    }
}

impl Distribution<F8> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> F8 {
        F8(rng.gen_range(0, CARD) as u32)
    }
}

impl Debug for F8 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:b}", self.0)
    }
}

impl Display for F8 {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
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
    fn f8_add() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let b: F8 = rng.gen();
        let c: F8 = rng.gen();
        let z = F8::zero();
        assert_eq!(a + (b + c), (a + b) + c);
        assert_eq!(a + b, b + a);
        assert_eq!(a + z, a);
    }

    #[test]
    fn f8_characteristic() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let z = F8::zero();
        assert_eq!(a + a, z);
    }

    #[test]
    fn f8_sub() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let b: F8 = rng.gen();
        assert_eq!(a + b, a - b);
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
    fn f8_neg() {
        let mut rng = rand::thread_rng();
        let a: F8 = rng.gen();
        let b: F8 = rng.gen();
        let z = F8::zero();
        assert_eq!(-z, z);
        assert_eq!(--a, a);
        assert_eq!(a + -b, a - b);
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
}