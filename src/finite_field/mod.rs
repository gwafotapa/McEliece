use rand::{
    distributions::{Distribution, Standard},
    Rng,
};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub use finite_field_2::F2;

mod finite_field_2;

pub trait FiniteFieldElement: FieldElement {
    fn characteristic_exponent() -> u32;
    fn order() -> u32 {
        Self::characteristic().pow(Self::characteristic_exponent())
    }
    fn exp(n: u32) -> Self;
    fn log(self) -> Option<u32>;
    fn to_canonical_basis(self) -> u32;
}

pub trait FieldElement:
    Add<Output = Self>
    + AddAssign
    + Copy
    + Eq
    + Inv<Output = Self>
    + Mul<Output = Self>
    + MulAssign
    + Neg<Output = Self>
    + Sub<Output = Self>
    + SubAssign
{
    fn zero() -> Self;
    fn one() -> Self;
    fn characteristic() -> u32;
}

pub trait Inv {
    type Output;

    fn inv(self) -> Option<Self::Output>;
}

pub trait CharacteristicTwo: FieldElement {
    fn from(elt: F2) -> Self {
        match elt {
            F2::Zero => Self::zero(),
            F2::One => Self::one(),
        }
    }
}

static mut CARD: u32 = 0;
static mut EXP: Vec<FFElt> = Vec::new();
static mut LOG: Vec<u32> = Vec::new();

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct FFElt(pub u32);

impl CharacteristicTwo for FFElt {}

impl FiniteFieldElement for FFElt {
    fn characteristic_exponent() -> u32 {
        3
    }

    fn exp(i: u32) -> Self {
        unsafe { EXP[i as usize] }
    }

    fn log(self) -> Option<u32> {
        if self == Self(0) {
            None
        } else {
            unsafe { Some(LOG[self.0 as usize]) }
        }
    }

    fn to_canonical_basis(self) -> u32 {
        self.0
    }
}

impl FieldElement for FFElt {
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

impl Add for FFElt {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for FFElt {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for FFElt {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + other
    }
}

impl SubAssign for FFElt {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Mul for FFElt {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulo = |a| unsafe {
            if a >= CARD {
                a - (CARD - 1)
            } else {
                a
            }
        };

        if self == Self(0) || other == Self(0) {
            Self(0)
        } else {
            unsafe { EXP[modulo(LOG[self.0 as usize] + LOG[other.0 as usize]) as usize] }
        }
    }
}

impl MulAssign for FFElt {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Neg for FFElt {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Inv for FFElt {
    type Output = Self;

    fn inv(self) -> Option<Self::Output> {
        match self {
            Self(0) => None,
            _ => unsafe { Some(EXP[(CARD - 1 - LOG[self.0 as usize]) as usize]) },
        }
    }
}

impl Distribution<FFElt> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FFElt {
        unsafe { FFElt(rng.gen_range(0, CARD) as u32) }
    }
}

impl Debug for FFElt {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:b}", self.0)
    }
}

impl Display for FFElt {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self(0) => write!(f, "0"),
            Self(1) => write!(f, "1"),
            Self(2) => write!(f, "a"),
            _ => write!(f, "a^{}", self.log().unwrap()),
        }
    }
}

impl FFElt {
    fn primitive_poly(m: u32) -> u32 {
        match m {
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
        }
    }

    pub fn generate_finite_field(order: u32) -> () {
        unsafe {
            if CARD == order {
                return;
            }
        }
        let mut prime_factors = trial_division(order);
        let m = prime_factors.len() as u32;
        prime_factors.dedup();
        if prime_factors.len() != 1 {
            panic!("Order must be a power of a prime");
        }
        let p = prime_factors[0];
        if p != 2 {
            panic!("Only characteristic 2 is supported");
        }

        unsafe {
            CARD = order;
            EXP.resize(order as usize, Self(0));
            LOG.resize(order as usize, 0);
            EXP[0] = Self(1);
            LOG[1] = 0;
            let mut elt = 1_u32;
            for i in 1..order {
                elt *= 2;
                if elt >= order {
                    elt ^= Self::primitive_poly(m);
                }
                EXP[i as usize] = Self(elt);
                LOG[elt as usize] = i;
            }
        }
    }
}

/// Computes the prime factors of (non zero) integer n by trial division
/// https://en.wikipedia.org/wiki/Trial_division
/// ```
/// # use mceliece::finite_field::trial_division;
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
    fn ffelt_add() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let b: FFElt = rng.gen();
        let c: FFElt = rng.gen();
        let z = FFElt::zero();
        assert_eq!(a + (b + c), (a + b) + c);
        assert_eq!(a + b, b + a);
        assert_eq!(a + z, a);
    }

    #[test]
    fn ffelt_characteristic() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let z = FFElt::zero();
        assert_eq!(a + a, z);
    }

    #[test]
    fn ffelt_sub() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let b: FFElt = rng.gen();
        assert_eq!(a + b, a - b);
    }

    #[test]
    fn ffelt_mul() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let b: FFElt = rng.gen();
        let c: FFElt = rng.gen();
        let i = FFElt::one();
        let z = FFElt::zero();
        assert_eq!(a * (b * c), (a * b) * c);
        assert_eq!(a * b, b * a);
        assert_eq!(a * i, a);
        assert_eq!(a * z, z);
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a * (b - c), a * b - a * c);
    }

    #[test]
    fn ffelt_neg() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let b: FFElt = rng.gen();
        let z = FFElt::zero();
        assert_eq!(-z, z);
        assert_eq!(--a, a);
        assert_eq!(a + -b, a - b);
    }

    #[test]
    fn ffelt_inv() {
        FFElt::generate_finite_field(256);

        let mut rng = rand::thread_rng();
        let a: FFElt = rng.gen();
        let i = FFElt::one();
        let z = FFElt::zero();
        assert_eq!(z.inv(), None);
        assert_eq!(i.inv(), Some(i));
        if a != FFElt::zero() {
            FFElt::generate_finite_field(256);
            assert_eq!(a.inv().unwrap().inv().unwrap(), a);
            assert_eq!(a * a.inv().unwrap(), i);
        }
    }
}
