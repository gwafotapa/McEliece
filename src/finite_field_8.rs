use crate::finite_field;
use finite_field::Exp;
use finite_field::Log;

use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, Neg, SubAssign};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F8(pub usize);

const CARD: usize = 8;

const EXP: [F8; CARD] = [F8(1), F8(2), F8(4), F8(3), F8(6), F8(7), F8(5), F8(1)];

const LOG: [usize; CARD] = [CARD, 0, 1, 3, 2, 6, 4, 5];

fn modulo(a: usize) -> usize {
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

impl finite_field::Zero for F8 {
    fn zero() -> F8 {
        F8(0)
    }
}

impl finite_field::One for F8 {
    fn one() -> F8 {
        F8(1)
    }
}

impl finite_field::Exp for F8 {
    fn exp(i: usize) -> F8 {
        EXP[i]
    }
}

impl finite_field::Log for F8 {
    fn log(self) -> Option<usize> {
        if self.0 == 0 {
            None
        } else {
            Some(LOG[self.0 as usize])
        }
    }
}

impl Neg for F8 {
    type Output = F8;
    
    fn neg(self) -> F8 {
        self
    }
}

impl finite_field::Inv for F8 {
    fn inv(self) -> Option<F8> {
        match self {
            F8(0) => None,
            _ => Some(F8::exp(CARD - 1 - F8::log(self).unwrap())),
        }
    }
}

impl Add for F8 {
    type Output = F8;

    fn add(self, other: F8) -> F8 {
        F8(self.0 ^ other.0)
    }
}

impl Sub for F8 {
    type Output = F8;

    fn sub(self, other: F8) -> F8 {
        F8(self.0 ^ other.0)
    }
}

impl AddAssign for F8 {
    fn add_assign(&mut self, other: F8) {
        *self = F8(self.0 ^ other.0);
    }
}

impl SubAssign for F8 {
    fn sub_assign(&mut self, other: F8) {
        *self += other;
    }
}

impl Mul for F8 {
    type Output = F8;

    fn mul(self, other: F8) -> F8 {
        if self == F8(0) || other == F8(0) {
            F8(0)
        } else {
            F8::exp(modulo(F8::log(self).unwrap() + F8::log(other).unwrap()))
        }
    }
}

impl MulAssign for F8 {
    fn mul_assign(&mut self, other: F8) {
        if *self == F8(0) || other == F8(0) {
            *self = F8(0);
        } else {
            *self = F8::exp(modulo(F8::log(*self).unwrap() + F8::log(other).unwrap()))
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
        F8(rng.gen_range(0, CARD) as usize)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::finite_field::Inv;
    use crate::finite_field::One;
    use crate::finite_field::Zero;

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
