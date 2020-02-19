use crate::finite_field;

use rand::distributions;
use rand::Rng;
use std::fmt;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F7(usize);

const CARD: usize = 7;

const EXP: [F7; CARD] = [F7(1), F7(3), F7(2), F7(6), F7(4), F7(5), F7(1)];

const LOG: [usize; CARD] = [CARD, 0, 2, 1, 4, 5, 3];

fn modulo(a: usize) -> usize {
    if a >= CARD {
        a - (CARD - 1)
    } else {
        a
    }
}

pub fn log(a: F7) -> usize {
    LOG[a.0 as usize]
}

pub fn exp(i: usize) -> F7 {
    EXP[i]
}

impl finite_field::Zero for F7 {
    fn zero() -> F7 {
        F7(0)
    }
}

impl finite_field::One for F7 {
    fn one() -> F7 {
        F7(1)
    }
}

impl finite_field::Inverse for F7 {
    fn inv(&self) -> Option<F7> {
        match self {
            F7(0) => None,
            _ => Some(exp(CARD - 1 - log(*self))),
        }
    }
}

impl Add for F7 {
    type Output = F7;

    fn add(self, other: F7) -> F7 {
        F7((self.0 + other.0) % CARD)
    }
}

impl Sub for F7 {
    type Output = F7;

    fn sub(self, other: F7) -> F7 {
        F7(
            ((((self.0 as isize - other.0 as isize) % CARD as isize) + CARD as isize)
                % CARD as isize) as usize,
        )
    }
}

impl AddAssign for F7 {
    fn add_assign(&mut self, other: F7) {
        *self = F7((self.0 + other.0) % CARD);
    }
}

impl Mul for F7 {
    type Output = F7;

    fn mul(self, other: F7) -> F7 {
        if self == F7(0) || other == F7(0) {
            F7(0)
        } else {
            exp(modulo(log(self) + log(other)))
        }
    }
}

impl MulAssign for F7 {
    fn mul_assign(&mut self, other: F7) {
        if *self == F7(0) || other == F7(0) {
            *self = F7(0);
        } else {
            *self = exp(modulo(log(*self) + log(other)))
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
        F7(rng.gen_range(0, CARD) as usize)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::finite_field::Inverse;
    use crate::finite_field::One;
    use crate::finite_field::Zero;

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
}
