use crate::finite_field;

use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, Neg, SubAssign};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F2 {
    Zero,
    One,
}

impl finite_field::Zero for F2 {
    fn zero() -> F2 {
        F2::Zero
    }
}

impl finite_field::One for F2 {
    fn one() -> F2 {
        F2::One
    }
}

impl finite_field::Exp for F2 {
    fn exp(_i: usize) -> Self {
        F2::One
    }
}

impl finite_field::Log for F2 {
    fn log(self) -> Option<usize> {
        match self {
            F2::Zero => None,
            F2::One => Some(0),
        }
    }
}

impl finite_field::AsU32 for F2 {
    fn as_u32(self) -> u32 {
        match self {
            F2::Zero => 0,
            F2::One => 1,
        }
    }
}

impl Neg for F2 {
    type Output = F2;
    
    fn neg(self) -> F2 {
        self
    }
}

impl finite_field::Inv for F2 {
    fn inv(self) -> Option<F2> {
        match self {
            F2::Zero => None,
            F2::One => Some(F2::One),
        }
    }
}

impl Add for F2 {
    type Output = F2;

    fn add(self, other: F2) -> F2 {
        match self {
            F2::Zero => other,
            F2::One => match other {
                F2::Zero => F2::One,
                F2::One => F2::Zero,
            },
        }
    }
}

impl Sub for F2 {
    type Output = F2;

    fn sub(self, other: F2) -> F2 {
        match self {
            F2::Zero => other,
            F2::One => match other {
                F2::Zero => F2::One,
                F2::One => F2::Zero,
            },
        }
    }
}

impl AddAssign for F2 {
    fn add_assign(&mut self, other: F2) {
        *self = match self {
            F2::Zero => other,
            F2::One => match other {
                F2::Zero => F2::One,
                F2::One => F2::Zero,
            },
        };
    }
}

impl SubAssign for F2 {
    fn sub_assign(&mut self, other: F2) {
        // *self = match self {
        //     F2::Zero => other,
        //     F2::One => match other {
        //         F2::Zero => F2::One,
        //         F2::One => F2::Zero,
        //     },
        // };
        *self += other;
    }
}

impl Mul for F2 {
    type Output = F2;

    fn mul(self, other: F2) -> F2 {
        match self {
            F2::Zero => F2::Zero,
            F2::One => other,
        }
    }
}

impl MulAssign for F2 {
    fn mul_assign(&mut self, other: F2) {
        *self = match self {
            F2::Zero => F2::Zero,
            F2::One => other,
        };
    }
}

impl fmt::Display for F2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                F2::Zero => 0,
                F2::One => 1,
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
    use crate::finite_field::Inv;
    use crate::finite_field::One;
    use crate::finite_field::Zero;

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
