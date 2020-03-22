use rand::{rngs::ThreadRng, Rng};

use super::{Field, FiniteField};

const CARD: u32 = 7;

const EXP: [F7Elt; CARD as usize] = [1, 3, 2, 6, 4, 5, 1];

const LOG: [u32; CARD as usize] = [CARD as u32, 0, 2, 1, 4, 5, 3];

type F7Elt = u32;

#[derive(Eq)]
pub struct F7 {}

impl PartialEq for F7 {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl Field for F7 {
    type FElt = F7Elt;

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        CARD
    }

    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        (a + b) % CARD
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        (CARD + a - b) % CARD
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        let modulo = |x| {
            if x >= CARD {
                x - (CARD - 1)
            } else {
                x
            }
        };

        if a == 0 || b == 0 {
            0
        } else {
            EXP[modulo(LOG[a as usize] + LOG[b as usize]) as usize]
        }
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        (CARD - a) % CARD
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            Some(EXP[(CARD - 1 - LOG[a as usize]) as usize])
        }
    }

    fn random_element(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, CARD)
    }
}

impl FiniteField for F7 {
    fn characteristic_exponent(&self) -> u32 {
        1
    }

    fn exp(&self, n: u32) -> Self::FElt {
        EXP[n as usize]
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            Some(LOG[a as usize])
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn f7_add() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let c = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.add(a, f.add(b, c)), f.add(f.add(a, b), c));
        assert_eq!(f.add(a, b), f.add(b, a));
        assert_eq!(f.add(a, z), a);
    }

    #[test]
    fn f7_characteristic() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let z = f.zero();
        let mut s = f.zero();
        for _i in 0..7 {
            s = f.add(s, a);
        }

        assert_eq!(s, z);
    }

    #[test]
    fn f7_sub() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.sub(a, z), a);
        assert_eq!(f.sub(a, a), z);
    }

    #[test]
    fn f7_mul() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let c = f.random_element(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.mul(a, f.mul(b, c)), f.mul(f.mul(a, b), c));
        assert_eq!(f.mul(a, b), f.mul(b, a));
        assert_eq!(f.mul(a, i), a);
        assert_eq!(f.mul(a, z), z);
        assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
        assert_eq!(f.mul(a, f.sub(b, c)), f.sub(f.mul(a, b), f.mul(a, c)));
    }

    #[test]
    fn f7_neg() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.neg(z), z);
        assert_eq!(f.neg(f.neg(a)), a);
        assert_eq!(f.add(a, f.neg(b)), f.sub(a, b));
    }

    #[test]
    fn f7_inv() {
        let f = &F7 {};
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.inv(z), None);
        assert_eq!(f.inv(i), Some(i));
        if a != f.zero() {
            assert_eq!(f.inv(f.inv(a).unwrap()), Some(a));
            assert_eq!(f.mul(a, f.inv(a).unwrap()), i);
        }
    }
}
