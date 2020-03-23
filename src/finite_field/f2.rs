use rand::{rngs::ThreadRng, Rng};

use super::{CharacteristicTwo, F2FiniteExtension, Field, FiniteField};

#[derive(Eq, PartialEq)]
pub struct F2 {}

impl Field for F2 {
    type FElt = u32;

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        2
    }

    /// Adds two finite field elements
    /// ```
    /// # use mceliece::finite_field::{Field, F2};
    ///
    /// let f2 = F2 {};
    /// assert_eq!(f2.add(0, 0), 0);
    /// assert_eq!(f2.add(0, 1), 1);
    /// assert_eq!(f2.add(1, 0), 1);
    /// assert_eq!(f2.add(1, 1), 0);
    /// ```    
    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a ^ b
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        self.add(a, b)
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a & b
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        a
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            Some(1)
        }
    }

    fn random_element(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, 2)
    }
}

impl FiniteField for F2 {
    fn characteristic_exponent(&self) -> u32 {
        1
    }

    fn exp(&self, _n: u32) -> Self::FElt {
        1
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            Some(0)
        }
    }

    fn elt_to_str(&self, a: Self::FElt) -> String {
        a.to_string()
    }
}

impl CharacteristicTwo for F2 {
    fn from(&self, _f2: &F2, elt: <F2 as Field>::FElt) -> Self::FElt {
        elt
    }
}

impl F2FiniteExtension for F2 {
    fn elt_to_u32(&self, a: Self::FElt) -> u32 {
        a
    }

    fn u32_to_elt(&self, n: u32) -> Self::FElt {
        n
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f2_add() {
        let f2 = F2 {};
        assert_eq!(f2.add(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.add(f2.zero(), f2.one()), f2.one());
        assert_eq!(f2.add(f2.one(), f2.zero()), f2.one());
        assert_eq!(f2.add(f2.one(), f2.one()), f2.zero());
    }

    #[test]
    fn f2_sub() {
        let f2 = F2 {};
        assert_eq!(f2.sub(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.sub(f2.zero(), f2.one()), f2.one());
        assert_eq!(f2.sub(f2.one(), f2.zero()), f2.one());
        assert_eq!(f2.sub(f2.one(), f2.one()), f2.zero());
    }

    #[test]
    fn f2_mul() {
        let f2 = F2 {};
        assert_eq!(f2.mul(f2.zero(), f2.zero()), f2.zero());
        assert_eq!(f2.mul(f2.zero(), f2.one()), f2.zero());
        assert_eq!(f2.mul(f2.one(), f2.zero()), f2.zero());
        assert_eq!(f2.mul(f2.one(), f2.one()), f2.one());
    }

    #[test]
    fn f2_inv() {
        let f2 = F2 {};
        assert_eq!(f2.inv(f2.zero()), None);
        assert_eq!(f2.inv(f2.one()), Some(f2.one()));
    }

    #[test]
    fn f2_neg() {
        let f2 = F2 {};
        assert_eq!(f2.neg(f2.zero()), f2.zero());
        assert_eq!(f2.neg(f2.one()), f2.one());
    }
}
