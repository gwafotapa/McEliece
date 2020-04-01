//! Finite fields of characteristic 2

use rand::rngs::ThreadRng;

pub use f2::F2;
pub use f2m::F2m;
pub use f7::F7;

pub trait Field {
    /// Field Element
    type FElt: Copy + Eq;

    /// Returns identity element of field addition
    fn zero(&self) -> Self::FElt;

    /// Returns identity element of field multiplication
    fn one(&self) -> Self::FElt;

    /// Returns field characteristic
    fn characteristic(&self) -> usize;

    /// Adds two field elements
    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;

    /// Substracts two field elements
    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;

    /// Multiplies two field elements
    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;

    /// Returns additive inverse of an element
    fn neg(&self, a: Self::FElt) -> Self::FElt;

    /// Returns multiplicative inverse of an element
    fn inv(&self, a: Self::FElt) -> Option<Self::FElt>;

    /// Returns a random element of the field
    fn random_element(&self, rng: &mut ThreadRng) -> Self::FElt;
}

pub trait FiniteField: Field {
    /// Returns m where field order is p<sup>m</sup> with p prime
    fn characteristic_exponent(&self) -> u32;

    /// Returns field order
    fn order(&self) -> usize {
        self.characteristic().pow(self.characteristic_exponent())
    }

    /// Returns the nth power of a primitive element
    fn exp(&self, n: u32) -> Self::FElt;

    /// Returns, if it exists, the discrete logarithm of an element
    fn log(&self, a: Self::FElt) -> Option<u32>;

    /// Converts element to a string as a power of a primitive element
    fn elt_to_str(&self, a: Self::FElt) -> String {
        if a == self.zero() {
            "0".to_owned()
        } else if a == self.one() {
            "1".to_owned()
        } else if self.log(a) == Some(1) {
            "a".to_owned()
        } else {
            format!("a^{}", self.log(a).unwrap())
        }
    }
}

pub trait CharacteristicTwo: Field {}

pub trait F2FiniteExtension: CharacteristicTwo + FiniteField {
    /// Converts an element to an u32
    ///
    /// The conversion is such that the binary representation of the u32
    /// matches the decomposition of the element on the canonical basis.
    fn elt_to_u32(&self, a: Self::FElt) -> u32;

    /// Converts an u32 to an element
    ///
    /// The conversion is such that the binary representation of the u32
    /// matches the decomposition of the element on the canonical basis.
    ///
    /// # Panics
    ///
    /// Panics if the u32 is greater than or equal to field order.
    fn u32_to_elt(&self, n: u32) -> Self::FElt;
}

pub mod f2;
pub mod f2m;
pub mod f7;
