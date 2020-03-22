use rand::rngs::ThreadRng;

pub use f2::F2;
pub use f2m::F2m;
pub use f7::F7;

pub mod f2;
pub mod f2m;
pub mod f7;

pub trait Field {
    type FElt: Copy + Eq; // Field Element

    fn zero(&self) -> Self::FElt;
    fn one(&self) -> Self::FElt;
    fn characteristic(&self) -> u32;
    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;
    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;
    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt;
    fn neg(&self, a: Self::FElt) -> Self::FElt;
    fn inv(&self, a: Self::FElt) -> Option<Self::FElt>;
    fn random_element(&self, rng: &mut ThreadRng) -> Self::FElt;
}

pub trait FiniteField: Field {
    fn characteristic_exponent(&self) -> u32;
    fn order(&self) -> u32 {
        self.characteristic().pow(self.characteristic_exponent())
    }
    fn exp(&self, n: u32) -> Self::FElt;
    fn log(&self, a: Self::FElt) -> Option<u32>;
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

pub trait CharacteristicTwo: Field {
    // "Converts" the 0 and 1 of F2 to the 0 and 1 of an extension of F2
    fn from(&self, f2: &F2, elt: <F2 as Field>::FElt) -> Self::FElt {
        if elt == f2.zero() {
            self.zero()
        } else {
            self.one()
        }
    }
}

pub trait F2FiniteExtension: CharacteristicTwo + FiniteField {
    // The u32 binary representation matches the decomposition of 'a' on the canonical basis
    fn elt_to_u32(&self, a: Self::FElt) -> u32;
    fn u32_to_elt(&self, n: u32) -> Self::FElt;
}
