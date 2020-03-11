// use rand::{rngs::ThreadRng, Rng};
use rand::rngs::ThreadRng;
// use rand::{
//     distributions::{Distribution, Standard},
//     rngs::ThreadRng,
//     Rng,
// };
use std::fmt::{Formatter, Result};
//     ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
// };

pub use ff2::{F2Elt, F2};
pub use ff2m::F2m;
pub use ff7::F7;

pub mod ff2;
pub mod ff2m;
pub mod ff7;

// pub trait FiniteField: Field {
//     fn characteristic_exponent() -> u32;
//     fn order() -> u32 {
//         Self::characteristic().pow(Self::characteristic_exponent())
//     }
//     // fn exp(n: u32) -> Self;
//     // fn log(self) -> Option<u32>;
//     // fn to_canonical_basis(self) -> u32;
// }

// type F2mElt = u32; // Finite Field of characteristic 2 Element

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
    fn random(&self, rng: &mut ThreadRng) -> Self::FElt;
    fn to_string_debug(&self, a: Self::FElt) -> String;
    fn to_string_display(&self, a: Self::FElt) -> String;
}

pub trait FiniteField: Field {
    fn characteristic_exponent(&self) -> u32;
    fn order(&self) -> u32 {
        self.characteristic().pow(self.characteristic_exponent())
    }
    fn exp(&self, n: u32) -> Self::FElt;
    fn log(&self, a: Self::FElt) -> Option<u32>;
}

// pub trait CharacteristicTwo: FieldElement {
//     fn from(elt: F2) -> Self {
//         match elt {
//             F2::Zero => Self::zero(),
//             F2::One => Self::one(),
//         }
//     }
// }
