use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// pub use self::{finite_field_2::F2, finite_field_7::F7};
pub use self::{
    finite_field_1024::F1024, finite_field_2::F2, finite_field_256::F256, finite_field_7::F7,
    finite_field_8::F8,
};

mod finite_field_1024;
mod finite_field_2;
mod finite_field_256;
mod finite_field_7;
mod finite_field_8;

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
