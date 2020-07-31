use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Neg, Sub, SubAssign},
};

use super::{ColVec, Mat};
use crate::finite_field::{F2FiniteExtension, FieldTrait, FiniteField};

impl<F> From<Mat<F>> for ColVec<F>
where
    F: FieldTrait,
{
    fn from(mat: Mat<F>) -> Self {
        ColVec(mat)
    }
}

impl<F> Clone for ColVec<F>
where
    F: FieldTrait,
{
    fn clone(&self) -> Self {
        ColVec(self.0.clone())
    }
}

impl<F> Add for ColVec<F>
where
    F: FieldTrait,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        self + &other
    }
}

impl<F> Add<&ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        ColVec(self.0 + &other.0)
    }
}

impl<F> Add<ColVec<F>> for &ColVec<F>
where
    F: FieldTrait,
{
    type Output = ColVec<F>;

    fn add(self, other: ColVec<F>) -> Self::Output {
        self + &other
    }
}

impl<F> Add for &ColVec<F>
where
    F: FieldTrait,
{
    type Output = ColVec<F>;

    fn add(self, other: Self) -> Self::Output {
        ColVec(&self.0 + &other.0)
    }
}

impl<F> AddAssign<ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<F> AddAssign<&ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    fn add_assign(&mut self, other: &Self) {
        self.0 += &other.0;
    }
}

impl<F> Sub for ColVec<F>
where
    F: FieldTrait,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        self - &other
    }
}

impl<F> Sub<&ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        ColVec(self.0 - &other.0)
    }
}

impl<F> Sub<ColVec<F>> for &ColVec<F>
where
    F: FieldTrait,
{
    type Output = ColVec<F>;

    fn sub(self, other: ColVec<F>) -> Self::Output {
        self - &other
    }
}

impl<F> Sub for &ColVec<F>
where
    F: FieldTrait,
{
    type Output = ColVec<F>;

    fn sub(self, other: Self) -> Self::Output {
        ColVec(&self.0 - &other.0)
    }
}

impl<F> SubAssign<ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<F> SubAssign<&ColVec<F>> for ColVec<F>
where
    F: FieldTrait,
{
    fn sub_assign(&mut self, other: &Self) {
        self.0 -= &other.0;
    }
}

impl<F> Neg for ColVec<F>
where
    F: FieldTrait,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        ColVec(-self.0)
    }
}

impl<F> Neg for &ColVec<F>
where
    F: FieldTrait,
{
    type Output = ColVec<F>;

    fn neg(self) -> Self::Output {
        ColVec(-&self.0)
    }
}

impl<F> Index<usize> for ColVec<F>
where
    F: FieldTrait,
{
    type Output = F::FieldElement;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[(index, 0)]
    }
}

impl<F> IndexMut<usize> for ColVec<F>
where
    F: FieldTrait,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[(index, 0)]
    }
}

impl<F> Debug for ColVec<F>
where
    F: F2FiniteExtension,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}", self.0)
    }
}

impl<F> Display for ColVec<F>
where
    F: FiniteField,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", self.0)
    }
}
