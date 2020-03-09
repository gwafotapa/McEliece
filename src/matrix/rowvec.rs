// use crate::finite_field::FieldElement;
use super::{CharacteristicTwo, FieldElement, FiniteFieldElement, Mat, F2};

use rand::{
    distributions::{Distribution, Standard},
    rngs::ThreadRng,
    Rng,
};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Clone, Eq, PartialEq)]
pub struct RowVec<T>(Mat<T>);

impl<T: FieldElement> Add for RowVec<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<T: FieldElement> Add<&RowVec<T>> for RowVec<T> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<T: FieldElement> Add<RowVec<T>> for &RowVec<T> {
    type Output = RowVec<T>;

    fn add(self, other: RowVec<T>) -> Self::Output {
        self + &other
    }
}

impl<T: FieldElement> Add for &RowVec<T> {
    type Output = RowVec<T>;

    fn add(self, other: Self) -> Self::Output {
        RowVec(&self.0 + &other.0)
    }
}

impl<T: FieldElement> AddAssign<RowVec<T>> for RowVec<T> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<T: FieldElement> AddAssign<&RowVec<T>> for RowVec<T> {
    fn add_assign(&mut self, other: &Self) {
        self.0 += &other.0;
    }
}

impl<T: FieldElement> Sub for RowVec<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<T: FieldElement> Sub<&RowVec<T>> for RowVec<T> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<T: FieldElement> Sub<RowVec<T>> for &RowVec<T> {
    type Output = RowVec<T>;

    fn sub(self, other: RowVec<T>) -> Self::Output {
        self - &other
    }
}

impl<T: FieldElement> Sub for &RowVec<T> {
    type Output = RowVec<T>;

    fn sub(self, other: Self) -> Self::Output {
        RowVec(&self.0 - &other.0)
    }
}

impl<T: FieldElement> SubAssign<RowVec<T>> for RowVec<T> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<T: FieldElement> SubAssign<&RowVec<T>> for RowVec<T> {
    fn sub_assign(&mut self, other: &Self) {
        self.0 -= &other.0;
    }
}

impl<T: FieldElement> Mul<Mat<T>> for RowVec<T> {
    type Output = Self;

    fn mul(self, other: Mat<T>) -> Self::Output {
        &self * &other
    }
}

impl<T: FieldElement> Mul<&Mat<T>> for RowVec<T> {
    type Output = Self;

    fn mul(self, other: &Mat<T>) -> Self::Output {
        &self * other
    }
}

impl<T: FieldElement> Mul<Mat<T>> for &RowVec<T> {
    type Output = RowVec<T>;

    fn mul(self, other: Mat<T>) -> Self::Output {
        self * &other
    }
}

impl<T: FieldElement> Mul<&Mat<T>> for &RowVec<T> {
    type Output = RowVec<T>;

    fn mul(self, other: &Mat<T>) -> Self::Output {
        if self.cols() != other.rows() {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let mut prod = RowVec::zero(other.cols());
        for j in 0..prod.cols() {
            prod[j] = T::zero();
            for k in 0..self.cols() {
                prod[j] += self[k] * other[(k, j)];
            }
        }
        prod
    }
}

impl<T: FieldElement> MulAssign<Mat<T>> for RowVec<T> {
    fn mul_assign(&mut self, other: Mat<T>) {
        *self *= &other;
    }
}

impl<T: FieldElement> MulAssign<&Mat<T>> for RowVec<T> {
    fn mul_assign(&mut self, other: &Mat<T>) {
        if self.cols() != other.rows() || self.cols() != other.cols() {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        for j in 0..self.cols() {
            self[j] = T::zero();
            for k in 0..tmp.cols() {
                self[j] += tmp[k] * other[(k, j)];
            }
        }
    }
}

impl<T: FieldElement> Neg for RowVec<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<T: FieldElement> Neg for &RowVec<T> {
    type Output = RowVec<T>;

    fn neg(self) -> Self::Output {
        RowVec(-&self.0)
    }
}

impl<T> Index<usize> for RowVec<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[(0, index)]
    }
}

impl<T> IndexMut<usize> for RowVec<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[(0, index)]
    }
}

impl<T: Debug + Display + FiniteFieldElement> Debug for RowVec<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}", self.0)
    }
}

impl<T: Display + FiniteFieldElement> Display for RowVec<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", self.0)
    }
}

impl<T: FieldElement> RowVec<T> {
    pub fn new(data: Vec<T>) -> Self {
        if data.len() == 0 {
            panic!("Empty row vector");
        }
        RowVec(Mat::new(1, data.len(), data))
    }

    pub fn zero(cols: usize) -> Self {
        RowVec(Mat::zero(1, cols))
    }

    pub fn rows(&self) -> usize {
        1
    }

    pub fn cols(&self) -> usize {
        self.0.cols()
    }

    pub fn data(&self) -> Vec<T> {
        self.0.data()
    }

    pub fn weight(&self) -> usize {
        let mut weight = 0;
        for j in 0..self.cols() {
            if self[j] != T::zero() {
                weight += 1;
            }
        }
        weight
    }

    pub fn random(rng: &mut ThreadRng, n: usize) -> Self
    where Standard: Distribution<T>,
    {
        let mut vec = RowVec::zero(n);
        for i in 0..n {
            vec[i] = rng.gen();
        }
        vec
    }
    
    pub fn random_with_weight(rng: &mut ThreadRng, n: usize, w: usize) -> Self
    where
        Standard: Distribution<T>,
    {
        let mut vec = RowVec::zero(n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..w {
            // Draw a random column index
            let nbr = rng.gen_range(0, n - i);
            let mut elt = rng.gen();
            while elt == T::zero() {
                elt = rng.gen();
            }
            vec[cols[nbr]] = elt;

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        vec
    }

    pub fn sum(&mut self, vec1: &Self, vec2: &Self) {
        if self.rows() != vec1.rows() || self.rows() != vec2.rows() {
            panic!("Cannot add vectors: dimensions don't match");
        }

        self.0.sum(&vec1.0, &vec2.0);
    }

    pub fn from(a: &RowVec<F2>) -> Self
    where
        T: CharacteristicTwo,
    {
        RowVec(Mat::from(&a.0))
    }

    pub fn transpose(&self) -> Mat<T> {
        self.0.transpose()
    }

    // pub fn transpose(&self) -> ColVec<T> {
    //     ColVec::new(self.cols(), self.data())
    // }
}
