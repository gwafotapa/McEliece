use rand::{rngs::ThreadRng, Rng};
use std::{
    convert::TryInto,
    error::Error,
    fmt::{self, Debug, Display, Formatter},
    fs::File,
    io::{Read, Write},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    rc::Rc,
};

use super::{Mat, Perm};
use crate::finite_field::{F2FiniteExtension, Field, FiniteField, F2};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Eq, PartialEq)]
pub struct RowVec<F: Eq + Field>(Mat<F>);

impl<F: Eq + Field> Clone for RowVec<F> {
    fn clone(&self) -> Self {
        RowVec(self.0.clone())
    }
}

impl<F: Eq + Field> Add for RowVec<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<F: Eq + Field> Add<&RowVec<F>> for RowVec<F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<F: Eq + Field> Add<RowVec<F>> for &RowVec<F> {
    type Output = RowVec<F>;

    fn add(self, other: RowVec<F>) -> Self::Output {
        self + &other
    }
}

impl<F: Eq + Field> Add for &RowVec<F> {
    type Output = RowVec<F>;

    fn add(self, other: Self) -> Self::Output {
        RowVec(&self.0 + &other.0)
    }
}

impl<F: Eq + Field> AddAssign<RowVec<F>> for RowVec<F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<F: Eq + Field> AddAssign<&RowVec<F>> for RowVec<F> {
    fn add_assign(&mut self, other: &Self) {
        self.0 += &other.0;
    }
}

impl<F: Eq + Field> Sub for RowVec<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<F: Eq + Field> Sub<&RowVec<F>> for RowVec<F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<F: Eq + Field> Sub<RowVec<F>> for &RowVec<F> {
    type Output = RowVec<F>;

    fn sub(self, other: RowVec<F>) -> Self::Output {
        self - &other
    }
}

impl<F: Eq + Field> Sub for &RowVec<F> {
    type Output = RowVec<F>;

    fn sub(self, other: Self) -> Self::Output {
        RowVec(&self.0 - &other.0)
    }
}

impl<F: Eq + Field> SubAssign<RowVec<F>> for RowVec<F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<F: Eq + Field> SubAssign<&RowVec<F>> for RowVec<F> {
    fn sub_assign(&mut self, other: &Self) {
        self.0 -= &other.0;
    }
}

impl<F: Eq + Field> Mul<Mat<F>> for RowVec<F> {
    type Output = Self;

    fn mul(self, other: Mat<F>) -> Self::Output {
        &self * &other
    }
}

impl<F: Eq + Field> Mul<&Mat<F>> for RowVec<F> {
    type Output = Self;

    fn mul(self, other: &Mat<F>) -> Self::Output {
        &self * other
    }
}

impl<F: Eq + Field> Mul<Mat<F>> for &RowVec<F> {
    type Output = RowVec<F>;

    fn mul(self, other: Mat<F>) -> Self::Output {
        self * &other
    }
}

impl<F: Eq + Field> Mul<&Mat<F>> for &RowVec<F> {
    type Output = RowVec<F>;

    fn mul(self, other: &Mat<F>) -> Self::Output {
        RowVec(&self.0 * other)
    }
}

impl<F: Eq + Field> MulAssign<Mat<F>> for RowVec<F> {
    fn mul_assign(&mut self, other: Mat<F>) {
        *self *= &other;
    }
}

impl<F: Eq + Field> MulAssign<&Mat<F>> for RowVec<F> {
    fn mul_assign(&mut self, other: &Mat<F>) {
        self.0 *= other;
    }
}

impl<F: Eq + Field> Mul<Perm> for RowVec<F> {
    type Output = Self;

    fn mul(self, other: Perm) -> Self::Output {
        &self * &other
    }
}

impl<F: Eq + Field> Mul<&Perm> for RowVec<F> {
    type Output = Self;

    fn mul(self, other: &Perm) -> Self::Output {
        &self * other
    }
}

impl<F: Eq + Field> Mul<Perm> for &RowVec<F> {
    type Output = RowVec<F>;

    fn mul(self, other: Perm) -> Self::Output {
        self * &other
    }
}

impl<F: Eq + Field> Mul<&Perm> for &RowVec<F> {
    type Output = RowVec<F>;

    fn mul(self, other: &Perm) -> Self::Output {
        self.extract_cols(other.data())
    }
}

impl<F: Eq + Field> Neg for RowVec<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F: Eq + Field> Neg for &RowVec<F> {
    type Output = RowVec<F>;

    fn neg(self) -> Self::Output {
        RowVec(-&self.0)
    }
}

impl<F: Eq + Field> Index<usize> for RowVec<F> {
    type Output = F::FieldElement;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[(0, index)]
    }
}

impl<F: Eq + Field> IndexMut<usize> for RowVec<F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[(0, index)]
    }
}

impl<F: Eq + F2FiniteExtension> Debug for RowVec<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.0)
    }
}

impl<F: Eq + FiniteField> Display for RowVec<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<F: Eq + Field> RowVec<F> {
    pub fn new(field: &Rc<F>, data: Vec<F::FieldElement>) -> Self {
        if data.len() == 0 {
            panic!("Empty row vector");
        }
        RowVec(Mat::new(field, 1, data.len(), data))
    }

    pub fn zero(field: &Rc<F>, cols: usize) -> Self {
        RowVec(Mat::zero(field, 1, cols))
    }

    pub fn field(&self) -> &Rc<F> {
        self.0.field()
    }

    pub fn rows(&self) -> usize {
        1
    }

    pub fn cols(&self) -> usize {
        self.0.cols()
    }

    pub fn data(&self) -> &Vec<F::FieldElement> {
        self.0.data()
    }

    pub fn weight(&self) -> usize {
        let mut weight = 0;
        for j in 0..self.cols() {
            if self[j] != self.0.field.zero() {
                weight += 1;
            }
        }
        weight
    }

    pub fn random(rng: &mut ThreadRng, f: &Rc<F>, n: usize) -> Self {
        RowVec(Mat::random(rng, f, 1, n))
    }

    pub fn random_with_weight(rng: &mut ThreadRng, f: &Rc<F>, n: usize, w: usize) -> Self {
        let mut vec = RowVec::zero(f, n);
        let mut cols = Vec::with_capacity(n);
        for i in 0..n {
            cols.push(i);
        }

        for _i in 0..w {
            let mut elt = f.random_element(rng);
            while elt == f.zero() {
                elt = f.random_element(rng);
            }
            let index = rng.gen_range(0, cols.len());
            vec[cols[index]] = elt;
            cols.swap_remove(index);
        }
        vec
    }

    pub fn transpose(&self) -> Mat<F> {
        self.0.transpose()
    }

    pub fn is_zero(&self) -> bool {
        for i in 0..self.cols() {
            if self[i] != self.field().zero() {
                return false;
            }
        }
        true
    }

    pub fn extract_cols(&self, perm: &Vec<usize>) -> Self {
        RowVec(self.0.extract_cols(perm))
    }
}

impl<'a> RowVec<F2> {
    pub fn random_f2(rng: &mut ThreadRng, n: usize) -> Self {
        let f2 = &Rc::new(F2::generate(()));
        RowVec(Mat::random(rng, f2, 1, n))
    }

    pub fn write(&self, file_name: &str) -> Result<()> {
        let mut f = File::create(file_name)?;
        let len = 4 + crate::div_ceil(self.cols(), 8);
        let mut vec = Vec::with_capacity(len);
        vec.extend_from_slice(&(self.cols() as u32).to_be_bytes());
        let mut byte = 0;
        let mut shift = 7;
        for i in 0..self.cols() {
            byte |= (self[i] as u8) << shift;
            if shift == 0 {
                vec.push(byte);
                byte = 0;
                shift = 7;
            } else {
                shift -= 1;
            }
        }
        if self.cols() % 8 != 0 {
            vec.push(byte);
        }
        f.write_all(&vec)?;
        Ok(())
    }

    pub fn read_vector(file_name: &str) -> Result<RowVec<F2>> {
        let mut f = File::open(file_name)?;
        let mut vec = Vec::new();
        f.read_to_end(&mut vec)?;
        let cols = u32::from_be_bytes(vec[0..4].try_into()?) as usize;
        let f2 = &Rc::new(F2 {});
        let mut rowvec = RowVec::zero(f2, cols);
        let mut k = 4;
        let mut shift = 7;
        for i in 0..cols {
            rowvec[i] = ((vec[k] >> shift) & 1).into();
            if shift == 0 {
                k += 1;
                shift = 7;
            } else {
                shift -= 1;
            }
        }
        Ok(rowvec)
    }
}
