use rand::{rngs::ThreadRng, Rng};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    fs::File,
    io::{Read, Write},
};

use super::Mat;
use crate::finite_field::{CharacteristicTwo, Field, FiniteField, F2, F2FiniteExtension};

#[derive(Eq, PartialEq)]
pub struct RowVec<'a, F: Eq + Field>(Mat<'a, F>);

impl<'a, F: Eq + Field> Clone for RowVec<'a, F> {
    fn clone(&self) -> Self {
        RowVec(self.0.clone())
    }
}

impl<'a, F: Eq + Field> Add for RowVec<'a, F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<'a, F: Eq + Field> Add<&RowVec<'a, F>> for RowVec<'a, F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<'a, F: Eq + Field> Add<RowVec<'a, F>> for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn add(self, other: RowVec<'a, F>) -> Self::Output {
        self + &other
    }
}

impl<'a, F: Eq + Field> Add for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn add(self, other: Self) -> Self::Output {
        RowVec(&self.0 + &other.0)
    }
}

impl<'a, F: Eq + Field> AddAssign<RowVec<'a, F>> for RowVec<'a, F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<'a, F: Eq + Field> AddAssign<&RowVec<'a, F>> for RowVec<'a, F> {
    fn add_assign(&mut self, other: &Self) {
        self.0 += &other.0;
    }
}

impl<'a, F: Eq + Field> Sub for RowVec<'a, F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<'a, F: Eq + Field> Sub<&RowVec<'a, F>> for RowVec<'a, F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<'a, F: Eq + Field> Sub<RowVec<'a, F>> for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn sub(self, other: RowVec<'a, F>) -> Self::Output {
        self - &other
    }
}

impl<'a, F: Eq + Field> Sub for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn sub(self, other: Self) -> Self::Output {
        RowVec(&self.0 - &other.0)
    }
}

impl<'a, F: Eq + Field> SubAssign<RowVec<'a, F>> for RowVec<'a, F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<'a, F: Eq + Field> SubAssign<&RowVec<'a, F>> for RowVec<'a, F> {
    fn sub_assign(&mut self, other: &Self) {
        self.0 -= &other.0;
    }
}

impl<'a, F: Eq + Field> Mul<Mat<'a, F>> for RowVec<'a, F> {
    type Output = Self;

    fn mul(self, other: Mat<'a, F>) -> Self::Output {
        &self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Mat<'a, F>> for RowVec<'a, F> {
    type Output = Self;

    fn mul(self, other: &Mat<'a, F>) -> Self::Output {
        &self * other
    }
}

impl<'a, F: Eq + Field> Mul<Mat<'a, F>> for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn mul(self, other: Mat<'a, F>) -> Self::Output {
        self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Mat<'a, F>> for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn mul(self, other: &Mat<'a, F>) -> Self::Output {
        RowVec(&self.0 * other)
    }
}

impl<'a, F: Eq + Field> MulAssign<Mat<'a, F>> for RowVec<'a, F> {
    fn mul_assign(&mut self, other: Mat<'a, F>) {
        *self *= &other;
    }
}

impl<'a, F: Eq + Field> MulAssign<&Mat<'a, F>> for RowVec<'a, F> {
    fn mul_assign(&mut self, other: &Mat<'a, F>) {
        self.0 *= other;
    }
}

impl<'a, F: Eq + Field> Neg for RowVec<'a, F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, F: Eq + Field> Neg for &RowVec<'a, F> {
    type Output = RowVec<'a, F>;

    fn neg(self) -> Self::Output {
        RowVec(-&self.0)
    }
}

impl<'a, F: Eq + Field> Index<usize> for RowVec<'a, F> {
    type Output = F::FElt;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[(0, index)]
    }
}

impl<'a, F: Eq + Field> IndexMut<usize> for RowVec<'a, F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[(0, index)]
    }
}

impl<'a, F: Eq + F2FiniteExtension> Debug for RowVec<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}", self.0)
    }
}

impl<'a, F: Eq + FiniteField> Display for RowVec<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", self.0)
    }
}

impl<'a, F: Eq + Field> RowVec<'a, F> {
    pub fn new(field: &'a F, data: Vec<F::FElt>) -> Self {
        if data.len() == 0 {
            panic!("Empty row vector");
        }
        RowVec(Mat::new(field, 1, data.len(), data))
    }

    pub fn zero(field: &'a F, cols: usize) -> Self {
        RowVec(Mat::zero(field, 1, cols))
    }

    pub fn field(&self) -> &'a F {
        self.0.field()
    }

    pub fn rows(&self) -> usize {
        1
    }

    pub fn cols(&self) -> usize {
        self.0.cols()
    }

    pub fn data(&self) -> Vec<F::FElt> {
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

    pub fn random(rng: &mut ThreadRng, f: &'a F, n: usize) -> Self {
        RowVec(Mat::random(rng, f, 1, n))
    }

    pub fn random_with_weight(rng: &mut ThreadRng, f: &'a F, n: usize, w: usize) -> Self {
        let mut vec = RowVec::zero(f, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..w {
            // Draw a random column index
            let nbr = rng.gen_range(0, n - i);
            let mut elt = f.random_element(rng);
            while elt == f.zero() {
                elt = f.random_element(rng);
            }
            vec[cols[nbr]] = elt;

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        vec
    }

    // TODO: Is this function useful ?
    pub fn sum(&mut self, vec1: &Self, vec2: &Self) {
        let f = self.field();
        if f != vec1.field() || f != vec2.field() {
            panic!("Cannot add vectors: fields don't match");
        }
        if self.rows() != vec1.rows() || self.rows() != vec2.rows() {
            panic!("Cannot add vectors: dimensions don't match");
        }
        self.0.sum(&vec1.0, &vec2.0);
    }

    pub fn from<'b>(f: &'a F, vec_f2: &RowVec<'b, F2>) -> Self
    where
        F: CharacteristicTwo,
    {
        RowVec(Mat::from(f, &vec_f2.0))
    }

    pub fn transpose(&self) -> Mat<'a, F> {
        self.0.transpose()
    }

    // TODO: Is this function useful ?
    // pub fn transpose(&self) -> ColVec<T> {
    //     ColVec::new(self.cols(), self.data())
    // }

    pub fn is_zero(&self) -> bool {
        for i in 0..self.cols() {
            if self[i] != self.field().zero() {
                return false;
            }
        }
        true
    }

    pub fn permute_cols(&mut self, cols: &Vec<usize>) {
        self.0.permute_cols(cols);
    }
}

impl<'a> RowVec<'a, F2> {
    pub fn save_vector(&self, file_name: &str) {
        let mut f = File::create(file_name).expect("Unable to create file");
        let mut s = format!("{:x}#", self.cols());
        let mut byte: u8 = 0;
        let mut cnt_mod_8 = 7;
        for i in 0..self.cols() {
            byte |= (self[i] << cnt_mod_8) as u8;
            if cnt_mod_8 == 0 {
                s.push_str(format!("{:02x}", byte).as_str());
                cnt_mod_8 = 7;
                byte = 0;
            } else {
                cnt_mod_8 -= 1;
            }
        }
        if self.cols() % 8 != 0 {
            s.push_str(format!("{:02x}", byte).as_str());
        }
        f.write_all(s.as_bytes());
    }

    pub fn load_vector(file_name: &str, f2: &'a F2) -> RowVec<'a, F2> {
        let mut f = File::open(file_name).expect("Unable to open file");
        let mut buf = String::new();
        f.read_to_string(&mut buf).expect("Unable to read data");
        let v: Vec<&str> = buf.split('#').collect();
        let cols = u32::from_str_radix(v[0], 16).unwrap() as usize;
        let mut vec = RowVec::zero(f2, cols);
        let mut cnt_mod_8 = 7;
        let s = hex::decode(v[1]).unwrap();
        let mut iter = s.iter();
        let mut byte = iter.next().unwrap();
        for i in 0..cols {
            vec[i] = (*byte as u32 >> cnt_mod_8) & 1_u32;
            if cnt_mod_8 == 0 && i != cols - 1 {
                byte = iter.next().unwrap();
                cnt_mod_8 = 7;
            } else {
                cnt_mod_8 -= 1;
            }
        }
        vec
    }
}
