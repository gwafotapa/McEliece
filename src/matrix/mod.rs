use rand::{
    distributions::{Distribution, Standard},
    rngs::ThreadRng,
    Rng,
};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    string::ToString,
};

use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, F2};

pub use rowvec::RowVec;

// Type T must represent an element from a field, meaning all elements except 0 are inversible.
#[derive(Clone, Eq, PartialEq)]
pub struct Mat<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T: FieldElement> Add for Mat<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<T: FieldElement> Add<&Mat<T>> for Mat<T> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<T: FieldElement> Add<Mat<T>> for &Mat<T> {
    type Output = Mat<T>;

    fn add(self, other: Mat<T>) -> Self::Output {
        self + &other
    }
}

impl<T: FieldElement> Add for &Mat<T> {
    type Output = Mat<T>;

    fn add(self, other: Self) -> Self::Output {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        let mut sum = Mat::zero(self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            sum.data[i] = self.data[i] + other.data[i];
        }
        sum
    }
}

impl<T: FieldElement> AddAssign<Mat<T>> for Mat<T> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<T: FieldElement> AddAssign<&Mat<T>> for Mat<T> {
    fn add_assign(&mut self, other: &Self) {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] += other.data[i];
        }
    }
}

impl<T: FieldElement> Sub for Mat<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<T: FieldElement> Sub<&Mat<T>> for Mat<T> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<T: FieldElement> Sub<Mat<T>> for &Mat<T> {
    type Output = Mat<T>;

    fn sub(self, other: Mat<T>) -> Self::Output {
        self - &other
    }
}

impl<T: FieldElement> Sub for &Mat<T> {
    type Output = Mat<T>;

    fn sub(self, other: Self) -> Self::Output {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        let mut diff = Mat::zero(self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            diff.data[i] = self.data[i] - other.data[i];
        }
        diff
    }
}

impl<T: FieldElement> SubAssign<Mat<T>> for Mat<T> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<T: FieldElement> SubAssign<&Mat<T>> for Mat<T> {
    fn sub_assign(&mut self, other: &Self) {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] -= other.data[i];
        }
    }
}

impl<T: FieldElement> Mul for Mat<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<T: FieldElement> Mul<&Mat<T>> for Mat<T> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<T: FieldElement> Mul<Mat<T>> for &Mat<T> {
    type Output = Mat<T>;

    fn mul(self, other: Mat<T>) -> Self::Output {
        self * &other
    }
}

impl<T: FieldElement> Mul for &Mat<T> {
    type Output = Mat<T>;

    fn mul(self, other: Self) -> Self::Output {
        if self.cols != other.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let mut prod = Mat::zero(self.rows, other.cols);
        for i in 0..prod.rows {
            for j in 0..prod.cols {
                for k in 0..self.cols {
                    prod[(i, j)] += self[(i, k)] * other[(k, j)];
                }
            }
        }
        prod
    }
}

impl<T: FieldElement> MulAssign<Mat<T>> for Mat<T> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<T: FieldElement> MulAssign<&Mat<T>> for Mat<T> {
    fn mul_assign(&mut self, other: &Self) {
        if self.cols != other.rows || self.cols != other.cols {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = T::zero();
                for k in 0..tmp.cols {
                    self[(i, j)] += tmp[(i, k)] * other[(k, j)];
                }
            }
        }
    }
}

impl<T: FieldElement> Neg for Mat<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<T: FieldElement> Neg for &Mat<T> {
    type Output = Mat<T>;

    fn neg(self) -> Self::Output {
        let mut opp = Mat::zero(self.rows, self.cols);

        for i in 0..self.rows * self.cols {
            opp.data[i] = -self.data[i];
        }
        opp
    }
}

impl<T> Index<(usize, usize)> for Mat<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for Mat<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}

impl<T: ToString> Debug for Mat<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let mut s = String::new();

        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                s.push_str(&self[(i, j)].to_string());
                s.push(' ');
            }
            s.push_str(&self[(i, self.cols - 1)].to_string());
            s.push('\n');
        }

        write!(f, "\n{:?}", s)
    }
}

impl<T: ToString> Display for Mat<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let mut s = String::new();

        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                s.push_str(&self[(i, j)].to_string());
                s.push(' ');
            }
            s.push_str(&self[(i, self.cols - 1)].to_string());
            s.push('\n');
        }

        write!(f, "\n{}", s)
    }
}

impl<T: FieldElement> Mat<T> {
    pub fn new(rows: usize, cols: usize, data: Vec<T>) -> Self {
        if data.len() != rows * cols {
            panic!("Wrong dimensions");
        }
        Mat { rows, cols, data }
    }

    pub fn zero(rows: usize, cols: usize) -> Self {
        Mat {
            rows,
            cols,
            data: vec![T::zero(); rows * cols],
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn data(&self) -> Vec<T> {
        self.data.clone()
    }

    // pub fn set(&mut self, row: usize, col: usize, val: T) {
    //     self[(row, col)] = val;
    // }

    // pub fn get(&self, row: usize, col: usize) -> T {
    //     self[(row, col)]
    // }

    pub fn random(rng: &mut ThreadRng, n: usize, m: usize) -> Self
    where
        Standard: Distribution<T>,
    {
        let mut mat = Mat::zero(n, m);
        for i in 0..n {
            for j in 0..m {
                mat[(i, j)] = rng.gen();
            }
        }
        mat
    }

    pub fn is_permutation(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }

        let n = self.rows;
        let mut cols = Vec::with_capacity(n); // indices of columns containing a '1'
        cols.resize(n, 0);

        for i in 0..n {
            // loop on rows
            let mut has_one = false;
            for j in 0..n {
                if self[(i, j)] == T::zero() {
                    continue;
                } else if self[(i, j)] == T::one() {
                    if cols[j] == 1 {
                        return false;
                    }
                    cols[j] = 1;
                    has_one = true;
                } else {
                    return false;
                }
            }
            if !has_one {
                // row contains zeroes only
                return false;
            }
        }
        true
    }

    pub fn permutation_random(rng: &mut ThreadRng, n: usize) -> Self {
        let mut mat = Mat::zero(n, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..n {
            // Draw a random column index j to set the '1' on this row
            let nbr = rng.gen_range(0, n - i);
            mat[(i, cols[nbr])] = T::one();

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        mat
    }

    pub fn identity(n: usize) -> Self {
        let mut id = Mat::zero(n, n);
        for i in 0..n {
            id[(i, i)] = T::one();
        }
        id
    }

    pub fn sum(&mut self, mat1: &Self, mat2: &Self) {
        if self.rows != mat1.rows
            || self.rows != mat2.rows
            || self.cols != mat1.cols
            || self.cols != mat2.cols
        {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = mat1[(i, j)] + mat2[(i, j)];
            }
        }
    }

    pub fn prod(&mut self, mat1: &Self, mat2: &Self) {
        if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = T::zero();
                for k in 0..mat1.cols {
                    sum += mat1[(i, k)] * mat2[(k, j)];
                }
                self[(i, j)] = sum;
            }
        }
    }

    // Generates a random row vector of length n and weight t
    pub fn weighted_vector_random(rng: &mut ThreadRng, n: usize, t: usize) -> Self
    where
        Standard: Distribution<T>,
    {
        let mut vec = Mat::zero(1, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..t {
            // Draw a random column index
            let nbr = rng.gen_range(0, n - i);
            let mut elt = rng.gen();
            while elt == T::zero() {
                elt = rng.gen();
            }
            vec[(0, cols[nbr])] = elt;

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        vec
    }

    pub fn weight(&self) -> Option<usize> {
        if self.rows != 1 {
            return None;
        }

        let mut cnt = 0;
        for j in 0..self.cols {
            if self[(0, j)] != T::zero() {
                cnt += 1;
            }
        }
        Some(cnt)
    }

    pub fn transpose(&self) -> Self {
        let mut t = Mat::zero(self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }
}

impl<T: CharacteristicTwo> Mat<T> {
    pub fn from(a: &Mat<F2>) -> Self {
        let mut b = Mat::zero(a.rows(), a.cols());
        for i in 0..a.rows() {
            for j in 0..a.cols() {
                b[(i, j)] = CharacteristicTwo::from(a[(i, j)]);
            }
        }
        b
    }
}

impl<T: CharacteristicTwo + FiniteFieldElement> Mat<T> {
    pub fn binary_form(&self) -> Mat<F2> {
        let m = T::characteristic_exponent();
        let mut bin = Mat::zero(m as usize * self.rows, self.cols);
        for j in 0..self.cols {
            for i in 0..self.rows {
                for k in 0..m as usize {
                    bin[(m as usize * i + k, j)] =
                        match (self[(i, j)].to_canonical_basis() >> k) & 1 {
                            0 => F2::zero(),
                            1 => F2::one(),
                            _ => panic!("Unexpected value"),
                        }
                }
            }
        }
        bin
    }
}

mod gauss;
mod rowvec;
