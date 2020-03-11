use rand::{rngs::ThreadRng, Rng};
// use rand::{
//     distributions::{Distribution, Standard},
//     rngs::ThreadRng,
//     Rng,
// };
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::{Field, FiniteField};
// use crate::finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, F2};

pub use rowvec::RowVec;

// Type T must represent an element from a field, meaning all elements except 0 are inversible.
#[derive(Eq, PartialEq)]
pub struct Mat<'a, F: Eq + Field> {
    field: &'a F,
    rows: usize,
    cols: usize,
    data: Vec<F::FElt>,
}

impl<'a, F: Eq + Field> Clone for Mat<'a, F> {
    fn clone(&self) -> Self {
        Mat {
            field: self.field, // shallow copy
            rows: self.rows,
            cols: self.cols,
            data: self.data.clone(),
        }
    }
}

impl<'a, F: Eq + Field> Add for Mat<'a, F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<'a, F: Eq + Field> Add<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<'a, F: Eq + Field> Add<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn add(self, other: Mat<'a, F>) -> Self::Output {
        self + &other
    }
}

impl<'a, F: Eq + Field> Add for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn add(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        let mut sum = Mat::zero(self.field, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            sum.data[i] = sum.field.add(self.data[i], other.data[i]);
        }
        sum
    }
}

impl<'a, F: Eq + Field> AddAssign<Mat<'a, F>> for Mat<'a, F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<'a, F: Eq + Field> AddAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn add_assign(&mut self, other: &Self) {
        if self.field != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = self.field.add(self.data[i], other.data[i]);
        }
    }
}

impl<'a, F: Eq + Field> Sub for Mat<'a, F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<'a, F: Eq + Field> Sub<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<'a, F: Eq + Field> Sub<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn sub(self, other: Mat<'a, F>) -> Self::Output {
        self - &other
    }
}

impl<'a, F: Eq + Field> Sub for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn sub(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        let mut diff = Mat::zero(self.field, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            diff.data[i] = diff.field.sub(self.data[i], other.data[i]);
        }
        diff
    }
}

impl<'a, F: Eq + Field> SubAssign<Mat<'a, F>> for Mat<'a, F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<'a, F: Eq + Field> SubAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn sub_assign(&mut self, other: &Self) {
        if self.field != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = self.field.sub(self.data[i], other.data[i]);
        }
    }
}

impl<'a, F: Eq + Field> Mul for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<'a, F: Eq + Field> Mul<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, other: Mat<'a, F>) -> Self::Output {
        self * &other
    }
}

impl<'a, F: Eq + Field> Mul for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let mut prod = Mat::zero(self.field, self.rows, other.cols);
        for i in 0..prod.rows {
            for j in 0..prod.cols {
                for k in 0..self.cols {
                    prod[(i, j)] = prod
                        .field
                        .add(prod[(i, j)], prod.field.mul(self[(i, k)], other[(k, j)]));
                }
            }
        }
        prod
    }
}

impl<'a, F: Eq + Field> MulAssign<Mat<'a, F>> for Mat<'a, F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<'a, F: Eq + Field> MulAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn mul_assign(&mut self, other: &Self) {
        if self.field != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows || self.cols != other.cols {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = self.field.zero();
                for k in 0..tmp.cols {
                    self[(i, j)] = self
                        .field
                        .add(self[(i, j)], self.field.mul(tmp[(i, k)], other[(k, j)]));
                }
            }
        }
    }
}

impl<'a, F: Eq + Field> Neg for Mat<'a, F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, F: Eq + Field> Neg for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn neg(self) -> Self::Output {
        let mut opp = Mat::zero(self.field, self.rows, self.cols);

        for i in 0..self.rows * self.cols {
            opp.data[i] = self.field.neg(self.data[i]);
        }
        opp
    }
}

impl<'a, F: Eq + Field> Index<(usize, usize)> for Mat<'a, F> {
    type Output = F::FElt;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

impl<'a, F: Eq + Field> IndexMut<(usize, usize)> for Mat<'a, F> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}

impl<'a, F: Eq + FiniteField> Debug for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(
                    f,
                    "{:>width$}",
                    self.field.to_string_debug(self[(i, j)]),
                    width = 1 + self.field.characteristic_exponent() as usize
                )?;
            }
            write!(
                f,
                "{:>width$}\n",
                self.field.to_string_debug(self[(i, self.cols - 1)]),
                width = 1 + self.field.characteristic_exponent() as usize
            )?;
        }
        Ok(())
    }
}

impl<'a, F: Eq + FiniteField> Display for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        // number of digits of order
        let digits = ((32 - self.field.order().leading_zeros()) / 3 + 1) as usize;

        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(
                    f,
                    "{:>width$} ",
                    self.field.to_string_display(self[(i, j)]),
                    width = if self.field.order() == 2 {
                        1
                    } else {
                        2 + digits
                    }
                )?;
            }
            write!(
                f,
                "{:>width$}\n",
                self.field.to_string_display(self[(i, self.cols - 1)]),
                width = if self.field.order() == 2 {
                    1
                } else {
                    2 + digits
                }
            )?;
        }
        Ok(())
    }
}

impl<'a, F: Eq + Field> Mat<'a, F> {
    pub fn new(field: &'a F, rows: usize, cols: usize, data: Vec<F::FElt>) -> Self {
        if data.len() != rows * cols {
            panic!("Wrong dimensions");
        }
        Self {
            field,
            rows,
            cols,
            data,
        }
    }

    pub fn zero(field: &'a F, rows: usize, cols: usize) -> Self {
        Self {
            field,
            rows,
            cols,
            data: vec![field.zero(); rows * cols],
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn data(&self) -> Vec<F::FElt> {
        self.data.clone()
    }

    //     // pub fn set(&mut self, row: usize, col: usize, val: T) {
    //     //     self[(row, col)] = val;
    //     // }

    //     // pub fn get(&self, row: usize, col: usize) -> T {
    //     //     self[(row, col)]
    //     // }

    pub fn random(rng: &mut ThreadRng, f: &'a F, n: usize, m: usize) -> Self
// where
    //     Standard: Distribution<F::FElt>,
    {
        let mut mat = Self::zero(f, n, m);
        for i in 0..n {
            for j in 0..m {
                // mat[(i, j)] = rng.gen();
                mat[(i, j)] = f.random(rng);
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
                if self[(i, j)] == self.field.zero() {
                    continue;
                } else if self[(i, j)] == self.field.one() {
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

    pub fn permutation_random(rng: &mut ThreadRng, f: &'a F, n: usize) -> Self {
        let mut mat = Self::zero(f, n, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..n {
            // Draw a random column index j to set the '1' on this row
            let nbr = rng.gen_range(0, n - i);
            mat[(i, cols[nbr])] = f.one();

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        mat
    }

    pub fn identity(f: &'a F, n: usize) -> Self {
        let mut id = Self::zero(f, n, n);
        for i in 0..n {
            id[(i, i)] = f.one();
        }
        id
    }

    pub fn sum(&mut self, mat1: &Self, mat2: &Self) {
        if self.field != mat1.field || self.field != mat2.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != mat1.rows
            || self.rows != mat2.rows
            || self.cols != mat1.cols
            || self.cols != mat2.cols
        {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = self.field.add(mat1[(i, j)], mat2[(i, j)]);
            }
        }
    }

    pub fn prod(&mut self, mat1: &Self, mat2: &Self) {
        if self.field != mat1.field || self.field != mat2.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = self.field.zero();
                for k in 0..mat1.cols {
                    sum = self
                        .field
                        .add(sum, self.field.mul(mat1[(i, k)], mat2[(k, j)]));
                }
                self[(i, j)] = sum;
            }
        }
    }

    //     // Generates a random row vector of length n and weight t
    //     pub fn weighted_vector_random(rng: &mut ThreadRng, n: usize, t: usize) -> Self
    //     where
    //         Standard: Distribution<T>,
    //     {
    //         let mut vec = Self::zero(1, n);
    //         let mut cols = Vec::with_capacity(n); // remaining column indices
    //         for i in 0..n {
    //             cols.push(i);
    //         }

    //         for i in 0..t {
    //             // Draw a random column index
    //             let nbr = rng.gen_range(0, n - i);
    //             let mut elt = rng.gen();
    //             while elt == T::zero() {
    //                 elt = rng.gen();
    //             }
    //             vec[(0, cols[nbr])] = elt;

    //             // Remove the index from the list by putting it at the end
    //             cols.swap(nbr, n - 1 - i);
    //         }
    //         vec
    //     }

    //     pub fn weight(&self) -> Option<usize> {
    //         if self.rows != 1 {
    //             return None;
    //         }

    //         let mut cnt = 0;
    //         for j in 0..self.cols {
    //             if self[(0, j)] != T::zero() {
    //                 cnt += 1;
    //             }
    //         }
    //         Some(cnt)
    //     }

    pub fn transpose(&self) -> Self {
        let mut t = Self::zero(self.field, self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }
}

// impl<T: CharacteristicTwo> Mat<'a, F> {
//     pub fn from(a: &Mat<F2>) -> Self {
//         let mut b = Mat::zero(a.rows(), a.cols());
//         for i in 0..a.rows() {
//             for j in 0..a.cols() {
//                 b[(i, j)] = CharacteristicTwo::from(a[(i, j)]);
//             }
//         }
//         b
//     }
// }

// impl<T: CharacteristicTwo + FiniteFieldElement> Mat<'a, F> {
//     pub fn binary_form(&self) -> Mat<F2> {
//         let m = T::characteristic_exponent();
//         let mut bin = Mat::zero(m as usize * self.rows, self.cols);
//         for j in 0..self.cols {
//             for i in 0..self.rows {
//                 for k in 0..m as usize {
//                     bin[(m as usize * i + k, j)] =
//                         match (self[(i, j)].to_canonical_basis() >> k) & 1 {
//                             0 => F2::zero(),
//                             1 => F2::one(),
//                             _ => panic!("Unexpected value"),
//                         }
//                 }
//             }
//         }
//         bin
//     }
// }

mod gauss;
mod rowvec;
