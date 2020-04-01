use std::{
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use super::{Field, FiniteField, Mat, Perm};
use crate::finite_field::F2FiniteExtension;

impl<'a, F: Eq + Field> Clone for Mat<'a, F> {
    fn clone(&self) -> Self {
        Mat {
            field: self.field, // field is shared
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        let mut sum = Mat::zero(f, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            sum.data[i] = f.add(self.data[i], other.data[i]);
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = f.add(self.data[i], other.data[i]);
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        let mut diff = Mat::zero(f, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            diff.data[i] = f.sub(self.data[i], other.data[i]);
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = f.sub(self.data[i], other.data[i]);
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let mut prod = Mat::zero(f, self.rows, other.cols);
        for i in 0..prod.rows {
            for j in 0..prod.cols {
                for k in 0..self.cols {
                    prod[(i, j)] = f.add(prod[(i, j)], f.mul(self[(i, k)], other[(k, j)]));
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
        let f = self.field;
        if f != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows || self.cols != other.cols {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = f.zero();
                for k in 0..tmp.cols {
                    self[(i, j)] = f.add(self[(i, j)], f.mul(tmp[(i, k)], other[(k, j)]));
                }
            }
        }
    }
}

impl<'a, F: Eq + Field> Mul<Perm> for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: Perm) -> Self::Output {
        &self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Perm> for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: &Perm) -> Self::Output {
        &self * other
    }
}

impl<'a, F: Eq + Field> Mul<Perm> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, other: Perm) -> Self::Output {
        self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Perm> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, perm: &Perm) -> Self::Output {
        if self.cols() != perm.len() {
            panic!("Cannot multiply matrices: dimensions don't match");
        }
        self.extract_cols(perm.data())
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
        let f = self.field;
        let mut opp = Mat::zero(f, self.rows, self.cols);

        for i in 0..self.rows * self.cols {
            opp.data[i] = f.neg(self.data[i]);
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

impl<'a, F: Eq + F2FiniteExtension> Debug for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let k = self.field;
        let m = k.characteristic_exponent() as usize;
        let width = crate::div_ceil(m, 4);
        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(f, "{:>w$x} ", k.elt_to_u32(self[(i, j)]), w = width)?;
            }
            write!(
                f,
                "{:>w$x}\n",
                k.elt_to_u32(self[(i, self.cols - 1)]),
                w = width
            )?;
        }
        Ok(())
    }
}

impl<'a, F: Eq + FiniteField> Display for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let k = self.field;

        // Upper bound on the number of digits of order
        let digits =
            ((usize::min_value().leading_zeros() - k.order().leading_zeros()) / 3 + 1) as usize;

        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(
                    f,
                    "{:>width$} ",
                    k.elt_to_str(self[(i, j)]),
                    width = if k.order() == 2 { 1 } else { 2 + digits }
                )?;
            }
            write!(
                f,
                "{:>width$}\n",
                k.elt_to_str(self[(i, self.cols - 1)]),
                width = if k.order() == 2 { 1 } else { 2 + digits }
            )?;
        }
        Ok(())
    }
}
