use std::{
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    rc::Rc,
};

use super::{FieldTrait, FiniteField, Mat, Perm};
use crate::finite_field::{F2FiniteExtension, Field}; // TODO: make some order here and in the crate

impl<F> Clone for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn clone(&self) -> Self {
        Mat {
            field: Rc::clone(&self.field),
            rows: self.rows,
            cols: self.cols,
            data: self.data.clone(),
        }
    }
}

impl<F> Add for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<F> Add<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<F> Add<Mat<F>> for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn add(self, other: Mat<F>) -> Self::Output {
        self + &other
    }
}

impl<F> Add for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn add(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        let f = self.field();
        let mut sum = Mat::zero(Field::Some(f), self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            sum.data[i] = f.add(self.data[i], other.data[i]);
        }
        sum
    }
}

impl<F> AddAssign<Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<F> AddAssign<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
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

impl<F> Sub for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<F> Sub<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<F> Sub<Mat<F>> for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn sub(self, other: Mat<F>) -> Self::Output {
        self - &other
    }
}

impl<F> Sub for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn sub(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        let f = self.field();
        let mut diff = Mat::zero(Field::Some(f), self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            diff.data[i] = f.sub(self.data[i], other.data[i]);
        }
        diff
    }
}

impl<F> SubAssign<Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<F> SubAssign<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
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

impl<F> Mul for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<F> Mul<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<F> Mul<Mat<F>> for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn mul(self, other: Mat<F>) -> Self::Output {
        self * &other
    }
}

impl<F> Mul for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn mul(self, other: Self) -> Self::Output {
        if self.field != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let f = self.field();
        let mut prod = Mat::zero(Field::Some(f), self.rows, other.cols);
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

impl<F> MulAssign<Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<F> MulAssign<&Mat<F>> for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn mul_assign(&mut self, other: &Self) {
        if self.field != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows || self.cols != other.cols {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        let f = tmp.field();
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

impl<F> Mul<Perm> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn mul(self, other: Perm) -> Self::Output {
        &self * &other
    }
}

impl<F> Mul<&Perm> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn mul(self, other: &Perm) -> Self::Output {
        &self * other
    }
}

impl<F> Mul<Perm> for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn mul(self, other: Perm) -> Self::Output {
        self * &other
    }
}

impl<F> Mul<&Perm> for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn mul(self, perm: &Perm) -> Self::Output {
        if self.cols() != perm.len() {
            panic!("Cannot multiply matrices: dimensions don't match");
        }
        self.extract_cols(perm.data())
    }
}

impl<F> Neg for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F> Neg for &Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = Mat<F>;

    fn neg(self) -> Self::Output {
        let f = self.field();
        let mut opp = Mat::zero(Field::Some(f), self.rows, self.cols);

        for i in 0..self.rows * self.cols {
            opp.data[i] = f.neg(self.data[i]);
        }
        opp
    }
}

impl<F> Index<(usize, usize)> for Mat<F>
where
    F: Eq + FieldTrait,
{
    type Output = F::FieldElement;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

impl<F> IndexMut<(usize, usize)> for Mat<F>
where
    F: Eq + FieldTrait,
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}

impl<F> Debug for Mat<F>
where
    F: Eq + F2FiniteExtension,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let k = self.field();
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

impl<F> Display for Mat<F>
where
    F: Eq + FiniteField,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let k = self.field();

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
