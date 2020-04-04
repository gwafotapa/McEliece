use std::{
    cmp,
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    rc::Rc,
};

use super::Poly;
use crate::finite_field::{F2FiniteExtension, Field, FiniteField};

impl<F: Eq + Field> PartialEq for Poly<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.field != other.field || self.degree() != other.degree() {
            return false;
        }
        for i in 0..self.degree() + 1 {
            if self[i] != other[i] {
                return false;
            }
        }
        true
    }
}

impl<F: Eq + Field> Clone for Poly<F> {
    fn clone(&self) -> Self {
        Poly {
            field: Rc::clone(&self.field),
            data: self.data.clone(),
        }
    }
}

impl<F: Eq + Field> Add for Poly<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<F: Eq + Field> Add<&Poly<F>> for Poly<F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<F: Eq + Field> Add<Poly<F>> for &Poly<F> {
    type Output = Poly<F>;

    fn add(self, other: Poly<F>) -> Self::Output {
        self + &other
    }
}

impl<F: Eq + Field> Add for &Poly<F> {
    type Output = Poly<F>;

    fn add(self, other: Self) -> Self::Output {
        let mut sum = self.clone();
        sum += other;
        sum
    }
}

impl<F: Eq + Field> AddAssign<Poly<F>> for Poly<F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<F: Eq + Field> AddAssign<&Poly<F>> for Poly<F> {
    fn add_assign(&mut self, other: &Self) {
        self.data.resize(
            1 + cmp::max(self.degree(), other.degree()),
            self.field.zero(),
        );

        for i in 0..other.degree() + 1 {
            self[i] = self.field.add(self[i], other[i]);
        }

        self.update_len();
    }
}

impl<F: Eq + Field> Sub for Poly<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<F: Eq + Field> Sub<&Poly<F>> for Poly<F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<F: Eq + Field> Sub<Poly<F>> for &Poly<F> {
    type Output = Poly<F>;

    fn sub(self, other: Poly<F>) -> Self::Output {
        self - &other
    }
}

impl<F: Eq + Field> Sub for &Poly<F> {
    type Output = Poly<F>;

    fn sub(self, other: Self) -> Self::Output {
        let mut diff = self.clone();
        diff -= other;
        diff
    }
}

impl<F: Eq + Field> SubAssign<Poly<F>> for Poly<F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<F: Eq + Field> SubAssign<&Poly<F>> for Poly<F> {
    fn sub_assign(&mut self, other: &Self) {
        self.data.resize(
            1 + cmp::max(self.degree(), other.degree()),
            self.field.zero(),
        );

        for i in 0..other.degree() + 1 {
            self[i] = self.field.sub(self[i], other[i]);
        }

        self.update_len();
    }
}

impl<F: Eq + Field> Mul for Poly<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<F: Eq + Field> Mul<&Poly<F>> for Poly<F> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<F: Eq + Field> Mul<Poly<F>> for &Poly<F> {
    type Output = Poly<F>;

    fn mul(self, other: Poly<F>) -> Self::Output {
        self * &other
    }
}

impl<F: Eq + Field> Mul for &Poly<F> {
    type Output = Poly<F>;

    fn mul(self, other: Self) -> Self::Output {
        let f = self.field();
        let mut prod = Poly::zero(f, self.degree() + other.degree() + 1);

        for i in 0..self.degree() + 1 {
            for j in 0..other.degree() + 1 {
                prod[i + j] = f.add(prod[i + j], f.mul(self[i], other[j]));
            }
        }

        prod.update_len();
        prod
    }
}

impl<F: Eq + Field> MulAssign<Poly<F>> for Poly<F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<F: Eq + Field> MulAssign<&Poly<F>> for Poly<F> {
    fn mul_assign(&mut self, other: &Self) {
        let tmp = self.clone();
        let zero = self.field.zero();
        self.data.iter_mut().map(|x| *x = zero).count();
        self.data.resize(tmp.degree() + other.degree() + 1, zero);

        for i in 0..tmp.degree() + 1 {
            for j in 0..other.degree() + 1 {
                self[i + j] = self
                    .field
                    .add(self[i + j], self.field.mul(tmp[i], other[j]));
            }
        }

        self.update_len();
    }
}

impl<F: Eq + Field> Neg for Poly<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<F: Eq + Field> Neg for &Poly<F> {
    type Output = Poly<F>;

    fn neg(self) -> Self::Output {
        let f = self.field();
        let mut opp = Poly::zero(f, self.degree() + 1);

        for i in 0..self.degree() + 1 {
            opp[i] = f.neg(self[i]);
        }
        opp
    }
}

impl<F: Eq + Field> Index<usize> for Poly<F> {
    type Output = F::FieldElement;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<F: Eq + Field> IndexMut<usize> for Poly<F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<F: Eq + F2FiniteExtension> Debug for Poly<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let k = self.field();
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != k.zero() {
                if self[i] != k.one() || i == 0 {
                    write!(f, "{:X}", k.elt_to_u32(self[i]))?;
                }
                match i {
                    0 => (),
                    1 => write!(f, "x")?,
                    _ => write!(f, "x^{}", i)?,
                }
                if i < self.degree() {
                    write!(f, " + ")?;
                }
            }
            i += 1;
        }
        Ok(())
    }
}

impl<F: Eq + FiniteField> Display for Poly<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let k = self.field();
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != k.zero() {
                if self[i] != k.one() || i == 0 {
                    write!(f, "{}", k.elt_to_str(self[i]))?;
                }
                match i {
                    0 => (),
                    1 => write!(f, "x")?,
                    _ => write!(f, "x^{}", i)?,
                }
                if i < self.degree() {
                    write!(f, " + ")?;
                }
            }
            i += 1;
        }
        Ok(())
    }
}
