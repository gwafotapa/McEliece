use std::{
    cmp,
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::{F2FiniteExtension, Field, FiniteField};
use super::Poly;

impl<'a, F: Eq + Field> PartialEq for Poly<'a, F> {
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

impl<'a, F: Eq + Field> Clone for Poly<'a, F> {
    fn clone(&self) -> Self {
        Poly {
            field: self.field, // field is shared
            data: self.data.clone(),
        }
    }
}

impl<'a, F: Eq + Field> Add for Poly<'a, F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<'a, F: Eq + Field> Add<&Poly<'a, F>> for Poly<'a, F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<'a, F: Eq + Field> Add<Poly<'a, F>> for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn add(self, other: Poly<'a, F>) -> Self::Output {
        self + &other
    }
}

impl<'a, F: Eq + Field> Add for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn add(self, other: Self) -> Self::Output {
        let mut sum = self.clone();
        sum += other;
        sum
    }
}

impl<'a, F: Eq + Field> AddAssign<Poly<'a, F>> for Poly<'a, F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<'a, F: Eq + Field> AddAssign<&Poly<'a, F>> for Poly<'a, F> {
    fn add_assign(&mut self, other: &Self) {
        let f = self.field;
        self.data
            .resize(1 + cmp::max(self.degree(), other.degree()), f.zero());

        for i in 0..other.degree() + 1 {
            self[i] = f.add(self[i], other[i]);
        }

        self.update_len();
    }
}

impl<'a, F: Eq + Field> Sub for Poly<'a, F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<'a, F: Eq + Field> Sub<&Poly<'a, F>> for Poly<'a, F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<'a, F: Eq + Field> Sub<Poly<'a, F>> for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn sub(self, other: Poly<'a, F>) -> Self::Output {
        self - &other
    }
}

impl<'a, F: Eq + Field> Sub for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn sub(self, other: Self) -> Self::Output {
        let mut diff = self.clone();
        diff -= other;
        diff
    }
}

impl<'a, F: Eq + Field> SubAssign<Poly<'a, F>> for Poly<'a, F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<'a, F: Eq + Field> SubAssign<&Poly<'a, F>> for Poly<'a, F> {
    fn sub_assign(&mut self, other: &Self) {
        let f = self.field;
        self.data
            .resize(1 + cmp::max(self.degree(), other.degree()), f.zero());

        for i in 0..other.degree() + 1 {
            self[i] = f.sub(self[i], other[i]);
        }

        self.update_len();
    }
}

impl<'a, F: Eq + Field> Mul for Poly<'a, F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Poly<'a, F>> for Poly<'a, F> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<'a, F: Eq + Field> Mul<Poly<'a, F>> for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn mul(self, other: Poly<'a, F>) -> Self::Output {
        self * &other
    }
}

impl<'a, F: Eq + Field> Mul for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn mul(self, other: Self) -> Self::Output {
        let f = self.field;
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

impl<'a, F: Eq + Field> MulAssign<Poly<'a, F>> for Poly<'a, F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<'a, F: Eq + Field> MulAssign<&Poly<'a, F>> for Poly<'a, F> {
    fn mul_assign(&mut self, other: &Self) {
        let f = self.field;
        let tmp = self.clone();
        self.data.iter_mut().map(|x| *x = f.zero()).count();
        self.data
            .resize(tmp.degree() + other.degree() + 1, f.zero());

        for i in 0..tmp.degree() + 1 {
            for j in 0..other.degree() + 1 {
                self[i + j] = f.add(self[i + j], f.mul(tmp[i], other[j]));
            }
        }

        self.update_len();
    }
}

impl<'a, F: Eq + Field> Neg for Poly<'a, F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, F: Eq + Field> Neg for &Poly<'a, F> {
    type Output = Poly<'a, F>;

    fn neg(self) -> Self::Output {
        let f = self.field;
        let mut opp = Poly::zero(f, self.degree() + 1);

        for i in 0..self.degree() + 1 {
            opp[i] = f.neg(self[i]);
        }
        opp
    }
}

impl<'a, F: Eq + Field> Index<usize> for Poly<'a, F> {
    type Output = F::FElt;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<'a, F: Eq + Field> IndexMut<usize> for Poly<'a, F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<'a, F: Eq + F2FiniteExtension> Debug for Poly<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let k = self.field;
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

impl<'a, F: Eq + FiniteField> Display for Poly<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let k = self.field;
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
