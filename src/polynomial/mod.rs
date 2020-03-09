// use log::info;

use rand::{
    distributions::{Distribution, Standard},
    rngs::ThreadRng,
    Rng,
};
use std::{
    cmp,
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::FieldElement;

#[derive(Clone, Eq, PartialEq)]
pub struct Poly<T>(Vec<T>);

impl<T: FieldElement> Add for Poly<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<T: FieldElement> Add<&Poly<T>> for Poly<T> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<T: FieldElement> Add<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, other: Poly<T>) -> Self::Output {
        self + &other
    }
}

impl<T: FieldElement> Add for &Poly<T> {
    type Output = Poly<T>;

    fn add(self, other: Self) -> Self::Output {
        let mut sum = self.clone();
        sum += other;
        sum
    }
}

impl<T: FieldElement> AddAssign<Poly<T>> for Poly<T> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<T: FieldElement> AddAssign<&Poly<T>> for Poly<T> {
    fn add_assign(&mut self, other: &Self) {
        self.0
            .resize(1 + cmp::max(self.degree(), other.degree()), T::zero());

        for i in 0..other.degree() + 1 {
            self[i] += other[i];
        }

        self.update_len();
    }
}

impl<T: FieldElement> Sub for Poly<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<T: FieldElement> Sub<&Poly<T>> for Poly<T> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<T: FieldElement> Sub<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, other: Poly<T>) -> Self::Output {
        self - &other
    }
}

impl<T: FieldElement> Sub for &Poly<T> {
    type Output = Poly<T>;

    fn sub(self, other: Self) -> Self::Output {
        let mut diff = self.clone();
        diff -= other;
        diff
    }
}

impl<T: FieldElement> SubAssign<Poly<T>> for Poly<T> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<T: FieldElement> SubAssign<&Poly<T>> for Poly<T> {
    fn sub_assign(&mut self, other: &Self) {
        self.0
            .resize(1 + cmp::max(self.degree(), other.degree()), T::zero());

        for i in 0..other.degree() + 1 {
            self[i] -= other[i];
        }

        self.update_len();
    }
}

impl<T: FieldElement> Mul for Poly<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<T: FieldElement> Mul<&Poly<T>> for Poly<T> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<T: FieldElement> Mul<Poly<T>> for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, other: Poly<T>) -> Self::Output {
        self * &other
    }
}

impl<T: FieldElement> Mul for &Poly<T> {
    type Output = Poly<T>;

    fn mul(self, other: Self) -> Self::Output {
        let mut prod = Poly::zero(self.degree() + other.degree() + 1);

        for i in 0..self.degree() + 1 {
            for j in 0..other.degree() + 1 {
                prod[i + j] += self[i] * other[j];
            }
        }

        prod.update_len();
        prod
    }
}

impl<T: FieldElement> MulAssign<Poly<T>> for Poly<T> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<T: FieldElement> MulAssign<&Poly<T>> for Poly<T> {
    fn mul_assign(&mut self, other: &Self) {
        let tmp = self.clone();
        self.0.iter_mut().map(|x| *x = T::zero()).count();
        self.0.resize(tmp.degree() + other.degree() + 1, T::zero());

        for i in 0..tmp.degree() + 1 {
            for j in 0..other.degree() + 1 {
                self[i + j] += tmp[i] * other[j];
            }
        }

        self.update_len();
    }
}

impl<T: FieldElement> Neg for Poly<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<T: FieldElement> Neg for &Poly<T> {
    type Output = Poly<T>;

    fn neg(self) -> Self::Output {
        let mut opp = Poly::zero(self.degree() + 1);

        for i in 0..self.degree() + 1 {
            opp[i] = -self[i];
        }
        opp
    }
}

impl<T> Index<usize> for Poly<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> IndexMut<usize> for Poly<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T: FieldElement + Debug> Debug for Poly<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != T::zero() {
                if self[i] != T::one() || i == 0 {
                    write!(f, "{:?}", self[i])?;
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

impl<T: FieldElement + Display> Display for Poly<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != T::zero() {
                if self[i] != T::one() || i == 0 {
                    write!(f, "{}", self[i])?;
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

impl<T: FieldElement> Poly<T> {
    pub fn zero(capacity: usize) -> Self {
        if capacity == 0 {
            panic!("Vector capacity cannot be zero");
        }
        Poly(vec![T::zero(); capacity])
    }

    pub fn x_n(n: usize) -> Self {
        let mut v = vec![T::zero(); n + 1];
        v[n] = T::one();
        Poly(v)
    }

    pub fn degree(&self) -> usize {
        let mut degree = self.0.len() - 1;
        while self[degree] == T::zero() && degree != 0 {
            degree -= 1;
        }
        degree
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self[0] == T::zero()
    }

    pub fn update_len(&mut self) {
        self.0.truncate(self.degree() + 1);
    }

    pub fn random(rng: &mut ThreadRng, degree: usize) -> Self
    where
        Standard: Distribution<T>,
    {
        let mut p = Poly::zero(degree + 1);
        for i in 0..degree {
            p[i] = rng.gen();
        }
        while p[degree] == T::zero() {
            p[degree] = rng.gen();
        }
        p
    }

    pub fn eval(&self, point: T) -> T {
        let mut eval = self[0];
        for i in 1..self.degree() + 1 {
            let mut x = T::one();
            for _j in 0..i {
                x *= point;
            }
            eval += self[i] * x;
        }

        eval
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn euclidean_division(a: &Self, b: &Self) -> (Self, Self) {
        if b.is_zero() {
            panic!("Euclidean division by the null polynomial");
        }

        let mut q = Poly::zero(1);
        let mut r = a.clone();
        let d = b.degree();
        let c = b[d];
        while r.degree() >= d && !r.is_zero() {
            let mut s = Poly::x_n(r.degree() - d);
            s[r.degree() - d] = r[r.degree()] * c.inv().unwrap();
            q += &s;
            r -= &s * b;
        }
        (q, r)
    }

    pub fn modulo(&mut self, modulus: &Self) {
        let mut q = Poly::zero(1);
        let d = modulus.degree();
        let c = modulus[d];
        while self.degree() >= d && !self.is_zero() {
            let mut s = Poly::x_n(self.degree() - d);
            s[self.degree() - d] = self[self.degree()] * c.inv().unwrap();
            q += &s;
            *self -= &s * modulus;
        }

        self.update_len();
    }

    pub fn gcd(a: &Self, b: &Self) -> Self {
        if b.is_zero() {
            return a.clone();
        }

        let (_q, r) = Poly::euclidean_division(a, b);
        Poly::gcd(b, &r)
    }

    pub fn neg_mut(&mut self) {
        for i in 0..self.degree() + 1 {
            self[i] = -self[i];
        }
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn extended_gcd(a: &Self, b: &Self) -> (Self, Self, Self, Self, Self) {
        let mut r: Vec<Self> = Vec::new();
        let mut s: Vec<Self> = Vec::new();
        let mut t: Vec<Self> = Vec::new();
        r.push(a.clone());
        r.push(b.clone());
        let s0 = Poly::x_n(0);
        let s1 = Poly::zero(1);
        s.push(s0);
        s.push(s1);
        let t0 = Poly::zero(1);
        let t1 = Poly::x_n(0);
        t.push(t0);
        t.push(t1);
        let mut i = 1;
        while !r[i].is_zero() {
            let (q, _r) = Poly::euclidean_division(&r[i - 1], &r[i]);
            // println!("{}\n", q.to_str());
            r.push(&r[i - 1] - &q * &r[i]);
            s.push(&s[i - 1] - &q * &s[i]);
            t.push(&t[i - 1] - &q * &t[i]);
            i += 1;
            // println!("{}\n", r[i].to_str());
        }

        let mut a1 = t.remove(i);
        let mut b1 = s.remove(i);
        if i % 2 == 0 {
            a1.neg_mut();
        } else {
            b1.neg_mut();
        }
        let u = s.remove(i - 1);
        let v = t.remove(i - 1);
        let g = r.remove(i - 1);

        (g, u, v, a1, b1)
    }

    // TODO: is this algorithm ok outside of characteristic two ?
    // or should the goppa code be binary in which case it needs to be moved in the submodule
    pub fn goppa_extended_gcd(g: &Self, t: &Self) -> (Self, Self) {
        let mut i = 1;

        let mut a: Vec<Self> = Vec::new();
        let mut b: Vec<Self> = Vec::new();
        a.push(g.clone());
        a.push(t.clone());
        let b0 = Poly::zero(1);
        let b1 = Poly::x_n(0);
        b.push(b0);
        b.push(b1);

        while a[i].degree() > g.degree() / 2 {
            i += 1;
            let (q, r) = Poly::euclidean_division(&a[i - 2], &a[i - 1]);
            // println!("{}\n", q.to_str());
            b.push(&b[i - 2] + &q * &b[i - 1]);
            a.push(r);
            // println!("{}\n", r[i].to_str());

            // if (Poly::prod(&b[i], &t).modulo(&g) != a[i].modulo(&g)) {
            //     panic!("Broken loop invariant between a, b, t and g");
            // }
        }

        (a.remove(i), b.remove(i))
    }
}

pub mod characteristic_two;
