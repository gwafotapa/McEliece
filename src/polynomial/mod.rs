// use log::info;

use rand::rngs::ThreadRng;
use std::{
    cmp,
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::Field;

// #[derive(Clone, Eq, PartialEq)]
#[derive(Eq, PartialEq)]
pub struct Poly<'a, F: Eq + Field> {
    field: &'a F,
    data: Vec<F::FElt>,
}

impl<'a, F: Eq + Field> Clone for Poly<'a, F> {
    fn clone(&self) -> Self {
        Poly {
            field: self.field, // shallow copy
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
        let mut opp = Poly::zero(self.field, self.degree() + 1);

        for i in 0..self.degree() + 1 {
            opp[i] = self.field.neg(self[i]);
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

impl<'a, F: Eq + Field> Debug for Poly<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != self.field.zero() {
                if self[i] != self.field.one() || i == 0 {
                    write!(f, "{}", self.field.to_string_debug(self[i]))?;
                    // write!(f, "{:?}", self[i])?;
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

impl<'a, F: Eq + Field> Display for Poly<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != self.field.zero() {
                if self[i] != self.field.one() || i == 0 {
                    write!(f, "{}", self.field.to_string_display(self[i]))?;
                    // write!(f, "{}", self[i])?;
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

impl<'a, F: Eq + Field> Poly<'a, F> {
    pub fn zero(field: &'a F, capacity: usize) -> Self {
        if capacity == 0 {
            panic!("Vector capacity cannot be zero");
        }
        Self {
            field,
            data: vec![field.zero(); capacity],
        }
    }

    pub fn x_n(field: &'a F, n: usize) -> Self {
        let mut v = vec![field.zero(); n + 1];
        v[n] = field.one();
        Self { field, data: v }
    }

    pub fn support(f: &'a F, support: &[usize]) -> Self {
        if support.is_empty() {
            panic!("Support cannot be empty");
        }

        let mut vec = vec![f.zero(); support[support.len() - 1]];
        for &i in support.iter() {
            if i >= vec.len() {
                vec.resize(i + 1, f.zero());
            }
            vec[i] = f.one();
        }
        Self {
            field: f,
            data: vec,
        }
    }

    pub fn field(&self) -> &'a F {
        self.field
    }
    
    pub fn degree(&self) -> usize {
        let mut degree = self.data.len() - 1;
        while self[degree] == self.field.zero() && degree != 0 {
            degree -= 1;
        }
        degree
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self[0] == self.field.zero()
    }

    pub fn update_len(&mut self) {
        self.data.truncate(self.degree() + 1);
    }

    pub fn random(rng: &mut ThreadRng, f: &'a F, degree: usize) -> Self {
        let mut p = Self::zero(f, degree + 1);
        for i in 0..degree {
            p[i] = f.random_element(rng);
        }
        while p[degree] == f.zero() {
            p[degree] = f.random_element(rng);
        }
        p
    }

    pub fn eval(&self, point: F::FElt) -> F::FElt {
        let f = self.field;
        let mut eval = self[0];
        for i in 1..self.degree() + 1 {
            let mut x = f.one();
            for _j in 0..i {
                x = f.mul(x, point);
            }
            eval = f.add(eval, f.mul(self[i], x));
        }
        eval
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn euclidean_division(a: &Self, b: &Self) -> (Self, Self) {
        if a.field != b.field {
            panic!("Cannot compute euclidean division: fields differ")
        }
        if b.is_zero() {
            panic!("Euclidean division by the null polynomial");
        }

        let f = a.field;
        let mut q = Self::zero(f, 1);
        let mut r = a.clone();
        let d = b.degree();
        let c = b[d];
        while r.degree() >= d && !r.is_zero() {
            let mut s = Self::x_n(f, r.degree() - d);
            s[r.degree() - d] = f.mul(r[r.degree()], f.inv(c).unwrap());
            q += &s;
            r -= &s * b;
        }
        (q, r)
    }

    pub fn modulo(&mut self, modulus: &Self) {
        let f = self.field;
        let mut q = Self::zero(f, 1);
        let d = modulus.degree();
        let c = modulus[d];
        while self.degree() >= d && !self.is_zero() {
            let mut s = Self::x_n(f, self.degree() - d);
            s[self.degree() - d] = f.mul(self[self.degree()], f.inv(c).unwrap());
            q += &s;
            *self -= &s * modulus;
        }

        self.update_len();
    }

    pub fn gcd(a: &Self, b: &Self) -> Self {
        if a.field != b.field {
            panic!("Cannot compute euclidean division: fields differ")
        }
        if b.is_zero() {
            return a.clone();
        }

        let (_q, r) = Self::euclidean_division(a, b);
        Self::gcd(b, &r)
    }

    pub fn neg_mut(&mut self) {
        for i in 0..self.degree() + 1 {
            self[i] = self.field.neg(self[i]);
        }
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn extended_gcd(a: &Self, b: &Self) -> (Self, Self, Self, Self, Self) {
        if a.field != b.field {
            panic!("Cannot compute euclidean division: fields differ")
        }

        let f = a.field;
        let mut r: Vec<Self> = Vec::new();
        let mut s: Vec<Self> = Vec::new();
        let mut t: Vec<Self> = Vec::new();
        r.push(a.clone());
        r.push(b.clone());
        let s0 = Self::x_n(f, 0);
        let s1 = Self::zero(f, 1);
        s.push(s0);
        s.push(s1);
        let t0 = Self::zero(f, 1);
        let t1 = Self::x_n(f, 0);
        t.push(t0);
        t.push(t1);
        let mut i = 1;
        while !r[i].is_zero() {
            let (q, _r) = Self::euclidean_division(&r[i - 1], &r[i]);
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
        if g.field != t.field {
            panic!("Cannot compute euclidean division: fields differ")
        }

        let f = g.field;
        let mut i = 1;

        let mut a: Vec<Self> = Vec::new();
        let mut b: Vec<Self> = Vec::new();
        a.push(g.clone());
        a.push(t.clone());
        let b0 = Self::zero(f, 1);
        let b1 = Self::x_n(f, 0);
        b.push(b0);
        b.push(b1);

        while a[i].degree() > g.degree() / 2 {
            i += 1;
            let (q, r) = Self::euclidean_division(&a[i - 2], &a[i - 1]);
            // println!("{}\n", q.to_str());
            b.push(&b[i - 2] + &q * &b[i - 1]);
            a.push(r);
            // println!("{}\n", r[i].to_str());

            // if (Self::prod(&b[i], &t).modulo(&g) != a[i].modulo(&g)) {
            //     panic!("Broken loop invariant between a, b, t and g");
            // }
        }

        (a.remove(i), b.remove(i))
    }
}

pub mod characteristic_two;
