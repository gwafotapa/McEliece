// extern crate rand;

use crate::finite_field::{Inverse, One, Zero};
use rand::{distributions, Rng};
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
use std::{cmp, fmt};

#[derive(Clone, Eq, PartialEq)]
pub struct Poly<T>(Vec<T>);

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

impl<T> fmt::Debug for Poly<T>
where
    T: Copy
        + fmt::Display
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Neg<Output = T>
        + AddAssign
        + SubAssign
        + MulAssign
        + Inverse,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\n{}", self.to_str())
    }
}

impl<T> Poly<T>
where
    T: Copy
        + fmt::Display
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Neg<Output = T>
        + AddAssign
        + SubAssign
        + MulAssign
        + Inverse,
{
    pub fn new(n: usize) -> Poly<T> {
        let mut v = Vec::with_capacity(n);
        v.resize(n, T::zero());
        Poly(v)
    }

    pub fn x_n(n: usize) -> Poly<T> {
        let mut v = Vec::with_capacity(n);
        v.resize(n, T::zero());
        v.push(T::one());
        Poly(v)
    }

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn get(&self, i: usize) -> T {
        // self.0[i]
        self[i]
    }

    pub fn set(&mut self, i: usize, val: T) {
        self[i] = val;
    }

    pub fn to_str(&self) -> String {
        if self.degree() == 0 && self[0] == T::zero() {
            return String::from("0");
        }

        let mut s = String::new();
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != T::zero() {
                if self[i] != T::one() || i == 0 {
                    s.push_str(&self[i].to_string());
                }
                match i {
                    0 => (),
                    1 => s.push_str("x"),
                    _ => {
                        s.push_str("x^");
                        s.push_str(&i.to_string());
                    }
                }
                if i < self.degree() {
                    s.push_str(" + ");
                }
            }
            i += 1;
        }

        s
    }

    pub fn random(rng: &mut rand::rngs::ThreadRng, degree: usize) -> Poly<T>
    where
        distributions::Standard: distributions::Distribution<T>,
    {
        let mut p = Poly::new(degree + 1);
        for i in 0..degree {
            p[i] = rng.gen();
        }
        while p[degree] == T::zero() {
            p[degree] = rng.gen();
        }
        p
    }

    pub fn eval(&self, point: T) -> T {
        let mut eval = self.get(0);
        for i in 1..self.degree() + 1 {
            let mut x = T::one();
            for _j in 0..i {
                x *= point;
            }
            eval += self.get(i) * x;
        }

        eval
    }

    pub fn sum(p: &Poly<T>, q: &Poly<T>) -> Poly<T> {
        let mut s = p.clone();
        s.add(q);
        s
    }

    // pub fn as_sum(&mut self, p: &Poly<T>, q: &Poly<T>) {
    //     self.0
    //         .resize(1 + cmp::max(p.degree(), q.degree()), T::zero());

    //     for i in 0..self.degree() + 1 {
    //         self[i] = p[i] + q[i];
    //     }

    //     while self.degree() != 0 && self[self.degree()] == T::zero() {
    //         self.0.pop();
    //     }
    // }

    pub fn add(&mut self, p: &Poly<T>) {
        self.0
            .resize(1 + cmp::max(self.degree(), p.degree()), T::zero());

        for i in 0..p.degree() + 1 {
            self[i] += p[i];
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    pub fn diff(p: &Poly<T>, q: &Poly<T>) -> Poly<T> {
        let mut d = p.clone();
        d.sub(q);
        d
    }

    pub fn sub(&mut self, p: &Poly<T>) {
        self.0
            .resize(1 + cmp::max(self.degree(), p.degree()), T::zero());

        for i in 0..p.degree() + 1 {
            self[i] -= p[i];
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    pub fn prod(a: &Poly<T>, b: &Poly<T>) -> Poly<T> {
        let mut c = Poly::new(a.degree() + b.degree() + 1);

        // for i in 0..c.degree() + 1 {
        //     for j in 0..i + 1 {
        //         c[i] += a[j] * b[i - j];
        //     }
        // }

        for i in 0..a.degree() + 1 {
            for j in 0..b.degree() + 1 {
                c[i + j] += a[i] * b[j];
            }
        }

        while c.degree() != 0 && c[c.degree()] == T::zero() {
            c.0.pop();
        }
        c
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    // Case b = 0 ? Add a result type ?
    pub fn euclidean_division(a: &Poly<T>, b: &Poly<T>) -> (Poly<T>, Poly<T>) {
        let mut q = Poly::new(1);
        let mut r = a.clone();
        let d = b.degree();
        let c = b[d];
        while r.degree() >= d && (r.degree() != 0 || r[0] != T::zero()) {
            let mut s = Poly::x_n(r.degree() - d);
            s[r.degree() - d] = r[r.degree()] * c.inv().unwrap();
            q.add(&s);
            r.sub(&Poly::prod(&s, &b));
        }
        (q, r)
    }

    pub fn gcd(a: &Poly<T>, b: &Poly<T>) -> Poly<T> {
        if b.degree() == 0 && b[b.degree()] == T::zero() {
            return a.clone();
        }

        let (_q, r) = Poly::euclidean_division(a, b);
        Poly::gcd(b, &r)
    }

    pub fn neg(&mut self) {
        for i in 0..self.degree() + 1 {
            self[i] = -self[i];
        }
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn extended_gcd(a: &Poly<T>, b: &Poly<T>) -> (Poly<T>, Poly<T>, Poly<T>, Poly<T>, Poly<T>) {
        let mut r: Vec<Poly<T>> = Vec::new();
        let mut s: Vec<Poly<T>> = Vec::new();
        let mut t: Vec<Poly<T>> = Vec::new();
        r.push(a.clone());
        r.push(b.clone());
        let s0 = Poly::x_n(0);
        let s1 = Poly::new(1);
        s.push(s0);
        s.push(s1);
        let t0 = Poly::new(1);
        let t1 = Poly::x_n(0);
        t.push(t0);
        t.push(t1);
        let mut i = 1;
        while r[i].degree() != 0 || r[i][0] != T::zero() {
            let (q, _r) = Poly::euclidean_division(&r[i - 1], &r[i]);
            // println!("{}\n", q.to_str());
            r.push(Poly::diff(&r[i - 1], &Poly::prod(&q, &r[i])));
            s.push(Poly::diff(&s[i - 1], &Poly::prod(&q, &s[i])));
            t.push(Poly::diff(&t[i - 1], &Poly::prod(&q, &t[i])));
            i += 1;
            // println!("{}\n", r[i].to_str());
        }

        let mut a1 = t.remove(i); // add (-1)^(i-1)
        let mut b1 = s.remove(i); // add (-1)^i
        if i % 2 == 0 {
            a1.neg();
        } else {
            b1.neg();
        }
        let u = s.remove(i - 1);
        let v = t.remove(i - 1);
        let g = r.remove(i - 1);

        (g, u, v, a1, b1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_field_2::F2;

    #[test]
    fn print() {
        let mut p: Poly<F2> = Poly::new(6);
        p[0] = F2::one();
        p[1] = F2::zero();
        p[2] = F2::zero();
        p[3] = F2::zero();
        p[4] = F2::one();
        p[5] = F2::one();
        println!("{}", p.to_str());

        let mut q: Poly<F2> = Poly::x_n(3);
        q[2] = F2::one();
        println!("{}", q.to_str());
    }
}
