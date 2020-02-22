use crate::finite_field::{Inverse, One, Zero};
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub};
use std::{cmp, fmt};

#[derive(Clone, Eq, PartialEq)]
pub struct Poly<T>(Vec<T>);

impl<T> Index<usize> for Poly<T>
where
    T: Copy
        + fmt::Display
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + AddAssign
        + MulAssign
        + Inverse,
{
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> IndexMut<usize> for Poly<T>
where
    T: Copy
        + fmt::Display
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + AddAssign
        + MulAssign
        + Inverse,
{
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
        + AddAssign
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
        + AddAssign
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
        if self.degree() == 0 && self.get(0) == T::zero() {
            return String::from("0");
        }

        let mut s = String::new();
        let mut i = 0;
        while i <= self.degree() {
            if self.get(i) != T::zero() {
                if self.get(i) != T::one() || i == 0 {
                    s.push_str(&self.get(i).to_string());
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

    // TODO
    // Rename add function to sum
    // Add a add function mutating the original polynomial

    // pub fn add(&mut self, p: &Poly<T>, q: &Poly<T>) {
    //     self.0
    //         .resize(1 + cmp::max(p.degree(), q.degree()), T::zero());

    //     for i in 0..self.degree() + 1 {
    //         self.set(i, p.get(i) + q.get(i));
    //     }

    //     while self.get(self.degree()) == T::zero() {
    //         self.0.pop();
    //     }
    // }
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
