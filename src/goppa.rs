use crate::finite_field::{Inverse, One, Zero};
use crate::matrix::Mat;
use crate::polynomial::Poly;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

const GOPPA_N: usize = 8;
const GOPPA_T: usize = 2;

pub struct Goppa<T> {
    poly: Poly<T>,
    list: [T; GOPPA_N],
}

impl<T> Goppa<T>
where
    T: Copy
        + fmt::Display
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + AddAssign
        + Sub<Output = T>
        + Mul<Output = T>
        + MulAssign
        + Neg<Output = T>
        + Inverse,
{
    pub fn new(poly: &Poly<T>, list: [T; GOPPA_N]) -> Result<Goppa<T>, &'static str> {
        for i in 0..GOPPA_N {
            if poly.eval(list[i]) == T::zero() {
                println!("{} is a root of the polynomial", list[i]);
                return Err("List contains a root of the polynomial");
            }
        }
        Ok(Goppa::<T> {
            poly: poly.clone(),
            list,
        })
    }

    pub fn poly(&self) -> Poly<T> {
        self.poly.clone()
    }

    pub fn list(&self) -> [T; GOPPA_N] {
        self.list.clone()
    }

    pub fn parity_check_y(&self) -> Mat<T> {
        let mut y = Mat::new(GOPPA_T, GOPPA_N);
        for i in 0..GOPPA_N {
            y.set(0, i, T::one());
        }
        for i in 1..GOPPA_T {
            for j in 0..GOPPA_N {
                let elt = self.list[j];
                y.set(i, j, elt * y.get(i - 1, j));
            }
        }
        y
    }

    pub fn parity_check_z(&self) -> Mat<T> {
        let mut z = Mat::new(GOPPA_N, GOPPA_N);
        for i in 0..GOPPA_N {
            println!("{} {}", self.list[i], self.poly.eval(self.list[i]));
            z.set(i, i, self.poly.eval(self.list[i]).inv().unwrap());
        }
        z
    }

    pub fn parity_check(&self) -> Mat<T> {
        let y = self.parity_check_y();
        let z = self.parity_check_z();
        let mut h = Mat::new(GOPPA_T, GOPPA_N);
        h.mul(&y, &z);
        h
    }

    pub fn generator(&self) -> Mat<T> {
        let h = self.parity_check();
        let m = h.rows();
        let n = h.cols();
        let k = n - m;
        let (_u, hs, p) = h.standard_form().unwrap();
        let mut gs = Mat::new(k, n);
        for i in 0..k {
            gs.set(i, i, T::one());
            for j in k..n {
                gs.set(i, j, -hs.get(j - k, i));
            }
        }
        let mut g = Mat::new(k, n);
        g.mul(&gs, &p.inverse().unwrap().transpose());
        g
    }
}
