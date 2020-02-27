use crate::finite_field::{Inv, One, Zero, AsU32};
// use crate::finite_field_2::F2;
use crate::matrix::Mat;
use crate::polynomial::Poly;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

const M: u32 = 10;

pub struct Goppa<T> {
    len: usize,
    poly: Poly<T>,
    set: Vec<T>,
}

impl<T> Goppa<T>
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
        + Inv
        + AddAssign
        + SubAssign
    + MulAssign
    + AsU32
{
    pub fn new(poly: Poly<T>, set: Vec<T>) -> Result<Goppa<T>, &'static str> {
        for i in 0..set.len() {
            if poly.eval(set[i]) == T::zero() {
                println!("{} is a root of the polynomial", set[i]);
                return Err("Set contains a root of the polynomial");
            }
        }
        Ok(Goppa::<T> {
            len: set.len(),
            poly,
            set,
        })
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn poly(&self) -> &Poly<T> {
        &self.poly
    }

    pub fn set(&self) -> &Vec<T> {
        &self.set
    }

    pub fn parity_check_matrix_y(&self) -> Mat<T> {
        let n = self.len();
        let t = self.poly.degree();

        let mut y = Mat::new(t, n);
        for i in 0..n {
            y[(0, i)] = T::one();
        }
        for i in 1..t {
            for j in 0..n {
                let elt = self.set[j];
                y[(i, j)] = elt * y[(i - 1, j)];
            }
        }
        y
    }

    pub fn parity_check_matrix_z(&self) -> Mat<T> {
        let n = self.len();

        let mut z = Mat::new(n, n);
        for i in 0..n {
            // println!("{} {}", self.set[i], self.poly.eval(self.set[i]));
            z[(i, i)] = self.poly.eval(self.set[i]).inv().unwrap();
        }
        z
    }

    pub fn parity_check_matrix(&self) -> Mat<T> {
        let y = self.parity_check_matrix_y();
        let z = self.parity_check_matrix_z();
        let mut h = Mat::new(self.poly.degree(), self.len);
        h.as_prod(&y, &z);
        h
    }

    pub fn generator_matrix(&self) -> Mat<T> {
        let h = self.parity_check_matrix();
        let m = h.rows();
        let n = h.cols();
        let k = n - m;
        let (_u, hs, p) = h.standard_form().unwrap();
        let mut gs = Mat::new(k, n);
        for i in 0..k {
            gs[(i, i)] = T::one();
            for j in k..n {
                gs[(i, j)] = -hs[(j - k, i)];
            }
        }
        let mut g = Mat::new(k, n);
        g.as_prod(&gs, &p.inverse().unwrap().transpose());
        g
    }

    pub fn encode(&self, msg: &Mat<T>) -> Mat<T> {
        let g = self.generator_matrix();
        let cpt = Mat::prod(msg, &g.transpose());
        cpt
    }

    pub fn decode(&self, rcv: &Mat<T>) -> Mat<T> {
        let h = self.parity_check_matrix();
        let syndrome = Mat::prod(&h, &rcv.transpose());
        let zero = Mat::new(1, syndrome.cols());
        if syndrome == zero {
            return rcv.clone();
        }

        let mut s_x = Poly::new(syndrome.cols());
        for i in 0..syndrome.cols() {
            s_x[i] = syndrome[(1, i)];
        }

        let mut t = s_x.inverse_modulo(M, &self.poly);
        t.add(&Poly::x_n(1));
        t.square_root_modulo(M, &self.poly);
        let (mut a, mut b, _, _, _) = Poly::extended_gcd(&t, &self.poly);
        a.square();
        b.square();
        b.mul(&Poly::x_n(1));
        let n = self.set.len();
        // let aa = Mat::new(1, n);
        // for i in 0..a.degree() {
        //     aa[(1, i)] = a[i];
        // }
        // let bb = Mat::new(1, n);
        // for i in 0..b.degree() {
        //     bb[(1, i)] = b[i];
        // }
        let sigma = Poly::sum(&a, &b);
        let mut err = Mat::new(1, n);
        for i in 0..n { // use iterator and map ?
            err[(1, i)] = if sigma.eval(self.set[i]) == T::zero() { T::one() } else { T::zero() };
        }

        let msg = Mat::sum(&rcv, &err);
        msg
    }
}
