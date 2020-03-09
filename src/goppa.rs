use log::info;

use rand::{distributions, Rng};
use std::string::ToString;

use crate::{
    finite_field::{CharacteristicTwo, FieldElement, FiniteFieldElement, F2},
    matrix::{Mat, RowVec},
    polynomial::Poly,
};

pub struct Goppa<T> {
    len: usize,
    poly: Poly<T>,
    set: Vec<T>,
}

impl<T> Goppa<T>
where
    T: CharacteristicTwo + FiniteFieldElement + ToString,
{
    pub fn new(poly: Poly<T>, set: Vec<T>) -> Result<Goppa<T>, &'static str>
    where
        T: std::fmt::Debug,
    {
        if poly.is_irreducible() == false {
            return Err("Goppa polynomial is not irreducible");
        }

        Ok(Goppa::<T> {
            len: set.len(),
            poly,
            set,
        })
    }

    pub fn random(rng: &mut rand::rngs::ThreadRng, len: usize, t: usize) -> Goppa<T>
    where
        distributions::Standard: distributions::Distribution<T>,
        T: std::fmt::Debug,
    {
        let qm = T::order();
        let poly = Poly::random_irreducible(rng, t);
        let mut set = Vec::new();
        let mut list = Vec::new();
        for i in 0..qm {
            list.push(i);
        }
        for _i in 0..len {
            let j = list.swap_remove(rng.gen_range(0, list.len()));
            let elt = if j == qm - 1 {
                T::zero()
            } else {
                FiniteFieldElement::exp(j)
            };
            set.push(elt);
        }
        Goppa { len, poly, set }
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

    pub fn parity_check_matrix_x(&self) -> Mat<T> {
        let t = self.poly.degree();
        let mut x = Mat::zero(t, t);
        for i in 0..t {
            for j in 0..i + 1 {
                x[(i, j)] = self.poly[t - i + j];
            }
        }
        x
    }

    pub fn parity_check_matrix_y(&self) -> Mat<T> {
        let n = self.len();
        let t = self.poly.degree();

        let mut y = Mat::zero(t, n);
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

        let mut z = Mat::zero(n, n);
        for i in 0..n {
            // println!("{} {}", self.set[i], self.poly.eval(self.set[i]));
            z[(i, i)] = self.poly.eval(self.set[i]).inv().unwrap();
        }
        z
    }

    pub fn parity_check_matrix(&self) -> Mat<T> {
        let x = self.parity_check_matrix_x();
        let y = self.parity_check_matrix_y();
        let z = self.parity_check_matrix_z();
        // let xy = Mat::prod(&x, &y);
        // Mat::prod(&xy, &z)
        x * y * z
    }

    // pub fn syndrome(&self, h: Option<&Mat<T>>, x: &RowVec<T>) -> ColVec<T> {}

    // pub fn syndrome_poly(&self, h: Option<&Mat<T>>, x: &RowVec<T>) -> Poly<T> {}

    // TODO: should be changed because elements of the generator live in Fq and not in Fqm
    pub fn generator_matrix(h: &Mat<T>) -> Mat<T> {
        // Check that h is a parity check matrix
        let m = h.rows();
        let n = h.cols();
        let k = n - m;
        let (_u, hs, p) = h.standard_form().unwrap();
        let mut gs = Mat::zero(k, n);
        for i in 0..k {
            gs[(i, i)] = T::one();
            for j in k..n {
                gs[(i, j)] = -hs[(j - k, i)];
            }
        }
        // let mut g = Mat::zero(k, n);
        // g.as_prod(&gs, &p.inverse().unwrap().transpose());
        // g

        // Mat::prod(&gs, &p.inverse().unwrap().transpose())

        // Mat::prod(&gs, &p.transpose())
        gs * p.transpose()
    }

    // pub fn encode(&self, msg: &Mat<F2>) -> Mat<F2> {
    //     let h = self.parity_check_matrix();
    //     let hb = h.binary_form();
    //     let g = Goppa::generator_matrix(&hb);
    //     msg * g
    // }

    pub fn encode(&self, msg: &RowVec<F2>) -> RowVec<F2> {
        let h = self.parity_check_matrix();
        let hb = h.binary_form();
        let g = Goppa::generator_matrix(&hb);
        msg * g
    }

    // pub fn decode(&self, rcv: &Mat<F2>) -> Mat<F2> {
    pub fn decode(&self, rcv: &RowVec<F2>) -> RowVec<F2> {
        let h = self.parity_check_matrix();
        info!("parity check matrix: {:?}", h);
        // let rcv_fq = Mat::from(rcv);
        let rcv_fq = RowVec::from(rcv);
        info!("received word: {:?}", rcv_fq);
        // let syndrome = Mat::prod(&h, &rcv_fq.transpose());
        let syndrome = h * rcv_fq.transpose();
        info!("syndrome: {:?}", syndrome);

        let zero = Mat::zero(syndrome.rows(), 1);
        if syndrome == zero {
            return rcv.clone();
        }

        let mut deg_s_x = syndrome.rows() - 1;
        for i in 0..syndrome.rows() - 1 {
            if syndrome[(i, 0)] != T::zero() {
                break;
            }
            deg_s_x -= 1;
        }
        let mut s_x = Poly::zero(deg_s_x + 1);
        for i in 0..deg_s_x + 1 {
            s_x[i] = syndrome[(syndrome.rows() - 1 - i, 0)];
            // s_x[i] = syndrome[(i, 0)];
        }
        info!("S(x): {:?}", s_x);
        // debug!("this is a debug {}", "message");
        // error!("this is printed by default");

        let mut t = s_x.inverse_modulo(&self.poly);
        info!("T(x) = s(x)^-1 = {:?}", t);
        // t.add(&Poly::x_n(1));
        t += Poly::x_n(1);
        let s = t.square_root_modulo(&self.poly);
        info!("square root t(x) of T(x) + x: {:?}", s);
        // let (mut a, mut b, _, _, _) = Poly::extended_gcd(&s, &self.poly);
        let (mut a, mut b) = Poly::goppa_extended_gcd(&self.poly, &s);
        info!("a(x) = {:?}", a);
        info!("b(x) = {:?}", b);
        a.square();
        b.square();
        // b.mul(&Poly::x_n(1));
        b *= Poly::x_n(1);
        // let aa = Mat::zero(1, n);
        // for i in 0..a.degree() {
        //     aa[(1, i)] = a[i];
        // }
        // let bb = Mat::zero(1, n);
        // for i in 0..b.degree() {
        //     bb[(1, i)] = b[i];
        // }
        // let sigma = Poly::sum(&a, &b); // not necessarily monic
        let sigma = &a + &b;
        info!("sigma(x) = {:?}", sigma);
        let n = self.set.len();
        // let mut err = Mat::zero(1, n);
        let mut err = RowVec::zero(n);

        for i in 0..n {
            // use iterator and map ?
            // err[(0, i)] = if sigma.eval(self.set[i]) == T::zero() {
            err[i] = if sigma.eval(self.set[i]) == T::zero() {
                F2::one()
            } else {
                F2::zero()
            };
        }

        // let cdw = Mat::sum(rcv, &err);
        let cdw = rcv + err;
        cdw
    }
}
