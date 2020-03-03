use crate::finite_field::CharacteristicTwo;
use crate::finite_field::FiniteFieldElement;
use crate::finite_field_2::F2;
use crate::matrix::Mat;
use crate::polynomial::Poly;
use std::fmt;

pub struct Goppa<T> {
    len: usize,
    poly: Poly<T>,
    set: Vec<T>,
}

// trait BinaryCode<T> {
//     pub fn encode(g: &Mat<T>, msg: &Mat<T>) -> Mat<T> {
//         Mat::prod(msg, &g)
//     }

//     pub fn decode(&self, rcv: &Mat<F2>) -> Mat<F2>
//     where
//         T: CharacteristicTwo,
// }

impl<T> Goppa<T>
where
    T: Copy + fmt::Display + FiniteFieldElement,
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

    pub fn parity_check_matrix_x(&self) -> Mat<T> {
        let t = self.poly.degree();
        let mut x = Mat::new(t, t);
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
        let x = self.parity_check_matrix_x();
        let y = self.parity_check_matrix_y();
        let z = self.parity_check_matrix_z();
        let xy = Mat::prod(&x, &y);
        Mat::prod(&xy, &z)
        // Mat::prod(&y, &z)
    }

    pub fn generator_matrix(h: &Mat<T>) -> Mat<T> {
        // Check that h is a parity check matrix
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
}

pub trait BinaryIrreducible {
    fn encode(&self, msg: &Mat<F2>) -> Mat<F2>;
    fn decode(&self, rcv: &Mat<F2>) -> Mat<F2>;
}

impl<T> BinaryIrreducible for Goppa<T>
where
    T: Copy + fmt::Display + FiniteFieldElement + CharacteristicTwo,
{
    fn encode(&self, msg: &Mat<F2>) -> Mat<F2> {
        let h = self.parity_check_matrix();
        let hb = h.binary_form();
        let g = Goppa::generator_matrix(&hb);
        Mat::prod(msg, &g)
    }

    fn decode(&self, rcv: &Mat<F2>) -> Mat<F2> {
        let h = self.parity_check_matrix();
        let rcv_fq = Mat::from(rcv);
        let syndrome = Mat::prod(&h, &rcv_fq.transpose());
        let zero = Mat::new(syndrome.rows(), 1);
        if syndrome == zero {
            return rcv.clone();
        }

        println!("syndrome: {:?}", syndrome);
        let mut deg_s_x = syndrome.rows() - 1;
        for i in 0..syndrome.rows() - 1 {
            if syndrome[(i, 0)] != T::zero() {
                break;
            }
            deg_s_x -= 1;
        }
        let mut s_x = Poly::new(deg_s_x + 1);
        for i in 0..deg_s_x + 1 {
            s_x[i] = syndrome[(deg_s_x - i, 0)];
            // s_x[i] = syndrome[(i, 0)];
        }
        println!("s_x: {:?}", s_x);
        // debug!("this is a debug {}", "message");
        // error!("this is printed by default");

        let mut t = s_x.inverse_modulo(&self.poly);
        println!("T(x) = s(x)^-1 : {:?}", t);
        t.add(&Poly::x_n(1));
        let s = t.square_root_modulo(&self.poly);
        println!("square root t(x) of T(x) + x: {:?}", s);
        // let (mut a, mut b, _, _, _) = Poly::extended_gcd(&s, &self.poly);
        let (mut a, mut b) = Poly::goppa_extended_gcd(&self.poly, &s);
        println!("a(x): {:?}", a);
        println!("b(x): {:?}", b);
        a.square();
        b.square();
        b.mul(&Poly::x_n(1));
        // let aa = Mat::new(1, n);
        // for i in 0..a.degree() {
        //     aa[(1, i)] = a[i];
        // }
        // let bb = Mat::new(1, n);
        // for i in 0..b.degree() {
        //     bb[(1, i)] = b[i];
        // }
        let sigma = Poly::sum(&a, &b); // not necessarily monic
        println!("sigma(x): {:?}", sigma);
        let n = self.set.len();
        let mut err = Mat::new(1, n);
        for i in 0..n {
            // use iterator and map ?
            err[(0, i)] = if sigma.eval(self.set[i]) == T::zero() {
                F2::one()
            } else {
                F2::zero()
            };
        }

        let cdw = Mat::sum(rcv, &err);
        cdw
    }
}
