use log::info;

use rand::{rngs::ThreadRng, Rng};

use crate::{
    finite_field::{F2FiniteExtension, Field, F2},
    matrix::{Mat, RowVec},
    polynomial::Poly,
};

pub struct Goppa<'a, F: Eq + Field> {
    // len: usize,
    poly: Poly<'a, F>,
    set: Vec<F::FElt>,
}

impl<'a, F: Eq + Field> Goppa<'a, F>
where
    // F: CharacteristicTwo + FiniteField
    F: F2FiniteExtension,
{
    pub fn new(poly: Poly<'a, F>, set: Vec<F::FElt>) -> Result<Goppa<'a, F>, &'static str>
// where
    //     T: std::fmt::Debug,
    {
        if !poly.is_irreducible() {
            return Err("Goppa polynomial is not irreducible");
        }

        Ok(Self {
            // len: set.len(),
            poly,
            set,
        })
    }

    pub fn random(rng: &mut ThreadRng, field: &'a F, len: usize, degree: usize) -> Self
// where
    //     distributions::Standard: distributions::Distribution<T>,
    //     T: std::fmt::Debug,
    {
        let q = field.order();
        let poly = Poly::random_irreducible(rng, field, degree);
        let mut set = Vec::new();
        let mut list = Vec::new();
        for i in 0..q {
            list.push(i);
        }
        for _i in 0..len {
            let j = list.swap_remove(rng.gen_range(0, list.len()));
            let elt = if j == q - 1 {
                field.zero()
            } else {
                field.exp(j)
            };
            set.push(elt);
        }
        Self { poly, set }
    }

    pub fn len(&self) -> usize {
        self.set.len()
    }

    // TODO: Should I return references here or create new objects ?
    pub fn poly(&self) -> &Poly<'a, F> {
        &self.poly
    }

    // TODO: Should I return references here or create new objects ?
    pub fn set(&self) -> &Vec<F::FElt> {
        &self.set
    }

    pub fn field(&self) -> &'a F {
        self.poly.field()
    }

    pub fn parity_check_matrix_x(&self) -> Mat<'a, F> {
        let f = self.field();
        let t = self.poly.degree();
        let mut x = Mat::zero(f, t, t);
        for i in 0..t {
            for j in 0..i + 1 {
                x[(i, j)] = self.poly[t - i + j];
            }
        }
        x
    }

    pub fn parity_check_matrix_y(&self) -> Mat<'a, F> {
        let f = self.field();
        let n = self.len();
        let t = self.poly.degree();

        let mut y = Mat::zero(f, t, n);
        for i in 0..n {
            y[(0, i)] = f.one();
        }
        for i in 1..t {
            for j in 0..n {
                let elt = self.set[j];
                y[(i, j)] = f.mul(elt, y[(i - 1, j)]);
            }
        }
        y
    }

    pub fn parity_check_matrix_z(&self) -> Mat<'a, F> {
        let f = self.field();
        let n = self.len();

        let mut z = Mat::zero(f, n, n);
        for i in 0..n {
            // println!("{} {}", self.set[i], self.poly.eval(self.set[i]));
            z[(i, i)] = f.inv(self.poly.eval(self.set[i])).unwrap();
        }
        z
    }

    pub fn parity_check_matrix(&self) -> Mat<'a, F> {
        let x = self.parity_check_matrix_x();
        let y = self.parity_check_matrix_y();
        let z = self.parity_check_matrix_z();
        // let xy = Mat::prod(&x, &y);
        // Mat::prod(&xy, &z)
        x * y * z
    }

    // pub fn syndrome(&self, h: Option<&Mat<'a, F>>, x: &RowVec<'a, F>) -> ColVec<'a, F> {}

    // pub fn syndrome_poly(&self, h: Option<&Mat<'a, F>>, x: &RowVec<'a, F>) -> Poly<'a, F> {}

    pub fn encode<'b>(&self, msg: &RowVec<'b, F2>) -> RowVec<'b, F2>
    where
        F2: Eq + F2FiniteExtension,
    {
        let f2 = msg.field();
        let h = self.parity_check_matrix();
        let hb = h.binary_form(f2);
        let g = Mat::generator_matrix(&hb);
        msg * g
    }

    pub fn decode<'b>(&self, rcv: &RowVec<'b, F2>) -> RowVec<'b, F2>
    where
        F2: Eq + F2FiniteExtension,
    {
        let f = self.field();
        let f2 = rcv.field();
        let h = self.parity_check_matrix();
        info!("parity check matrix: {:?}", h);
        // let rcv_fq = Mat::from(rcv);
        let rcv_fq = RowVec::from(f, rcv);
        info!("received word: {:?}", rcv_fq);
        // let syndrome = Mat::prod(&h, &rcv_fq.transpose());
        let syndrome = h * rcv_fq.transpose();
        info!("syndrome: {:?}", syndrome);

        let zero = Mat::zero(f, syndrome.rows(), 1);
        if syndrome == zero {
            return rcv.clone();
        }

        let mut deg_s_x = syndrome.rows() - 1;
        for i in 0..syndrome.rows() - 1 {
            if syndrome[(i, 0)] != f.zero() {
                break;
            }
            deg_s_x -= 1;
        }
        let mut s_x = Poly::zero(f, deg_s_x + 1);
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
        t += Poly::x_n(f, 1);
        let s = t.square_root_modulo(&self.poly);
        info!("square root t(x) of T(x) + x: {:?}", s);
        // let (mut a, mut b, _, _, _) = Poly::extended_gcd(&s, &self.poly);
        let (mut a, mut b) = Poly::goppa_extended_gcd(&self.poly, &s);
        info!("a(x) = {:?}", a);
        info!("b(x) = {:?}", b);
        a.square();
        b.square();
        // b.mul(&Poly::x_n(1));
        b *= Poly::x_n(f, 1);
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
        let mut err = RowVec::zero(f2, n);

        for i in 0..n {
            // use iterator and map ?
            // err[(0, i)] = if sigma.eval(self.set[i]) == T::zero() {
            err[i] = if sigma.eval(self.set[i]) == f.zero() {
                f2.one()
            } else {
                f2.zero()
            };
        }

        // let cdw = Mat::sum(rcv, &err);
        let cdw = rcv + err;
        cdw
    }
}
