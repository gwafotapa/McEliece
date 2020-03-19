use log::info;

use rand::{rngs::ThreadRng, Rng};
use std::{error::Error, result};

use crate::{
    finite_field::{F2FiniteExtension, Field, F2},
    matrix::{Mat, RowVec},
    polynomial::Poly,
};

type Result<T> = result::Result<T, Box<dyn Error>>;

#[derive(Eq, PartialEq)]
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
    pub fn new(poly: Poly<'a, F>, set: Vec<F::FElt>) -> result::Result<Goppa<'a, F>, &'static str>
// where
    //     T: std::fmt::Debug,
    {
        if !poly.is_irreducible() {
            return Err("Goppa polynomial is not irreducible");
        }
        let f = poly.field();
        for i in 0..set.len() - 1 {
            if f.elt_to_u32(set[i]) == f.elt_to_u32(set[i + 1]) {
                return Err("Set elements must be different");
            }
            if f.elt_to_u32(set[i]) > f.elt_to_u32(set[i + 1]) {
                return Err("Set elements must be ordered according to their u32 representation");
            }
        }
        Ok(Self {
            // len: set.len(),
            poly,
            set,
        })
    }

    pub fn random(rng: &mut ThreadRng, field: &'a F, n: usize, degree: usize) -> Self
// where
    //     distributions::Standard: distributions::Distribution<T>,
    //     T: std::fmt::Debug,
    {
        let q = field.order() as usize;
        let poly = Poly::random_irreducible(rng, field, degree);
        let mut uints = Vec::with_capacity(q);
        for i in 0..q {
            uints.push(i);
        }
        for _i in 0..q - n {
            let uint = rng.gen_range(0, uints.len());
            uints.remove(uint);
        }
        let mut set = Vec::with_capacity(n);
        for i in 0..n {
            set.push(field.u32_to_elt(uints[i] as u32));
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

    // pub fn generator_matrix<'b, G>(h: &Mat<'b, G>) -> Mat<'b, G>
    // where G: Eq + F2FiniteExtension {
    pub fn generator_matrix(h: &Mat<'a, F>) -> (Mat<'a, F>, Vec<usize>) {
        // pub fn generator_matrix(h: &Mat<'a, F>) -> (Mat<'a, F>, Mat<'a, F>) {
        let f = h.field();
        let m = h.rows();
        let n = h.cols();
        let k = n - m;
        if let Some((_u, hs, p)) = h.standard_form() {
            let mut gs = Mat::zero(f, k, n);
            for i in 0..k {
                gs[(i, i)] = f.one();
                for j in k..n {
                    gs[(i, j)] = f.neg(hs[(j - k, i)]);
                }
            }
            let pt = p.transpose();
            let mut information_set = Vec::with_capacity(k);
            for j in 0..k {
                let mut i = 0;
                while p[(i, j)] == f.zero() {
                    i += 1
                }
                information_set.push(i);
            }
            (gs * pt, information_set)
        // (gs * pt, p)
        } else {
            panic!("Rows of the parity-check matrix aren't independant");
        }
    }

    // TODO: add a syndrome function ?
    // pub fn syndrome(&self, h: Option<&Mat<'a, F>>, x: &RowVec<'a, F>) -> ColVec<'a, F> {}
    // pub fn syndrome_poly(&self, h: Option<&Mat<'a, F>>, x: &RowVec<'a, F>) -> Poly<'a, F> {}

    // pub fn goppa_encode(g: &Mat<'a, F2>, msg: &RowVec<'a, F2>) -> RowVec<'a, F2> {
    //     msg * g
    // }

    pub fn syndrome<'b>(h: &Mat<'a, F>, x: &RowVec<'b, F2>) -> Mat<'a, F> {
        let f2 = x.field();
        let f = h.field();
        let mut s = Mat::zero(f, h.rows(), 1);
        for i in 0..h.rows() {
            for j in 0..x.cols() {
                if x[j] == f2.one() {
                    s[(i, 0)] = f.add(s[(i, 0)], h[(i, j)]);
                }
            }
        }
        s
    }

    pub fn encode<'b>(&self, msg: &RowVec<'b, F2>) -> RowVec<'b, F2> {
        let f2 = msg.field();
        let h = self.parity_check_matrix();
        let hb = h.binary_form(f2);
        let (g, _) = Goppa::generator_matrix(&hb);
        msg * g
        // Self::goppa_encode(&g, msg)
    }

    pub fn decode<'b>(&self, rcv: &RowVec<'b, F2>) -> RowVec<'b, F2> {
        let f = self.field();
        let f2 = rcv.field();
        let h = self.parity_check_matrix();
        info!("parity check matrix: {:?}", h);
        // let rcv_fq = Mat::from(rcv);
        // let rcv_fq = RowVec::from(f, rcv);
        // info!("received word: {:?}", rcv_fq);
        // let syndrome = h * rcv_fq.transpose();
        let syndrome = Goppa::syndrome(&h, &rcv);
        info!("syndrome: {:?}", syndrome);

        // let zero = Mat::zero(f, syndrome.rows(), 1);
        // if syndrome == zero {
        // if syndrome.is_zero() {
        //     return rcv.clone();
        // }

        let mut deg_s_x = syndrome.rows() - 1;
        for i in 0..syndrome.rows() - 1 {
            if syndrome[(i, 0)] != f.zero() {
                break;
            }
            deg_s_x -= 1;
        }
        // If syndrome is zero, we're done.
        if deg_s_x == 0 && syndrome[(syndrome.rows() - 1, 0)] == f.zero() {
            return rcv.clone();
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

// Should I generalise from F2m ?
impl<'a, F: Eq + F2FiniteExtension> Goppa<'a, F>
where
    F::FElt: std::fmt::Debug,
{
    // TODO: remove the trait ?
    pub fn to_hex_string(&self) -> String
// where
    //     F::FElt: std::fmt::LowerHex,
    {
        let f = self.field();
        // TODO: add 'use fmt::LowerHex'
        if f.order() > 1 << 16 {
            panic!("Cannot convert polynomial to hex string: field order must be at most 2^16");
        }
        // let mut s = String::new();
        // s.push_str(format!("{:04x}", f.order()).as_str());
        // s.push_str(format!("{:02x}", self.poly.degree()).as_str());
        println!("{:?}\n", self.poly);
        let mut s = self.poly.to_hex_string(); // TODO: Does it clone the string or not ?
        s.reserve(s.len() + 2 + 2 * (f.order() as usize / 8 + 1) + 1);
        s.push_str("##");
        let mut byte = 0;
        let mut cnt_mod_8 = 7;
        let mut j = 0;
        for i in 0..f.order() {
            if self.set[j] == f.u32_to_elt(i) {
                // if self.set.contains(&f.u32_to_elt(i)) {
                // TODO: rewrite using elt_to_u32
                byte |= 1 << cnt_mod_8;
                j += 1;
            }
            if cnt_mod_8 == 0 {
                s.push_str(format!("{:02x}", byte).as_str());
                byte = 0;
                cnt_mod_8 = 7;
            } else {
                cnt_mod_8 -= 1;
            }
        } // TODO: if field order is not divisible by 8, elements are missing !
        s
    }

    pub fn from_hex_string(s: &str, f: &'a F) -> Result<Self> {
        // let order = u32::from_str_radix(&s[0..4], 16).unwrap();
        // let t = u32::from_str_radix(&s[4..6], 16).unwrap();
        // let poly_len = 4 + 2 + 4 * (t as usize + 1);
        let v: Vec<&str> = s.split("##").collect();
        let poly = Poly::from_hex_string(v[0], f)?;
        println!("{:?}\n", poly);
        // let v: Vec<&str> = v[1].split(' ').collect();
        // let mut set = Vec::<<F as Field>:: FElt>::with_capacity(v.len());
        // for i in 0..v.len() {
        //     set.push(f.u32_to_elt(u32::from_str_radix(v[i], 16).unwrap()));
        // }

        let set_data = hex::decode(v[1])?;
        // println!("{:?}", set_data);
        let mut set = Vec::new();
        let mut iter = set_data.iter();
        let mut byte = iter.next().ok_or("Missing byte")?;
        let mut cnt_mod_8 = 7;
        for i in 0..f.order() {
            if (*byte >> cnt_mod_8) & 1 == 1 {
                set.push(f.u32_to_elt(i));
            }
            if cnt_mod_8 == 0 && i != f.order() - 1 {
                cnt_mod_8 = 7;
                byte = iter.next().ok_or("Missing byte")?;
            } else {
                cnt_mod_8 -= 1
            }
        }
        println!("{:?}\n", set);
        Ok(Goppa { poly, set })
    }
}
