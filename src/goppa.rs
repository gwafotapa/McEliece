use log::info;

use rand::{rngs::ThreadRng, Rng};
use std::{
    error::Error,
    fmt::{self, Debug, Display, Formatter},
    result,
};

use crate::{
    finite_field::{F2FiniteExtension, Field, F2},
    matrix::{Mat, RowVec},
    polynomial::Poly,
};

type Result<T> = result::Result<T, Box<dyn Error>>;

#[derive(Eq, PartialEq)]
pub struct Goppa<'a, F: Eq + Field> {
    poly: Poly<'a, F>,
    set: Vec<F::FElt>,
}

impl<'a, F> Debug for Goppa<'a, F>
where
    F: Eq + F2FiniteExtension,
    F::FElt: Debug,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let q = self.poly.field().order();
        let n = self.set.len();
        let t = self.poly.degree();
        write!(
            f,
            "Goppa code (L, g(x)) from F{} with n={} and t={}\n",
            q, n, t
        )?;
        write!(f, "g(x) = {:?}\n", self.poly)?;
        if n == q as usize {
            write!(f, "L = F{}\n", q)
        } else {
            write!(f, "L = {:X?}\n", self.set)
        }
    }
}

impl<'a, F> Display for Goppa<'a, F>
where
    F: Eq + F2FiniteExtension,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let fq = self.poly.field();
        let q = fq.order();
        let n = self.set.len();
        let t = self.poly.degree();
        write!(
            f,
            "Goppa code (L, g(x)) from F{} with n={} and t={}\n",
            q, n, t
        )?;
        write!(f, "g(x) = {}\n", self.poly)?;
        if n == q as usize {
            write!(f, "L = F{}\n", q)
        } else {
            // write!(f, "L = [")?;
            // for i in 0..n - 1 {
            //     write!(f, "{}, ", fq.elt_to_str(self.set[i]))?;
            // }
            // write!(f, "{}]\n", fq.elt_to_str(self.set[n - 1]))

            let set: Vec<String> = self
                .set
                .iter()
                .map(|x| fq.elt_to_str(*x))
                .collect();
            write!(f, "L = {:?}\n", set)
        }
    }
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
        if poly.degree() == 1 && set.contains(&f.mul(f.inv(poly[1]).unwrap(), poly[0])) {
            return Err("Set must not contain a root of the Goppa polynomial");
        }
        Ok(Self { poly, set })
    }

    pub fn random(rng: &mut ThreadRng, f: &'a F, n: usize, degree: usize) -> Self {
        let poly = Poly::random_monic_irreducible(rng, f, degree);

        // Build set L of n elements :
        // Start with all q elements of F(2^m) then remove n-q elements.
        let q = f.order() as usize;
        let mut uints = Vec::with_capacity(q);
        for i in 0..q {
            uints.push(i);
        }
        let mut elts_to_remove = q - n;

        if poly.degree() == 1 {
            let root = f.mul(f.inv(poly[1]).unwrap(), poly[0]);
            uints.remove(f.elt_to_u32(root) as usize);
            elts_to_remove -= 1;
        }

        for _i in 0..elts_to_remove {
            let index = rng.gen_range(0, uints.len());
            uints.remove(index);
        }
        let set = uints.iter().map(|x| f.u32_to_elt(*x as u32)).collect();

        Self { poly, set }
    }

    pub fn len(&self) -> usize {
        self.set.len()
    }

    pub fn poly(&self) -> &Poly<'a, F> {
        &self.poly
    }

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
            z[(i, i)] = f.inv(self.poly.eval(self.set[i])).unwrap();
        }
        z
    }

    pub fn parity_check_matrix(&self) -> Mat<'a, F> {
        let x = self.parity_check_matrix_x();
        let y = self.parity_check_matrix_y();
        let z = self.parity_check_matrix_z();
        x * y * z
    }

    pub fn generator_matrix(h: &Mat<'a, F>) -> (Mat<'a, F>, Vec<usize>) {
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
        } else {
            panic!("Rows of the parity-check matrix aren't independant");
        }
    }

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
        let mut hb = h.binary_form(f2);
        hb.remove_redundant_rows();
        let (g, _) = Goppa::generator_matrix(&hb);
        msg * g
    }

    pub fn decode<'b>(&self, rcv: &RowVec<'b, F2>) -> RowVec<'b, F2> {
        info!("Goppa decoding\n");
        let f = self.field();
        let f2 = rcv.field();
        let h = self.parity_check_matrix();
        // info!("parity check matrix:{:?}", h);
        let syndrome = Goppa::syndrome(&h, rcv);
        info!("syndrome:{:?}", syndrome);

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
        }
        info!("S(x) = {:?}", s_x);

        s_x.inverse_modulo(&self.poly);
        info!("T(x) = s(x)^-1 = {:?}", s_x);
        s_x += Poly::x_n(f, 1);
        s_x.square_root_modulo(&self.poly);
        info!("(T(x) + x)^(1/2) = {:?}", s_x);
        let (mut a, mut b) = Poly::goppa_extended_gcd(&self.poly, &s_x);
        info!("a(x) = {:?}", a);
        info!("b(x) = {:?}", b);
        a.square();
        b.square();
        info!("a(x)^2 = {:?}", a);
        info!("b(x)^2 = {:?}", b);
        b *= Poly::x_n(f, 1);
        let sigma = &a + &b;
        info!("sigma(x) = {:?}", sigma);

        let err = RowVec::new(
            f2,
            self.set
                .iter()
                .map(|x| {
                    if sigma.eval(*x) == f.zero() {
                        f2.one()
                    } else {
                        f2.zero()
                    }
                })
                .collect(),
        );

        // let n = self.set.len();
        // let mut err = RowVec::zero(f2, n);
        // for i in 0..n {
        //     err[i] = if sigma.eval(self.set[i]) == f.zero() {
        //         f2.one()
        //     } else {
        //         f2.zero()
        //     };
        // }
        info!("Error vector:{}", err);
        let cdw = rcv + err;
        cdw
    }
}

// Should I generalise from F2m ?
impl<'a, F: Eq + F2FiniteExtension> Goppa<'a, F>
// where
//     F::FElt: std::fmt::Debug,
{
    pub fn to_hex_string(&self) -> String {
        let f = self.field();
        if f.order() > 1 << 16 {
            panic!("Cannot convert polynomial to hex string: field order must be at most 2^16");
        }
        // println!("{:?}\n", self.poly);
        let mut s = self.poly.to_hex_string();
        s.reserve(s.len() + 2 + 2 * (f.order() as usize / 8 + 1) + 1);
        s.push_str("##");
        let mut byte = 0;
        let mut cnt_mod_8 = 7;
        let mut j = 0;
        for i in 0..f.order() {
            if f.u32_to_elt(i) == self.set[j] {
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
        }
        if f.order() % 8 != 0 {
            s.push_str(format!("{:02x}", byte).as_str());
        }
        s
    }

    pub fn from_hex_string(s: &str, f: &'a F) -> Result<Self> {
        let v: Vec<&str> = s.split("##").collect();
        let poly = Poly::from_hex_string(v[0], f)?;
        // println!("{:?}\n", poly);

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
        // println!("{:?}\n", set);
        Ok(Goppa { poly, set })
    }
}
