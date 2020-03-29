use rand::rngs::ThreadRng;
use std::error::Error;

use super::{Field, Poly};
use crate::finite_field::{f2m, CharacteristicTwo, F2FiniteExtension, FiniteField};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

impl<'a, F: CharacteristicTwo + Eq + Field> Poly<'a, F> {
    pub fn random_monic_irreducible(rng: &mut ThreadRng, f: &'a F, degree: usize) -> Self
    where
        F: FiniteField,
    {
        let mut p = Self::zero(f, degree + 1);
        p[degree] = f.one();
        for i in 0..degree {
            p[i] = f.random_element(rng);
        }
        while !p.is_irreducible() {
            for i in 0..degree {
                p[i] = f.random_element(rng);
            }
        }
        p
    }

    pub fn square(&mut self) {
        let f = self.field;
        let t = self.degree();
        self.data.resize(2 * t + 1, f.zero());

        for i in (1..t + 1).rev() {
            self[2 * i] = f.mul(self[i], self[i]);
            self[2 * i - 1] = f.zero();
        }
        self[0] = f.mul(self[0], self[0]);
    }

    pub fn square_root_modulo(&mut self, modulus: &Self)
    where
        F: FiniteField,
    {
        if self.field != modulus.field {
            panic!("Cannot compute square root modulo: fields don't match");
        }
        let m = self.field.characteristic_exponent();

        for _i in 0..m as usize * modulus.degree() - 1 {
            self.square();
            self.modulo(modulus);
        }
    }

    pub fn inverse_modulo_by_fast_exponentiation(&mut self, modulus: &Self)
    where
        F: FiniteField,
    {
        if self.field != modulus.field {
            panic!("Cannot compute inverse modulo: fields don't match");
        }
        if self.is_zero() {
            panic!("The null polynom has no inverse");
        }
        let m = self.field.characteristic_exponent();

        self.square();
        let mut tmp = self.clone();
        for _i in 0..m as usize * modulus.degree() - 2 {
            tmp.square();
            tmp.modulo(&modulus);
            *self *= &tmp;
            self.modulo(&modulus);
        }
    }

    /// Computes p<sup>n</sup> mod (modulus)
    pub fn pow_modulo(&mut self, mut n: u32, modulus: &Self) {
        if self.field != modulus.field {
            panic!("Cannot compute power modulo: fields don't match");
        }
        self.modulo(modulus);
        let mut tmp = self.clone();
        self.data.truncate(1);
        self[0] = self.field.one();
        if n & 1 == 1 {
            *self *= &tmp;
        }
        n >>= 1;
        while n != 0 {
            tmp.square();
            tmp.modulo(modulus);
            if n & 1 == 1 {
                *self *= &tmp;
                self.modulo(modulus);
            }
            n >>= 1;
        }
    }

    pub fn is_irreducible(&self) -> bool
    where
        F: FiniteField,
    {
        let n = self.degree() as u32;
        if n == 0 {
            return false;
        }

        let f = self.field;
        let q = f.order();
        let mut n_prime_factors = f2m::trial_division(n);
        n_prime_factors.dedup();
        let n_div_primes: Vec<u32> = n_prime_factors.iter().map(|x| n / x).collect();

        for i in 0..n_prime_factors.len() {
            let mut h = Self::x_n(f, 1);
            for _j in 0..n_div_primes[i] {
                h.pow_modulo(q, self);
            }
            h -= &Self::x_n(f, 1);
            let g = Self::gcd(self, &h);
            if g.degree() != 0 {
                return false;
            }
        }
        let mut g = Self::x_n(f, 1);
        for _i in 0..n {
            g.pow_modulo(q, self);
        }
        g -= &Self::x_n(f, 1);
        g.modulo(self);
        g.is_zero()
    }

    pub(crate) fn goppa_extended_gcd(g: &Self, t: &Self) -> (Self, Self) {
        if g.field != t.field {
            panic!("Cannot compute euclidean division: fields differ")
        }

        let f = g.field;
        let mut i = 1;

        let mut a = Vec::new();
        let mut b = Vec::new();
        a.push(g.clone());
        a.push(t.clone());
        let b0 = Self::zero(f, 1);
        let b1 = Self::x_n(f, 0);
        b.push(b0);
        b.push(b1);

        while a[i].degree() > g.degree() / 2 {
            i += 1;
            let (q, r) = Self::euclidean_division(&a[i - 2], &a[i - 1]);
            b.push(&b[i - 2] + &q * &b[i - 1]);
            a.push(r);
        }

        (a.pop().unwrap(), b.pop().unwrap())
    }
}

impl<'a, F: Eq + F2FiniteExtension> Poly<'a, F> {
    pub fn to_hex_string(&self) -> String {
        let f = self.field();
        if f.order() > 1 << 16 {
            panic!("Cannot convert polynomial to hex string: field order must be at most 2^16");
        }
        if self.degree() > 255 {
            panic!("Cannot convert polynomial to hex string: degree must be less than 256");
        }
        let len = 2 * (4 + 1) + 4 * (self.degree() + 1);
        let mut s = String::with_capacity(len);
        s.push_str(format!("{:x}#{:x}#", f.order(), self.degree()).as_str());
        for i in 0..self.degree() {
            s.push_str(format!("{:x} ", f.elt_to_u32(self[i])).as_str());
        }
        s.push_str(format!("{:x}", f.elt_to_u32(self[self.degree()])).as_str());
        s
    }

    pub fn from_hex_string(s: &str, f: &'a F) -> Result<Self> {
        let v: Vec<&str> = s.split('#').collect();
        let t = usize::from_str_radix(v[1], 16)?;
        let v: Vec<&str> = v[2].split(' ').collect();
        let mut poly = Self::zero(&f, t + 1);
        for i in 0..t + 1 {
            poly[i] = f.u32_to_elt(u32::from_str_radix(v[i], 16)?);
        }
        Ok(poly)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_field::{F2m, F2};

    #[test]
    fn is_irreducible() {
        env_logger::init();
        let f2 = &F2 {};
        let f1024 = &F2m::generate(1024);

        let zero = Poly::zero(f2, 1);
        assert!(!zero.is_irreducible());

        let mut constant_poly = Poly::zero(f1024, 1);
        constant_poly[0] = f1024.exp(533);
        assert!(!constant_poly.is_irreducible());

        let p = Poly::support(f2, &[0, 1, 2]);
        assert!(p.is_irreducible());

        let p = Poly::support(f2, &[0, 2]);
        assert!(!p.is_irreducible());

        let p = Poly::support(f1024, &[0, 2, 11]);
        assert!(p.is_irreducible());

        let p = Poly::support(f2, &[0, 3, 10]);
        assert!(p.is_irreducible());

        let p = Poly::support(f1024, &[0, 3, 10]);
        assert!(!p.is_irreducible());
    }
}
