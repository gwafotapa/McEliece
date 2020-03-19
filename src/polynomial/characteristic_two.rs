use log::info;

use rand::rngs::ThreadRng;
use std::error::Error;

use super::{Field, Poly};
use crate::finite_field::{CharacteristicTwo, F2FiniteExtension, FiniteField};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

impl<'a, F: CharacteristicTwo + Eq + Field> Poly<'a, F> {
    pub fn random_irreducible(rng: &mut ThreadRng, f: &'a F, degree: usize) -> Self
    where
        F: FiniteField,
    {
        let mut p = Self::zero(f, degree + 1);
        for i in 0..degree {
            p[i] = f.random_element(rng);
        }
        while p[degree] == f.zero() {
            p[degree] = f.random_element(rng);
        }
        while !p.is_irreducible() {
            for i in 0..degree {
                p[i] = f.random_element(rng);
            }
            while p[degree] == f.zero() {
                p[degree] = f.random_element(rng);
            }
        }
        p
    }

    // TODO: remove the characteristic 2 bound and modify function accordingly
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

    // TODO: shouldn't I simply use euclid algo that also works for odd characteristics
    pub fn inverse_modulo(&self, modulus: &Self) -> Self
    where
        F: FiniteField,
    {
        if self.field != modulus.field {
            panic!("Cannot compute inverse modulo: fields don't match");
        }
        let m = self.field.characteristic_exponent();
        if self.is_zero() {
            panic!("The null polynom has no inverse");
        }

        let mut res = self.clone();
        res.square();
        let mut tmp = res.clone();
        for _i in 0..m as usize * modulus.degree() - 2 {
            tmp.square();
            // println!("tmp.square() = {:?}\n", tmp);
            tmp.modulo(&modulus);
            // println!("tmp.modulo() = {:?}\n", tmp);
            res *= &tmp;
            // println!("res.mul(&tmp) = {:?}\n", res);
            res.modulo(&modulus);
            // println!("res.modulo(&modulus) = {:?}\n", res);
        }
        res
    }

    // Computes p^n mod (modulus)
    pub fn pow_modulo(&mut self, mut n: u32, modulus: &Self) {
        if self.field != modulus.field {
            panic!("Cannot compute power modulo: fields don't match");
        }
        let mut tmp = self.clone();
        tmp.modulo(modulus);
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

        // let q = T::characteristic();
        // let m = T::characteristic_exponent();
        let f = self.field;
        let qm = f.order();
        let mut n_prime_factors = trial_division(n);
        n_prime_factors.dedup();
        info!("decomposition of n in prime factors: {:?}", n_prime_factors);

        let k = n_prime_factors.len();
        let mut n_div_primes = Vec::new();
        for i in 0..k {
            n_div_primes.push(n / n_prime_factors[i]);
        }
        info!(
            "list of n/p where p is a prime factor of n: {:?}",
            n_div_primes
        );

        for i in 0..k {
            // let mut xqn_x = Self::x_n(q.pow(n_div_primes[i]) as usize);
            let mut h = Self::x_n(f, 1);
            for _j in 0..n_div_primes[i] {
                h.pow_modulo(qm, self);
            }
            h -= &Self::x_n(f, 1);
            let g = Self::gcd(self, &h);
            if g.degree() != 0 {
                // info!("{:?}", g);
                return false;
            }
        }
        let mut g = Self::x_n(f, 1);
        // info!("{:?}", g);
        for _i in 0..n {
            g.pow_modulo(qm, self);
            // info!("{:?}", g);
        }
        g -= &Self::x_n(f, 1);
        // info!("{:?}", g);
        g.modulo(self);
        // info!("{:?}", g);
        // info!("{}", g.degree());
        // info!("{:?}", g[0]);
        g.is_zero()
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

/// Computes the prime factors of (non zero) integer n by trial division
/// https://en.wikipedia.org/wiki/Trial_division
/// ```
/// # use mceliece::polynomial::characteristic_two::trial_division;
/// assert_eq!(trial_division(1), vec![]);
/// assert_eq!(trial_division(19), vec![19]);
/// assert_eq!(trial_division(77), vec![7, 11]);
/// assert_eq!(trial_division(12), vec![2, 2, 3]);
/// ```
pub fn trial_division(mut n: u32) -> Vec<u32> {
    if n == 0 {
        panic!("0 is an invalid input for trial division");
    }

    let mut prime_factors = Vec::new();
    while n % 2 == 0 {
        prime_factors.push(2);
        n /= 2;
    }
    let mut f = 3;
    while f * f <= n {
        if n % f == 0 {
            prime_factors.push(f);
            n /= f;
        } else {
            f += 2;
        }
    }
    if n != 1 {
        prime_factors.push(n);
    }
    prime_factors
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::finite_field::{F2m, F2};
    // use crate::finite_field::{FiniteFieldElement, F1024, F2};

    // #[test]
    // fn print() {
    //     env_logger::init();

    //     let f2 = &F2 {};
    //     let p = Poly::support(f, &[0, 4, 5]);
    //     println!("{:?}", p);

    //     let q = Poly::support(f, &[2, 3]);
    //     println!("{:?}", q);
    // }

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
