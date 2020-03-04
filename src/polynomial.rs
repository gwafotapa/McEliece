// extern crate rand;

use crate::finite_field::CharacteristicTwo;
use crate::finite_field::FiniteFieldElement;
// use crate::finite_field;

use rand::{distributions, Rng};
use std::ops::{Index, IndexMut};
use std::{cmp, fmt};

#[derive(Clone, Eq, PartialEq)]
pub struct Poly<T>(Vec<T>);

impl<T> Index<usize> for Poly<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> IndexMut<usize> for Poly<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T> fmt::Debug for Poly<T>
where
    T: Copy + fmt::Display + FiniteFieldElement,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.degree() == 0 && self[0] == T::zero() {
            return write!(f, "0");
        }

        let mut s = String::new();
        let mut i = 0;
        while i <= self.degree() {
            if self[i] != T::zero() {
                if self[i] != T::one() || i == 0 {
                    s.push_str(&self[i].to_string());
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

        write!(f, "{}", s)
    }
}

impl<T> Poly<T>
where
    T: Copy + fmt::Display + FiniteFieldElement,
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

    pub fn random(rng: &mut rand::rngs::ThreadRng, degree: usize) -> Poly<T>
    where
        distributions::Standard: distributions::Distribution<T>,
    {
        let mut p = Poly::new(degree + 1);
        for i in 0..degree {
            p[i] = rng.gen();
        }
        while p[degree] == T::zero() {
            p[degree] = rng.gen();
        }
        p
    }

    pub fn random_irreducible(rng: &mut rand::rngs::ThreadRng, degree: usize) -> Poly<T>
    where
        distributions::Standard: distributions::Distribution<T>,
        T: std::fmt::Debug + CharacteristicTwo,
    {
        let mut p = Poly::new(degree + 1);
        for i in 0..degree {
            p[i] = rng.gen();
        }
        while p[degree] == T::zero() {
            p[degree] = rng.gen();
        }
        while !p.is_irreducible() {
            for i in 0..degree {
                p[i] = rng.gen();
            }
            while p[degree] == T::zero() {
                p[degree] = rng.gen();
            }
        }
        p
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

    pub fn sum(p: &Poly<T>, q: &Poly<T>) -> Poly<T> {
        let mut s = p.clone();
        s.add(q);
        s
    }

    // pub fn as_sum(&mut self, p: &Poly<T>, q: &Poly<T>) {
    //     self.0
    //         .resize(1 + cmp::max(p.degree(), q.degree()), T::zero());

    //     for i in 0..self.degree() + 1 {
    //         self[i] = p[i] + q[i];
    //     }

    //     while self.degree() != 0 && self[self.degree()] == T::zero() {
    //         self.0.pop();
    //     }
    // }

    pub fn add(&mut self, p: &Poly<T>) {
        self.0
            .resize(1 + cmp::max(self.degree(), p.degree()), T::zero());

        for i in 0..p.degree() + 1 {
            self[i] += p[i];
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    pub fn diff(p: &Poly<T>, q: &Poly<T>) -> Poly<T> {
        let mut d = p.clone();
        d.sub(q);
        d
    }

    pub fn sub(&mut self, p: &Poly<T>) {
        self.0
            .resize(1 + cmp::max(self.degree(), p.degree()), T::zero());

        for i in 0..p.degree() + 1 {
            self[i] -= p[i];
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    pub fn prod(a: &Poly<T>, b: &Poly<T>) -> Poly<T> {
        let mut c = Poly::new(a.degree() + b.degree() + 1);

        // for i in 0..c.degree() + 1 {
        //     for j in 0..i + 1 {
        //         c[i] += a[j] * b[i - j];
        //     }
        // }

        for i in 0..a.degree() + 1 {
            for j in 0..b.degree() + 1 {
                c[i + j] += a[i] * b[j];
            }
        }

        while c.degree() != 0 && c[c.degree()] == T::zero() {
            c.0.pop();
        }
        c
    }

    pub fn mul(&mut self, a: &Poly<T>) {
        let tmp = self.clone();
        self.0.iter_mut().map(|x| *x = T::zero()).count();
        self.0.resize(self.degree() + a.degree() + 1, T::zero());

        for i in 0..tmp.degree() + 1 {
            for j in 0..a.degree() + 1 {
                self[i + j] += tmp[i] * a[j];
            }
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    // TODO: remove the characteristic 2 bound and modify function accordingly
    pub fn square(&mut self)
    where
        T: CharacteristicTwo,
    {
        // if T::one() + T::one() != T::zero() {
        //     panic!("Square root is only supported for characteristic 2");
        // }
        let t = self.degree();
        self.0.resize(2 * t + 1, T::zero());

        for i in (1..t + 1).rev() {
            self[2 * i] = self[i] * self[i];
            self[2 * i - 1] = T::zero();
        }
        self[0] = self[0] * self[0];
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn euclidean_division(a: &Poly<T>, b: &Poly<T>) -> (Poly<T>, Poly<T>) {
        if b.degree() == 0 && b[0] == T::zero() {
            panic!("Euclidean division by the null polynomial");
        }

        let mut q = Poly::new(1);
        let mut r = a.clone();
        let d = b.degree();
        let c = b[d];
        while r.degree() >= d && (r.degree() != 0 || r[0] != T::zero()) {
            let mut s = Poly::x_n(r.degree() - d);
            s[r.degree() - d] = r[r.degree()] * c.inv().unwrap();
            q.add(&s);
            r.sub(&Poly::prod(&s, &b));
        }
        (q, r)
    }

    pub fn modulo(&mut self, modulus: &Poly<T>) {
        let mut q = Poly::new(1);
        let d = modulus.degree();
        let c = modulus[d];
        while self.degree() >= d && (self.degree() != 0 || self[0] != T::zero()) {
            let mut s = Poly::x_n(self.degree() - d);
            s[self.degree() - d] = self[self.degree()] * c.inv().unwrap();
            q.add(&s);
            self.sub(&Poly::prod(&s, modulus));
        }

        while self.degree() != 0 && self[self.degree()] == T::zero() {
            self.0.pop();
        }
    }

    pub fn gcd(a: &Poly<T>, b: &Poly<T>) -> Poly<T> {
        if b.degree() == 0 && b[b.degree()] == T::zero() {
            return a.clone();
        }

        let (_q, r) = Poly::euclidean_division(a, b);
        Poly::gcd(b, &r)
    }

    pub fn neg(&mut self) {
        for i in 0..self.degree() + 1 {
            self[i] = -self[i];
        }
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
    pub fn extended_gcd(a: &Poly<T>, b: &Poly<T>) -> (Poly<T>, Poly<T>, Poly<T>, Poly<T>, Poly<T>) {
        let mut r: Vec<Poly<T>> = Vec::new();
        let mut s: Vec<Poly<T>> = Vec::new();
        let mut t: Vec<Poly<T>> = Vec::new();
        r.push(a.clone());
        r.push(b.clone());
        let s0 = Poly::x_n(0);
        let s1 = Poly::new(1);
        s.push(s0);
        s.push(s1);
        let t0 = Poly::new(1);
        let t1 = Poly::x_n(0);
        t.push(t0);
        t.push(t1);
        let mut i = 1;
        while r[i].degree() != 0 || r[i][0] != T::zero() {
            let (q, _r) = Poly::euclidean_division(&r[i - 1], &r[i]);
            // println!("{}\n", q.to_str());
            r.push(Poly::diff(&r[i - 1], &Poly::prod(&q, &r[i])));
            s.push(Poly::diff(&s[i - 1], &Poly::prod(&q, &s[i])));
            t.push(Poly::diff(&t[i - 1], &Poly::prod(&q, &t[i])));
            i += 1;
            // println!("{}\n", r[i].to_str());
        }

        let mut a1 = t.remove(i);
        let mut b1 = s.remove(i);
        if i % 2 == 0 {
            a1.neg();
        } else {
            b1.neg();
        }
        let u = s.remove(i - 1);
        let v = t.remove(i - 1);
        let g = r.remove(i - 1);

        (g, u, v, a1, b1)
    }

    pub fn goppa_extended_gcd(g: &Poly<T>, t: &Poly<T>) -> (Poly<T>, Poly<T>) {
        let mut i = 1;

        let mut a: Vec<Poly<T>> = Vec::new();
        let mut b: Vec<Poly<T>> = Vec::new();
        a.push(g.clone());
        a.push(t.clone());
        let b0 = Poly::new(1);
        let b1 = Poly::x_n(0);
        b.push(b0);
        b.push(b1);

        while a[i].degree() > g.degree() / 2 {
            i += 1;
            let (q, r) = Poly::euclidean_division(&a[i - 2], &a[i - 1]);
            // println!("{}\n", q.to_str());
            b.push(Poly::sum(&b[i - 2], &Poly::prod(&q, &b[i - 1])));
            a.push(r);
            // println!("{}\n", r[i].to_str());

            // if (Poly::prod(&b[i], &t).modulo(&g) != a[i].modulo(&g)) {
            //     panic!("Broken loop invariant between a, b, t and g");
            // }
        }

        (a.remove(i), b.remove(i))
    }

    // characteristic 2 only
    // finite field order is 2^m
    pub fn square_root_modulo(&self, modulus: &Poly<T>) -> Poly<T>
    where
        T: CharacteristicTwo,
    {
        let m = T::finite_field_m();
        // if T::one() + T::one() != T::zero() {
        //     panic!("Square root is only supported for characteristic 2");
        // }

        let mut res = self.clone();
        for _i in 0..m as usize * modulus.degree() - 1 {
            res.square();
            // println!("Square {}: {:?}\n", _i, res);
            res.modulo(modulus);
            // println!("Modulo {}: {:?}\n", _i, res);
        }
        res
    }

    // TODO: shouldn't I simply use euclid algo that also works for odd characteristics
    // characteristic 2 only
    // finite field order is 2^m
    pub fn inverse_modulo(&self, modulus: &Poly<T>) -> Poly<T>
    where
        T: CharacteristicTwo,
    {
        let m = T::finite_field_m();
        if self.degree() == 0 && self[0] == T::zero() {
            panic!("The null polynom has no inverse");
        }

        // if T::one() + T::one() != T::zero() {
        //     panic!("Square root is only supported for characteristic 2");
        // }

        let mut res = self.clone();
        res.square();
        let mut tmp = res.clone();
        for _i in 0..m as usize * modulus.degree() - 2 {
            tmp.square();
            // println!("tmp.square() = {:?}\n", tmp);
            tmp.modulo(&modulus);
            // println!("tmp.modulo() = {:?}\n", tmp);
            res.mul(&tmp);
            // println!("res.mul(&tmp) = {:?}\n", res);
            res.modulo(&modulus);
            // println!("res.modulo(&modulus) = {:?}\n", res);
        }
        res
    }

    // Computes p^n mod (modulus)
    pub fn pow_modulo(&mut self, mut n: u32, modulus: &Poly<T>)
    where
        T: CharacteristicTwo,
    {
        let mut tmp = self.clone();
        tmp.modulo(modulus);
        self.0.truncate(1);
        self[0] = T::one();
        if n & 1 == 1 {
            self.mul(&tmp)
        }
        n >>= 1;
        while n != 0 {
            tmp.square();
            tmp.modulo(modulus);
            if n & 1 == 1 {
                self.mul(&tmp);
                self.modulo(modulus);
            }
            n >>= 1;
        }
    }

    pub fn is_irreducible(&self) -> bool
    where
        T: std::fmt::Debug + CharacteristicTwo,
    {
        let n = self.degree() as u32;
        if n == 0 {
            return false;
        }

        // let q = T::finite_field_q();
        // let m = T::finite_field_m();
        let qm = T::finite_field_order();
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
            // let mut xqn_x = Poly::x_n(q.pow(n_div_primes[i]) as usize);
            let mut h = Poly::x_n(1);
            for _j in 0..n_div_primes[i] {
                h.pow_modulo(qm, self);
            }
            h.sub(&Self::x_n(1));
            let g = Poly::gcd(self, &h);
            if g.degree() != 0 {
                // info!("{:?}", g);
                return false;
            }
        }
        let mut g = Poly::x_n(1);
        info!("{:?}", g);
        for _i in 0..n {
            g.pow_modulo(qm, self);
            info!("{:?}", g);
        }
        g.sub(&Poly::x_n(1));
        info!("{:?}", g);
        g.modulo(self);
        info!("{:?}", g);
        // info!("{}", g.degree());
        // info!("{:?}", g[0]);
        g.degree() == 0 && g[0] == T::zero()
    }
}

/// Computes the prime factors of (non zero) integer n by trial division
/// https://en.wikipedia.org/wiki/Trial_division
/// ```
/// # use mceliece::polynomial::trial_division;
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
    use crate::finite_field_1024::F1024;
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
        println!("{:?}", p);

        let mut q: Poly<F2> = Poly::x_n(3);
        q[2] = F2::one();
        println!("{:?}", q);
    }

    #[test]
    fn is_irreducible() {
        env_logger::init();

        let zero: Poly<F2> = Poly::new(1);
        assert!(!zero.is_irreducible());

        let mut constant_poly: Poly<F1024> = Poly::new(1);
        constant_poly[0] = F1024::exp(533);
        assert!(!constant_poly.is_irreducible());

        let mut p: Poly<F2> = Poly::new(3);
        p[0] = F2::one();
        p[1] = F2::one();
        p[2] = F2::one();
        assert!(p.is_irreducible());

        let mut p: Poly<F2> = Poly::new(3);
        p[0] = F2::one();
        p[2] = F2::one();
        assert!(!p.is_irreducible());

        let mut p: Poly<F1024> = Poly::new(12);
        p[0] = F1024::one();
        p[2] = F1024::one();
        p[11] = F1024::one();
        assert!(p.is_irreducible());

        let mut p: Poly<F2> = Poly::new(11);
        p[0] = F2::one();
        p[3] = F2::one();
        p[10] = F2::one();
        assert!(p.is_irreducible());

        let mut p: Poly<F1024> = Poly::new(11);
        p[0] = F1024::one();
        p[3] = F1024::one();
        p[10] = F1024::one();
        assert!(!p.is_irreducible());
    }
}
