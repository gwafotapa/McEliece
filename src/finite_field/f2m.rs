use rand::{rngs::ThreadRng, Rng};

use super::{CharacteristicTwo, F2FiniteExtension, Field, FiniteField, F2};

#[derive(Eq)]
pub struct F2m {
    order: u32,
    m: u32, // characteristic exponent
    exp: Vec<<Self as Field>::FElt>,
    log: Vec<u32>,
}

impl PartialEq for F2m {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}

impl Field for F2m {
    type FElt = u32;

    fn generate(order: u32) -> Self {
        let (_, m) = match prime_power(order) {
            Ok(r) => r,
            Err(s) => panic!(s),
        };
        let mut f = Self {
            order,
            m,
            exp: vec![0; order as usize],
            log: vec![0; order as usize],
        };
        f.exp[0] = 1;
        f.log[1] = 0;
        let mut elt = 1;
        for i in 1..order {
            elt *= 2;
            if elt >= order {
                elt ^= Self::primitive_poly(order);
            }
            f.exp[i as usize] = elt;
            f.log[elt as usize] = i;
        }
        f
    }

    fn zero(&self) -> Self::FElt {
        0
    }

    fn one(&self) -> Self::FElt {
        1
    }

    fn characteristic(&self) -> u32 {
        2
    }

    fn add(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        a ^ b
    }

    fn sub(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        self.add(a, b)
    }

    fn mul(&self, a: Self::FElt, b: Self::FElt) -> Self::FElt {
        let modulo = |a| {
            if a >= self.order {
                a - (self.order - 1)
            } else {
                a
            }
        };

        if a == 0 || b == 0 {
            0
        } else {
            self.exp[modulo(self.log[a as usize] + self.log[b as usize]) as usize]
        }
    }

    fn neg(&self, a: Self::FElt) -> Self::FElt {
        a
    }

    fn inv(&self, a: Self::FElt) -> Option<Self::FElt> {
        if a == 0 {
            None
        } else {
            Some(self.exp[(self.order - 1 - self.log[a as usize]) as usize])
        }
    }

    fn random_element(&self, rng: &mut ThreadRng) -> Self::FElt {
        rng.gen_range(0, self.order)
    }
}

impl FiniteField for F2m {
    fn characteristic_exponent(&self) -> u32 {
        self.m
    }

    fn exp(&self, n: u32) -> Self::FElt {
        self.exp[n as usize]
    }

    fn log(&self, a: Self::FElt) -> Option<u32> {
        if a == 0 {
            None
        } else {
            Some(self.log[a as usize])
        }
    }
}

impl CharacteristicTwo for F2m {
    fn from(&self, _f2: &F2, elt: <F2 as Field>::FElt) -> Self::FElt {
        elt
    }
}

impl F2FiniteExtension for F2m {
    fn elt_to_u32(&self, a: Self::FElt) -> u32 {
        a
    }

    fn u32_to_elt(&self, n: u32) -> Self::FElt {
        n
    }
}

impl F2m {
    fn primitive_poly(order: u32) -> u32 {
        match prime_power(order) {
            Ok((_, m)) => match m {
                2 => 0x7,
                3 => 0xB,
                4 => 0x13,
                5 => 0x25,
                6 => 0x43,
                7 => 0x83,
                8 => 0x11D,
                9 => 0x211,
                10 => 0x409,
                11 => 0x805,
                12 => 0x1053,
                13 => 0x201B,
                14 => 0x4143,
                15 => 0x8003,
                16 => 0x110B,
                17 => 0x2009,
                _ => panic!("m must be at least 2 and at most 17"),
            },
            Err(s) => panic!(s),
        }
    }

    pub fn as_set(&self) -> Vec<<F2m as Field>::FElt> {
        let mut set = Vec::with_capacity(self.order as usize);
        for i in 0..self.order {
            set.push(i);
        }
        set
    }
}

/// Determines if number is a prime power.
/// ```
/// # use mceliece::finite_field::f2m::prime_power;
/// assert!(prime_power(2*3).is_err());
/// assert_eq!(prime_power(512), Ok((2, 9)));
/// assert_eq!(prime_power(81), Ok((3, 4)));
/// assert_eq!(prime_power(2_u32.pow(13) - 1), Ok((8191, 1)));
/// ```
pub fn prime_power(q: u32) -> std::result::Result<(u32, u32), &'static str> {
    let mut prime_factors = trial_division(q);
    let m = prime_factors.len() as u32;
    prime_factors.dedup();
    if prime_factors.len() != 1 {
        return Err("Number is not prime");
    }
    let p = prime_factors[0];
    Ok((p, m))
}

/// Computes the prime factors of a nonzero integer by trial division.  
/// https://en.wikipedia.org/wiki/Trial_division
/// ```
/// # use mceliece::finite_field::f2m::trial_division;
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
mod test {
    use super::*;

    #[test]
    fn f256_add() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let c = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.add(a, f.add(b, c)), f.add(f.add(a, b), c));
        assert_eq!(f.add(a, b), f.add(b, a));
        assert_eq!(f.add(a, z), a);
    }

    #[test]
    fn f256_characteristic() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.add(a, a), z);
    }

    #[test]
    fn f256_sub() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);

        assert_eq!(f.add(a, b), f.sub(a, b));
    }

    #[test]
    fn f256_mul() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let c = f.random_element(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.mul(a, f.mul(b, c)), f.mul(f.mul(a, b), c));
        assert_eq!(f.mul(a, b), f.mul(b, a));
        assert_eq!(f.mul(a, i), a);
        assert_eq!(f.mul(a, z), z);
        assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
        assert_eq!(f.mul(a, f.sub(b, c)), f.sub(f.mul(a, b), f.mul(a, c)));
    }

    #[test]
    fn f256_neg() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let b = f.random_element(&mut rng);
        let z = f.zero();

        assert_eq!(f.neg(z), z);
        assert_eq!(f.neg(f.neg(a)), a);
        assert_eq!(f.add(a, f.neg(b)), f.sub(a, b));
    }

    #[test]
    fn f256_inv() {
        let f = F2m::generate(256);
        let mut rng = rand::thread_rng();
        let a = f.random_element(&mut rng);
        let i = f.one();
        let z = f.zero();

        assert_eq!(f.inv(z), None);
        assert_eq!(f.inv(i), Some(i));
        if a != f.zero() {
            assert_eq!(f.inv(f.inv(a).unwrap()), Some(a));
            assert_eq!(f.mul(a, f.inv(a).unwrap()), i);
        }
    }
}
