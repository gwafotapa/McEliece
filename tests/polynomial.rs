extern crate mceliece;

use rand::Rng;
use mceliece::finite_field::*;
use mceliece::finite_field_2::*;
use mceliece::finite_field_7::*;
use mceliece::finite_field_1024::*;
// use mceliece::goppa::*;
// use mceliece::matrix::*;
use mceliece::polynomial::*;

#[test]
fn polynomial_f2_division() {
    let mut a: Poly<F2> = Poly::new(5);
    a[0] = F2::one();
    a[2] = F2::one();
    a[4] = F2::one();
    println!("a(x) = {}\n", a.to_str());
    let mut b: Poly<F2> = Poly::new(3);
    b[0] = F2::one();
    b[1] = F2::one();
    b[2] = F2::one();
    println!("b(x) = {}\n", b.to_str());
    let (q, r) = Poly::euclidean_division(&a, &b);
    println!("q(x) = {}\n", q.to_str());
    println!("r(x) = {}\n", r.to_str());
    assert_eq!(b, q);
    assert_eq!(r, Poly::new(1));
}

#[test]
fn polynomial_f2_gcd() {
    let mut a: Poly<F2> = Poly::new(5);
    a[0] = F2::one();
    a[1] = F2::one();
    a[4] = F2::one();
    println!("a(x) = {}\n", a.to_str());
    let mut b: Poly<F2> = Poly::new(3);
    b[0] = F2::one();
    b[1] = F2::one();
    b[2] = F2::one();
    println!("b(x) = {}\n", b.to_str());
    let mut c: Poly<F2> = Poly::new(2);
    c[0] = F2::one();
    c[1] = F2::one();
    println!("c(x) = {}\n", c.to_str());
    let mut ac = Poly::prod(&a, &c);
    println!("a(x)c(x) = {}\n", ac.to_str());
    let mut bc = Poly::prod(&b, &c);
    println!("b(x)c(x) = {}\n", bc.to_str());
    let d = Poly::gcd(&ac, &bc);
    println!("d(x) = gcd(a(x)c(x), b(x)c(x)) = {}\n", d.to_str());
    assert_eq!(c, d);
}

#[test]
fn polynomial_f2_extended_gcd() {
    // let mut a: Poly<F2> = Poly::new(5);
    // a[0] = F2::one();
    // a[1] = F2::one();
    // a[4] = F2::one();
    // println!("a(x) = {}\n", a.to_str());
    // let mut b: Poly<F2> = Poly::new(3);
    // b[0] = F2::one();
    // b[1] = F2::one();
    // b[2] = F2::one();
    // println!("b(x) = {}\n", b.to_str());
    // let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    // println!("d(x) = {}", d.to_str());
    // println!("u(x) = {}", u.to_str());
    // println!("v(x) = {}", v.to_str());
    // println!("a1(x) = {}", a1.to_str());
    // println!("b1(x) = {}", b1.to_str());
    // assert_eq!(a, Poly::prod(&d, &a1));
    // assert_eq!(b, Poly::prod(&d, &b1));
    // assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
    
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F2> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {}\n", a.to_str());
    let deg_b = rng.gen_range(0, deg_a);
    let mut b: Poly<F2> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {}\n", b.to_str());
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {}", d.to_str());
    println!("u(x) = {}", u.to_str());
    println!("v(x) = {}", v.to_str());
    println!("a1(x) = {}", a1.to_str());
    println!("b1(x) = {}", b1.to_str());
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}

#[test]
fn polynomial_f7_extended_gcd() {   
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F7> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {}\n", a.to_str());
    let deg_b = rng.gen_range(0, deg_a);
    let mut b: Poly<F7> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {}\n", b.to_str());
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {}", d.to_str());
    println!("u(x) = {}", u.to_str());
    println!("v(x) = {}", v.to_str());
    println!("a1(x) = {}", a1.to_str());
    println!("b1(x) = {}", b1.to_str());
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}

#[test]
fn polynomial_f1024_extended_gcd() {   
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {}\n", a.to_str());
    let deg_b = rng.gen_range(0, deg_a);
    let mut b: Poly<F1024> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {}\n", b.to_str());
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {}", d.to_str());
    println!("u(x) = {}", u.to_str());
    println!("v(x) = {}", v.to_str());
    println!("a1(x) = {}", a1.to_str());
    println!("b1(x) = {}", b1.to_str());
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}
