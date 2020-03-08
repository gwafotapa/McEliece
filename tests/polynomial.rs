// extern crate mceliece;

use mceliece::finite_field::*;
use mceliece::polynomial::*;
use rand::Rng;

#[test]
fn polynomial_f2_division() {
    let a: Poly<F2> = Poly::support(&[0, 2, 4]);
    println!("a(x) = {:?}\n", a);
    let b: Poly<F2> = Poly::support(&[0, 1, 2]);
    println!("b(x) = {:?}\n", b);
    let (q, r) = Poly::euclidean_division(&a, &b);
    println!("q(x) = {:?}\n", q);
    println!("r(x) = {:?}\n", r);
    assert_eq!(b, q);
    assert!(r.is_zero());
}

#[test]
fn polynomial_f2_gcd() {
    let a: Poly<F2> = Poly::support(&[0, 1, 4]);
    println!("a(x) = {:?}\n", a);
    let b: Poly<F2> = Poly::support(&[0, 1, 2]);
    println!("b(x) = {:?}\n", b);
    let c: Poly<F2> = Poly::support(&[0, 1]);
    println!("c(x) = {:?}\n", c);
    let ac = &a * &c;
    println!("a(x)c(x) = {:?}\n", ac);
    let bc = &b * &c;
    println!("b(x)c(x) = {:?}\n", bc);
    let d = Poly::gcd(&ac, &bc);
    println!("d(x) = gcd(a(x)c(x), b(x)c(x)) = {:?}\n", d);
    assert_eq!(c, d);
}

#[test]
fn polynomial_f2_extended_gcd() {
    // let mut a: Poly<F2> = Poly::new(5);
    // a[0] = F2::one();
    // a[1] = F2::one();
    // a[4] = F2::one();
    // println!("a(x) = {:?}\n", a);
    // let mut b: Poly<F2> = Poly::new(3);
    // b[0] = F2::one();
    // b[1] = F2::one();
    // b[2] = F2::one();
    // println!("b(x) = {:?}\n", b);
    // let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    // println!("d(x) = {:?}", d);
    // println!("u(x) = {:?}", u);
    // println!("v(x) = {:?}", v);
    // println!("a1(x) = {:?}", a1);
    // println!("b1(x) = {:?}", b1);
    // assert_eq!(a, Poly::prod(&d, &a1));
    // assert_eq!(b, Poly::prod(&d, &b1));
    // assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));

    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(1, 100);
    let a: Poly<F2> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let deg_b = rng.gen_range(0, deg_a);
    let b: Poly<F2> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {:?}\n", b);
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {:?}", d);
    println!("u(x) = {:?}", u);
    println!("v(x) = {:?}", v);
    println!("a1(x) = {:?}", a1);
    println!("b1(x) = {:?}", b1);
    assert_eq!(a, &d * &a1);
    assert_eq!(b, &d * &b1);
    assert_eq!(d, &a * &u + &b * &v);
}

#[test]
fn polynomial_f7_extended_gcd() {
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(1, 100);
    let a: Poly<F7> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let deg_b = rng.gen_range(0, deg_a);
    let b: Poly<F7> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {:?}\n", b);
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {:?}", d);
    println!("u(x) = {:?}", u);
    println!("v(x) = {:?}", v);
    println!("a1(x) = {:?}", a1);
    println!("b1(x) = {:?}", b1);
    assert_eq!(a, &d * &a1);
    assert_eq!(b, &d * &b1);
    assert_eq!(d, &a * &u + &b * &v);
}

#[test]
fn polynomial_f1024_extended_gcd() {
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(1, 100);
    let a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let deg_b = rng.gen_range(0, deg_a);
    let b: Poly<F1024> = Poly::random(&mut rng, deg_b);
    println!("b(x) = {:?}\n", b);
    let (d, u, v, a1, b1) = Poly::extended_gcd(&a, &b);
    println!("d(x) = {:?}", d);
    println!("u(x) = {:?}", u);
    println!("v(x) = {:?}", v);
    println!("a1(x) = {:?}", a1);
    println!("b1(x) = {:?}", b1);
    assert_eq!(a, &d * &a1);
    assert_eq!(b, &d * &b1);
    assert_eq!(d, &a * &u + &b * &v);
}

#[test]
fn polynomial_f1024_square() {
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let b = a.clone();
    a.square();
    assert_eq!(a, &b * &b);
}

#[test]
fn polynomial_f1024_modulo() {
    let mut rng = rand::thread_rng();
    // let deg_a = rng.gen_range(0, 100);
    // let a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    // println!("a(x) = {:?}\n", a);
    // let mut b = a.clone();
    // let g = Poly::x_n(deg_a + 1);
    // b.modulo(&g);
    // assert_eq!(a, b);

    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let deg_g = rng.gen_range(0, 100);
    let g: Poly<F1024> = Poly::random(&mut rng, deg_g);
    println!("g(x) = {:?}\n", g);
    let (_q, r) = Poly::euclidean_division(&a, &g);
    a.modulo(&g);
    assert_eq!(a, r);
}

#[test]
fn polynomial_f1024_sq_root_mod() {
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 11);
    let a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let aa = &a * &a;
    println!("a^2(x) = {:?}\n", aa);
    // let deg_g = rng.gen_range(0, 10);
    // let g = Poly::random(&mut rng, deg_g);
    let g = Poly::support(&[0, 2, 11]);
    // let s = aa.square_root_modulo(&Poly::x_n(100));
    let s = aa.square_root_modulo(&g);
    println!("s(x) = {:?}\n", s);
    assert_eq!(s, a);
}

#[test]
fn polynomial_f1024_inverse_mod() {
    let g = Poly::support(&[0, 2, 11]);

    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 11);
    let a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    // a.modulo(&g);
    println!("a(x) = {:?}\n", a);

    let inv = a.inverse_modulo(&g);
    println!("a^-1(x) = {:?}\n", inv);
    let mut p = &a * &inv;
    p.modulo(&g);
    let id = Poly::x_n(0);
    assert_eq!(p, id);
}
