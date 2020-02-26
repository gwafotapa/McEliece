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
    println!("a(x) = {:?}\n", a);
    let mut b: Poly<F2> = Poly::new(3);
    b[0] = F2::one();
    b[1] = F2::one();
    b[2] = F2::one();
    println!("b(x) = {:?}\n", b);
    let (q, r) = Poly::euclidean_division(&a, &b);
    println!("q(x) = {:?}\n", q);
    println!("r(x) = {:?}\n", r);
    assert_eq!(b, q);
    assert_eq!(r, Poly::new(1));
}

#[test]
fn polynomial_f2_gcd() {
    let mut a: Poly<F2> = Poly::new(5);
    a[0] = F2::one();
    a[1] = F2::one();
    a[4] = F2::one();
    println!("a(x) = {:?}\n", a);
    let mut b: Poly<F2> = Poly::new(3);
    b[0] = F2::one();
    b[1] = F2::one();
    b[2] = F2::one();
    println!("b(x) = {:?}\n", b);
    let mut c: Poly<F2> = Poly::new(2);
    c[0] = F2::one();
    c[1] = F2::one();
    println!("c(x) = {:?}\n", c);
    let ac = Poly::prod(&a, &c);
    println!("a(x)c(x) = {:?}\n", ac);
    let bc = Poly::prod(&b, &c);
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
    let deg_a = rng.gen_range(0, 100);
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
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}

#[test]
fn polynomial_f7_extended_gcd() {   
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
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
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}

#[test]
fn polynomial_f1024_extended_gcd() {   
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
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
    assert_eq!(a, Poly::prod(&d, &a1));
    assert_eq!(b, Poly::prod(&d, &b1));
    assert_eq!(d, Poly::sum(&Poly::prod(&a, &u), &Poly::prod(&b, &v)));
}

#[test]
fn polynomial_f1024_square() {
    let mut rng = rand::thread_rng();
    let deg_a = rng.gen_range(0, 100);
    let mut a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let b = a.clone();
    a.square();
    assert_eq!(a, Poly::prod(&b, &b));
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
    let deg_a = rng.gen_range(0, 10);
    let a: Poly<F1024> = Poly::random(&mut rng, deg_a);
    println!("a(x) = {:?}\n", a);
    let aa = Poly::prod(&a, &a);
    println!("a^2(x) = {:?}\n", aa);
    // let deg_g = rng.gen_range(0, 10);
    // let g = Poly::random(&mut rng, deg_g);
    let mut g = Poly::new(11);
    g[10] = F1024::one();
    g[3] = F1024::one();
    g[0] = F1024::one();
    // let s = aa.square_root_modulo(&Poly::x_n(100));
    let s = aa.square_root_modulo(&g);
    println!("s(x) = {:?}\n", s);
    assert_eq!(s, a);
}
