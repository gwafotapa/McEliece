use log::info;
use rand::Rng;

use mceliece::finite_field::*;
use mceliece::matrix::*;

mod common;

#[test]
fn matrix_f7_permutation_random() {
    common::setup();
    let f7 = &F7 {};
    let mut rng = rand::thread_rng();
    let p = Perm::random(&mut rng, 10);
    assert!(p.is_permutation());
}

#[test]
fn matrix_f7_is_invertible() {
    common::setup();
    let f7 = &F7 {};
    let id = Mat::identity(f7, 11);
    info!("{}", id);
    assert!(id.is_invertible());
}

#[test]
fn matrix_f7_inverse() {
    common::setup();
    let f7 = &F7 {};

    let id = Mat::identity(f7, 11);
    assert!(id.inverse().as_ref() == Some(&id));

    let mut rng = rand::thread_rng();
    let p = Perm::random(&mut rng, 11);
    assert_eq!(p, p.inverse().inverse());
}

#[test]
fn matrix_f7_invertible_random() {
    common::setup();
    let f7 = &F7 {};

    let mut rng = rand::thread_rng();
    let mat = Mat::invertible_random(&mut rng, f7, 15);
    assert!(mat.is_invertible());
    assert!(
        mat.inverse()
            .expect("Cannot inverse invertible matrix")
            .inverse()
            .expect("Cannot inverse invertible matrix")
            == mat
    );

    let prod = &mat * &mat.inverse().expect("Cannot inverse invertible matrix");
    let id = Mat::identity(f7, 15);
    assert!(prod == id);
}

#[test]
fn matrix_f7_add() {
    common::setup();
    let f7 = &F7 {};

    let mut rng = rand::thread_rng();
    let a = Mat::random(&mut rng, f7, 11, 11);
    let b = Mat::random(&mut rng, f7, 11, 11);
    let c = Mat::random(&mut rng, f7, 11, 11);
    let z = Mat::zero(f7, 11, 11);
    info!("{}", a);

    // Associativity
    assert!((&a + &b) + &c == &a + (&b + &c));

    // Commutativity
    assert!(&a + &b == &b + &a);

    // Neutral element
    assert!(&a + &z == a);

    // Characteristic
    let mut s = Mat::zero(f7, 11, 11);
    for _i in 0..7 {
        s += &a;
    }
    assert!(s == z);
}

#[test]
#[should_panic]
fn matrix_f7_mul_wrong_dimensions() {
    common::setup();
    let f7 = &F7 {};

    let mut prod = Mat::zero(f7, 5, 5);
    let mat1 = Mat::zero(f7, 5, 4);
    let mat2 = Mat::zero(f7, 3, 5);
    prod.prod(&mat1, &mat2);
}

#[test]
fn matrix_f7_mul() {
    common::setup();
    let f7 = &F7 {};

    let mut rng = rand::thread_rng();
    let a = Mat::random(&mut rng, f7, 10, 8);
    let b = Mat::random(&mut rng, f7, 8, 13);
    let c = Mat::random(&mut rng, f7, 13, 4);

    // Associativity
    assert!((&a * &b) * &c == &a * (&b * &c));

    // Neutral element
    let i8 = Mat::identity(f7, 8);
    let i10 = Mat::identity(f7, 10);
    assert!(&a * &i8 == a);
    assert!(&i10 * &a == a);

    // Zero case
    let z8 = Mat::zero(f7, 8, 8);
    let z10 = Mat::zero(f7, 10, 10);
    let z10_8 = Mat::zero(f7, 10, 8);
    assert!(&a * &z8 == z10_8);
    assert!(&z10 * &a == z10_8);

    // Distributivity
    let a = Mat::random(&mut rng, f7, 10, 12);
    let b = Mat::random(&mut rng, f7, 10, 12);
    let c = Mat::random(&mut rng, f7, 12, 9);
    let d = Mat::random(&mut rng, f7, 12, 9);

    // Left: (a + b)c = ac + bc
    assert!((&a + &b) * &c == &a * &c + &b * &c);

    // Right: a(b + c) = ab + ac
    assert!(&a * (&c + &d) == &a * &c + &a * &d);
}

#[test]
fn matrix_f7_rank() {
    common::setup();
    let f7 = &F7 {};
    let mat = Mat::zero(f7, 23, 4);
    assert!(mat.rank() == 0);

    let mat = Mat::identity(f7, 19);
    assert!(mat.rank() == 19);
}

#[test]
fn matrix_f7_standard_form() {
    common::setup();
    let f7 = &F7 {};

    let mat = Mat::identity(f7, 19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert!(h == mat);
    assert!(p.is_permutation());

    let mut h = Mat::random(&mut rng, f7, 13, 31);
    let inv = Mat::invertible_random(&mut rng, f7, 13);
    for i in 0..13 {
        for j in 0..13 {
            h[(i, j)] = inv[(i, j)];
        }
    }
    let (u, s, p) = h
        .standard_form()
        .expect("Failed to put a full rank matrix in standard form");

    info!("{}", u);
    info!("{}", s);
    info!("{}", p);

    assert!(u.is_invertible());
    assert!(s.is_standard_form());
    assert!(p.is_permutation());
    assert!(s == u * h * p);
}

#[test]
fn matrix_f7_transpose() {
    common::setup();
    let f7 = &F7 {};
    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat = Mat::random(&mut rng, f7, rows, cols);
    let tmat = mat.transpose();
    for i in 0..rows {
        for j in 0..cols {
            assert!(mat[(i, j)] == tmat[(j, i)]);
        }
    }
}

#[test]
fn matrix_f7_rowvec_weight() {
    common::setup();
    let f7 = &F7 {};
    let vec = RowVec::zero(f7, 4);
    assert!(vec.weight() == 0);

    let mut rng = rand::thread_rng();
    let vec = RowVec::random_with_weight(&mut rng, f7, 35, 13);
    assert!(vec.weight() == 13);
}
