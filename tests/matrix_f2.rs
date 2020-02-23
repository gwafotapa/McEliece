extern crate mceliece;

use rand::Rng;

use mceliece::finite_field::*;
use mceliece::finite_field_2::*;
use mceliece::matrix::*;

#[test]
fn matrix_f2_new() {
    let mut mat = Mat::new(3, 4);
    mat.set(2, 2, F2::one());
    mat.set(1, 0, F2::one());
    mat.set(2, 3, F2::one());

    let v: Vec<F2> = vec![
        F2::zero(),
        F2::zero(),
        F2::zero(),
        F2::zero(),
        F2::one(),
        F2::zero(),
        F2::zero(),
        F2::zero(),
        F2::zero(),
        F2::zero(),
        F2::one(),
        F2::one(),
    ];

    assert_eq!(3, mat.rows());
    assert_eq!(4, mat.cols());
    assert_eq!(mat.data(), v);
}

#[test]
fn matrix_f2_is_permutation() {
    let mut mat = Mat::new(5, 6);
    mat.set(0, 4, F2::one());
    mat.set(1, 2, F2::one());
    mat.set(2, 3, F2::one());
    mat.set(3, 5, F2::one());
    mat.set(4, 0, F2::one());
    println!();
    mat.print();
    assert!(!mat.is_permutation());

    let mut mat = Mat::new(6, 6);
    mat.set(0, 4, F2::one());
    mat.set(1, 2, F2::one());
    mat.set(3, 3, F2::one());
    mat.set(3, 5, F2::one());
    mat.set(4, 0, F2::one());
    mat.set(5, 1, F2::one());
    println!();
    mat.print();
    assert!(!mat.is_permutation());

    let mut mat = Mat::new(6, 6);
    mat.set(0, 4, F2::one());
    mat.set(1, 3, F2::one());
    mat.set(2, 3, F2::one());
    mat.set(3, 5, F2::one());
    mat.set(4, 0, F2::one());
    mat.set(5, 1, F2::one());
    println!();
    mat.print();
    assert!(!mat.is_permutation());

    let mut mat = Mat::new(6, 6);
    mat.set(0, 4, F2::one());
    mat.set(1, 2, F2::one());
    mat.set(2, 3, F2::one());
    mat.set(3, 5, F2::one());
    mat.set(4, 0, F2::one());
    mat.set(5, 1, F2::one());
    println!();
    mat.print();
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f2_permutation_random() {
    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::permutation_random(&mut rng, 10);
    println!();
    mat.print();
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f2_is_invertible() {
    let id: Mat<F2> = Mat::identity(11);
    println!();
    id.print();
    assert!(id.is_invertible());

    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::permutation_random(&mut rng, 20);
    println!();
    mat.print();
    assert!(mat.is_invertible());
}

#[test]
fn matrix_f2_inverse() {
    let id: Mat<F2> = Mat::identity(11);
    assert_eq!(id.inverse().as_ref(), Some(&id));

    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::permutation_random(&mut rng, 11);
    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse permutation matrix")
            .inverse()
            .expect("Cannot inverse permutation matrix"),
        mat
    );
}

#[test]
fn matrix_f2_invertible_random() {
    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::invertible_random(&mut rng, 15);
    assert!(mat.is_invertible());

    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse invertible matrix")
            .inverse()
            .expect("Cannot inverse invertible matrix"),
        mat
    );

    let mut prod = Mat::new(15, 15);
    prod.as_prod(
        &mat,
        &mat.inverse().expect("Cannot inverse invertible matrix"),
    );
    let id = Mat::identity(15);
    assert_eq!(prod, id);
}

#[test]
fn matrix_f2_as_sum() {
    let mut rng = rand::thread_rng();
    let a: Mat<F2> = Mat::random(&mut rng, 11, 11);
    let b: Mat<F2> = Mat::random(&mut rng, 11, 11);
    let c: Mat<F2> = Mat::random(&mut rng, 11, 11);
    let z = Mat::new(11, 11);
    let mut s = Mat::new(11, 11);
    a.print();

    // Associativity
    let mut ab = Mat::new(11, 11);
    ab.as_sum(&a, &b);
    let mut abc1 = Mat::new(11, 11);
    abc1.as_sum(&ab, &c);
    let mut bc = Mat::new(11, 11);
    bc.as_sum(&b, &c);
    let mut abc2 = Mat::new(11, 11);
    abc2.as_sum(&a, &bc);
    assert_eq!(abc1, abc2);

    // Commutativity
    let mut ba = Mat::new(11, 11);
    ba.as_sum(&b, &a);
    assert_eq!(ab, ba);

    // Neutral element
    s.as_sum(&a, &z);
    assert_eq!(s, a);

    // Characteristic
    s.as_sum(&a, &a);
    assert_eq!(s, z);
}

#[test]
#[should_panic]
fn matrix_f2_as_prod_wrong_dimensions() {
    let mut prod: Mat<F2> = Mat::new(5, 5);
    let mat1 = Mat::new(5, 4);
    let mat2 = Mat::new(3, 5);
    prod.as_prod(&mat1, &mat2);
}

#[test]
fn matrix_f2_as_prod() {
    let mut rng = rand::thread_rng();
    let a: Mat<F2> = Mat::random(&mut rng, 10, 8);
    let b: Mat<F2> = Mat::random(&mut rng, 8, 13);
    let c: Mat<F2> = Mat::random(&mut rng, 13, 4);

    // Associativity
    let mut ab = Mat::new(10, 13);
    ab.as_prod(&a, &b);
    let mut abc1 = Mat::new(10, 4);
    abc1.as_prod(&ab, &c);
    let mut bc = Mat::new(8, 4);
    bc.as_prod(&b, &c);
    let mut abc2 = Mat::new(10, 4);
    abc2.as_prod(&a, &bc);
    assert_eq!(abc1, abc2);

    // Neutral element
    let i8 = Mat::identity(8);
    let i10 = Mat::identity(10);
    let mut p: Mat<F2> = Mat::random(&mut rng, 10, 8);
    p.as_prod(&a, &i8);
    assert_eq!(p, a);
    p.as_prod(&i10, &a);
    assert_eq!(p, a);

    // Zero case
    let z8 = Mat::new(8, 8);
    let z10 = Mat::new(10, 10);
    let z10_8 = Mat::new(10, 8);
    p.as_prod(&a, &z8);
    assert_eq!(p, z10_8);
    p.as_prod(&z10, &a);
    assert_eq!(p, z10_8);

    // Distributivity
    let a: Mat<F2> = Mat::random(&mut rng, 10, 12);
    let b: Mat<F2> = Mat::random(&mut rng, 10, 12);
    let c: Mat<F2> = Mat::random(&mut rng, 12, 9);
    let d: Mat<F2> = Mat::random(&mut rng, 12, 9);
    let mut p: Mat<F2> = Mat::random(&mut rng, 10, 12);
    let mut q: Mat<F2> = Mat::random(&mut rng, 10, 9);
    let mut r: Mat<F2> = Mat::random(&mut rng, 10, 9);
    let mut s: Mat<F2> = Mat::random(&mut rng, 10, 9);
    let mut t: Mat<F2> = Mat::random(&mut rng, 10, 9);
    let mut u: Mat<F2> = Mat::random(&mut rng, 12, 9);

    // Left: (a + b)c = ac + bc
    p.as_sum(&a, &b);
    q.as_prod(&p, &c);
    r.as_prod(&a, &c);
    s.as_prod(&b, &c);
    t.as_sum(&r, &s);
    assert_eq!(q, t);

    // Right: a(c + d) = ac + ad
    u.as_sum(&c, &d);
    t.as_prod(&a, &u);
    q.as_prod(&a, &c);
    r.as_prod(&a, &d);
    s.as_sum(&q, &r);
    assert_eq!(s, t);
}

#[test]
fn matrix_f2_rank() {
    let mat: Mat<F2> = Mat::new(23, 4);
    assert_eq!(mat.rank(), 0);

    let mat: Mat<F2> = Mat::identity(19);
    assert_eq!(mat.rank(), 19);

    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::permutation_random(&mut rng, 34);
    assert_eq!(mat.rank(), 34);
}

#[test]
fn matrix_f2_weighted_vector_random() {
    let mat: Mat<F2> = Mat::new(3, 4);
    assert!(mat.weight() == None);

    let mut rng = rand::thread_rng();
    let vec: Mat<F2> = Mat::weighted_vector_random(&mut rng, 35, 13);
    assert_eq!(vec.weight().expect("Cannot compute vector's weight"), 13);
}

#[test]
fn matrix_f2_standard_form() {
    let mat: Mat<F2> = Mat::identity(19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert_eq!(h, mat);
    assert!(p.is_permutation());

    let mut h: Mat<F2> = Mat::random(&mut rng, 13, 31);
    let inv: Mat<F2> = Mat::invertible_random(&mut rng, 13);
    for i in 0..13 {
        for j in 0..13 {
            h.set(i, j, inv.get(i, j));
        }
    }
    let (u, s, p) = h
        .standard_form()
        .expect("Failed to put a full rank matrix in standard form");

    u.print();
    println!();
    s.print();
    println!();
    p.print();
    println!();
    
    assert!(u.is_invertible());
    assert!(s.is_standard_form());
    assert!(p.is_permutation());
    let mut uh = Mat::new(13, 31);
    uh.as_prod(&u, &h);
    let mut uhp = Mat::new(13, 31);
    uhp.as_prod(&uh, &p);
    assert_eq!(s, uhp);
}

#[test]
fn matrix_f2_transpose() {
    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat: Mat<F2> = Mat::random(&mut rng, rows, cols);
    assert_eq!(mat, mat.transpose().transpose());
}
