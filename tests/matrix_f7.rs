extern crate mceliece;

use rand::Rng;

use mceliece::finite_field_7::*;
use mceliece::matrix::*;

#[test]
fn matrix_f7_permutation_random() {
    let mut rng = rand::thread_rng();
    let mat: Mat<F7> = Mat::permutation_random(&mut rng, 10);
    println!("{:?}", mat);
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f7_is_invertible() {
    let id: Mat<F7> = Mat::identity(11);
    println!("{:?}", id);
    assert!(id.is_invertible());

    let mut rng = rand::thread_rng();
    let mat: Mat<F7> = Mat::permutation_random(&mut rng, 20);
    println!("{:?}", mat);
    assert!(mat.is_invertible());
}

#[test]
fn matrix_f7_inverse() {
    let id: Mat<F7> = Mat::identity(11);
    assert_eq!(id.inverse().as_ref(), Some(&id));

    let mut rng = rand::thread_rng();
    let mat: Mat<F7> = Mat::permutation_random(&mut rng, 11);
    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse permutation matrix")
            .inverse()
            .expect("Cannot inverse permutation matrix"),
        mat
    );
}

#[test]
fn matrix_f7_invertible_random() {
    let mut rng = rand::thread_rng();
    let mat: Mat<F7> = Mat::invertible_random(&mut rng, 15);
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
fn matrix_f7_as_sum() {
    let mut rng = rand::thread_rng();
    let a: Mat<F7> = Mat::random(&mut rng, 11, 11);
    let b: Mat<F7> = Mat::random(&mut rng, 11, 11);
    let c: Mat<F7> = Mat::random(&mut rng, 11, 11);
    let z = Mat::new(11, 11);
    let mut s = Mat::new(11, 11);
    println!("{:?}", a);

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
}

#[test]
#[should_panic]
fn matrix_f7_as_prod_wrong_dimensions() {
    let mut prod: Mat<F7> = Mat::new(5, 5);
    let mat1 = Mat::new(5, 4);
    let mat2 = Mat::new(3, 5);
    prod.as_prod(&mat1, &mat2);
}

#[test]
fn matrix_f7_as_prod() {
    let mut rng = rand::thread_rng();
    let a: Mat<F7> = Mat::random(&mut rng, 10, 8);
    let b: Mat<F7> = Mat::random(&mut rng, 8, 13);
    let c: Mat<F7> = Mat::random(&mut rng, 13, 4);

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
    let mut p: Mat<F7> = Mat::random(&mut rng, 10, 8);
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
    let a: Mat<F7> = Mat::random(&mut rng, 10, 12);
    let b: Mat<F7> = Mat::random(&mut rng, 10, 12);
    let c: Mat<F7> = Mat::random(&mut rng, 12, 9);
    let d: Mat<F7> = Mat::random(&mut rng, 12, 9);
    let mut p: Mat<F7> = Mat::random(&mut rng, 10, 12);
    let mut q: Mat<F7> = Mat::random(&mut rng, 10, 9);
    let mut r: Mat<F7> = Mat::random(&mut rng, 10, 9);
    let mut s: Mat<F7> = Mat::random(&mut rng, 10, 9);
    let mut t: Mat<F7> = Mat::random(&mut rng, 10, 9);
    let mut u: Mat<F7> = Mat::random(&mut rng, 12, 9);

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
fn matrix_f7_rank() {
    let mat: Mat<F7> = Mat::new(23, 4);
    assert_eq!(mat.rank(), 0);

    let mat: Mat<F7> = Mat::identity(19);
    assert_eq!(mat.rank(), 19);

    let mut rng = rand::thread_rng();
    let mat: Mat<F7> = Mat::permutation_random(&mut rng, 34);
    assert_eq!(mat.rank(), 34);
}

#[test]
fn matrix_f7_weighted_vector_random() {
    let mat: Mat<F7> = Mat::new(3, 4);
    assert!(mat.weight() == None);

    let mut rng = rand::thread_rng();
    let vec: Mat<F7> = Mat::weighted_vector_random(&mut rng, 35, 13);
    assert_eq!(vec.weight().expect("Cannot compute vector's weight"), 13);
}

#[test]
fn matrix_f7_standard_form() {
    let mat: Mat<F7> = Mat::identity(19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert_eq!(h, mat);
    assert!(p.is_permutation());

    let mut h: Mat<F7> = Mat::random(&mut rng, 13, 31);
    let inv: Mat<F7> = Mat::invertible_random(&mut rng, 13);
    for i in 0..13 {
        for j in 0..13 {
            h[(i, j)] = inv[(i, j)];
        }
    }
    let (u, s, p) = h
        .standard_form()
        .expect("Failed to put a full rank matrix in standard form");

    println!("{:?}", u);
    println!("{:?}", s);
    println!("{:?}", p);
    
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
fn matrix_f7_transpose() {
    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat: Mat<F7> = Mat::random(&mut rng, rows, cols);
    assert_eq!(mat, mat.transpose().transpose());
}
