use rand::Rng;

use mceliece::finite_field::*;
use mceliece::matrix::*;

#[test]
fn matrix_f1024_permutation_random() {
    FFElt::generate_finite_field(1024);
    
    let mut rng = rand::thread_rng();
    let mat: Mat<FFElt> = Mat::permutation_random(&mut rng, 10);
    println!("{:?}", mat);
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f1024_is_invertible() {
    FFElt::generate_finite_field(1024);
    
    let id: Mat<FFElt> = Mat::identity(11);
    println!("{:?}", id);
    assert!(id.is_invertible());

    let mut rng = rand::thread_rng();
    let mat: Mat<FFElt> = Mat::permutation_random(&mut rng, 20);
    println!("{:?}", mat);
    assert!(mat.is_invertible());
}

#[test]
fn matrix_f1024_inverse() {
    FFElt::generate_finite_field(1024);
    
    let id: Mat<FFElt> = Mat::identity(11);
    assert_eq!(id.inverse().as_ref(), Some(&id));

    let mut rng = rand::thread_rng();
    let mat: Mat<FFElt> = Mat::permutation_random(&mut rng, 11);
    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse permutation matrix")
            .inverse()
            .expect("Cannot inverse permutation matrix"),
        mat
    );
}

#[test]
fn matrix_f1024_invertible_random() {
    FFElt::generate_finite_field(1024);
    
    let mut rng = rand::thread_rng();
    let mat: Mat<FFElt> = Mat::invertible_random(&mut rng, 15);
    assert!(mat.is_invertible());

    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse invertible matrix")
            .inverse()
            .expect("Cannot inverse invertible matrix"),
        mat
    );

    let prod = &mat * &mat.inverse().expect("Cannot inverse invertible matrix");
    let id = Mat::identity(15);
    assert_eq!(prod, id);
}

#[test]
fn matrix_f1024_add() {
    FFElt::generate_finite_field(1024);
    
    let mut rng = rand::thread_rng();
    let a: Mat<FFElt> = Mat::random(&mut rng, 11, 11);
    let b: Mat<FFElt> = Mat::random(&mut rng, 11, 11);
    let c: Mat<FFElt> = Mat::random(&mut rng, 11, 11);
    let z = Mat::zero(11, 11);
    println!("{:?}", a);

    // Associativity
    assert_eq!((&a + &b) + &c, &a + (&b + &c));

    // Commutativity
    assert_eq!(&a + &b, &b + &a);

    // Neutral element
    assert_eq!(&a + &z, a);
    
    // Characteristic
    assert_eq!(&a + &a, z);
}

#[test]
#[should_panic]
fn matrix_f1024_mul_wrong_dimensions() {
    FFElt::generate_finite_field(1024);
    
    let mut prod: Mat<FFElt> = Mat::zero(5, 5);
    let mat1 = Mat::zero(5, 4);
    let mat2 = Mat::zero(3, 5);
    prod.prod(&mat1, &mat2);
}

#[test]
fn matrix_f1024_mul() {
    FFElt::generate_finite_field(1024);
    
    let mut rng = rand::thread_rng();
    let a: Mat<FFElt> = Mat::random(&mut rng, 10, 8);
    let b: Mat<FFElt> = Mat::random(&mut rng, 8, 13);
    let c: Mat<FFElt> = Mat::random(&mut rng, 13, 4);

    // Associativity
    assert_eq!((&a * &b) * &c, &a * (&b * &c));

    // Neutral element
    let i8 = Mat::identity(8);
    let i10 = Mat::identity(10);
    assert_eq!(&a * &i8, a);
    assert_eq!(&i10 * &a, a);

    // Zero case
    let z8 = Mat::zero(8, 8);
    let z10 = Mat::zero(10, 10);
    let z10_8 = Mat::zero(10, 8);
    assert_eq!(&a * &z8, z10_8);
    assert_eq!(&z10 * &a, z10_8);

    // Distributivity
    let a: Mat<F2> = Mat::random(&mut rng, 10, 12);
    let b: Mat<F2> = Mat::random(&mut rng, 10, 12);
    let c: Mat<F2> = Mat::random(&mut rng, 12, 9);
    let d: Mat<F2> = Mat::random(&mut rng, 12, 9);

    // Left: (a + b)c = ac + bc
    assert_eq!((&a + &b) * &c, &a * &c + &b * &c);

    // Right: a(b + c) = ab + ac
    assert_eq!(&a * (&c + &d), &a * &c + &a * &d);
}

#[test]
fn matrix_f1024_rank() {
    FFElt::generate_finite_field(1024);
    
    let mat: Mat<FFElt> = Mat::zero(23, 4);
    assert_eq!(mat.rank(), 0);

    let mat: Mat<FFElt> = Mat::identity(19);
    assert_eq!(mat.rank(), 19);

    let mut rng = rand::thread_rng();
    let mat: Mat<FFElt> = Mat::permutation_random(&mut rng, 34);
    assert_eq!(mat.rank(), 34);
}

#[test]
fn matrix_f1024_weighted_vector_random() {
    FFElt::generate_finite_field(1024);
    
    let mat: Mat<FFElt> = Mat::zero(3, 4);
    assert!(mat.weight() == None);

    let mut rng = rand::thread_rng();
    let vec: Mat<FFElt> = Mat::weighted_vector_random(&mut rng, 35, 13);
    assert_eq!(vec.weight().expect("Cannot compute vector's weight"), 13);
}

#[test]
fn matrix_f1024_standard_form() {
    FFElt::generate_finite_field(1024);
    
    let mat: Mat<FFElt> = Mat::identity(19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert_eq!(h, mat);
    assert!(p.is_permutation());

    let mut h: Mat<FFElt> = Mat::random(&mut rng, 13, 31);
    let inv: Mat<FFElt> = Mat::invertible_random(&mut rng, 13);
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
    assert_eq!(s, u * h * p);
}

#[test]
fn matrix_f1024_transpose() {
    FFElt::generate_finite_field(1024);
    
    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat: Mat<FFElt> = Mat::random(&mut rng, rows, cols);
    assert_eq!(mat, mat.transpose().transpose());
}
