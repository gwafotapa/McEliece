use rand::Rng;

use mceliece::finite_field::*;
use mceliece::matrix::*;

// #[test]
// fn matrix_f1024_permutation_random() {
//     let f1024 = &F2m::generate(1024);
//     let mut rng = rand::thread_rng();
//     let mat = Mat::permutation_random(&mut rng, f1024, 10);
//     println!("{:?}", mat);
//     assert!(mat.is_permutation());
// }

#[test]
fn matrix_f1024_is_invertible() {
    let f1024 = &F2m::generate(1024);

    let id = Mat::identity(f1024, 11);
    println!("{:?}", id);
    assert!(id.is_invertible());

    // let mut rng = rand::thread_rng();
    // let mat = Mat::permutation_random(&mut rng, f1024, 20);
    // println!("{:?}", mat);
    // assert!(mat.is_invertible());
}

#[test]
fn matrix_f1024_inverse() {
    let f1024 = &F2m::generate(1024);

    let id = Mat::identity(f1024, 11);
    assert_eq!(id.inverse().as_ref(), Some(&id));

    // let mut rng = rand::thread_rng();
    // let mat = Mat::permutation_random(&mut rng, f1024, 11);
    // assert_eq!(
    //     mat.inverse()
    //         .expect("Cannot inverse permutation matrix")
    //         .inverse()
    //         .expect("Cannot inverse permutation matrix"),
    //     mat
    // );
}

#[test]
fn matrix_f1024_invertible_random() {
    let f1024 = &F2m::generate(1024);

    let mut rng = rand::thread_rng();
    let mat = Mat::invertible_random(&mut rng, f1024, 15);
    assert!(mat.is_invertible());

    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse invertible matrix")
            .inverse()
            .expect("Cannot inverse invertible matrix"),
        mat
    );

    let prod = &mat * &mat.inverse().expect("Cannot inverse invertible matrix");
    let id = Mat::identity(f1024, 15);
    assert_eq!(prod, id);
}

#[test]
fn matrix_f1024_add() {
    let f1024 = &F2m::generate(1024);

    let mut rng = rand::thread_rng();
    let a = Mat::random(&mut rng, f1024, 11, 11);
    let b = Mat::random(&mut rng, f1024, 11, 11);
    let c = Mat::random(&mut rng, f1024, 11, 11);
    let z = Mat::zero(f1024, 11, 11);
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
    let f1024 = &F2m::generate(1024);

    let mut prod = Mat::zero(f1024, 5, 5);
    let mat1 = Mat::zero(f1024, 5, 4);
    let mat2 = Mat::zero(f1024, 3, 5);
    prod.prod(&mat1, &mat2);
}

#[test]
fn matrix_f1024_mul() {
    let f1024 = &F2m::generate(1024);

    let mut rng = rand::thread_rng();
    let a = Mat::random(&mut rng, f1024, 10, 8);
    let b = Mat::random(&mut rng, f1024, 8, 13);
    let c = Mat::random(&mut rng, f1024, 13, 4);

    // Associativity
    assert_eq!((&a * &b) * &c, &a * (&b * &c));

    // Neutral element
    let i8 = Mat::identity(f1024, 8);
    let i10 = Mat::identity(f1024, 10);
    assert_eq!(&a * &i8, a);
    assert_eq!(&i10 * &a, a);

    // Zero case
    let z8 = Mat::zero(f1024, 8, 8);
    let z10 = Mat::zero(f1024, 10, 10);
    let z10_8 = Mat::zero(f1024, 10, 8);
    assert_eq!(&a * &z8, z10_8);
    assert_eq!(&z10 * &a, z10_8);

    // Distributivity
    let a = Mat::random(&mut rng, f1024, 10, 12);
    let b = Mat::random(&mut rng, f1024, 10, 12);
    let c = Mat::random(&mut rng, f1024, 12, 9);
    let d = Mat::random(&mut rng, f1024, 12, 9);

    // Left: (a + b)c = ac + bc
    assert_eq!((&a + &b) * &c, &a * &c + &b * &c);

    // Right: a(b + c) = ab + ac
    assert_eq!(&a * (&c + &d), &a * &c + &a * &d);
}

#[test]
fn matrix_f1024_rank() {
    let f1024 = &F2m::generate(1024);

    let mat = Mat::zero(f1024, 23, 4);
    assert_eq!(mat.rank(), 0);

    let mat = Mat::identity(f1024, 19);
    assert_eq!(mat.rank(), 19);

    // let mut rng = rand::thread_rng();
    // let mat = Mat::permutation_random(&mut rng, f1024, 34);
    // assert_eq!(mat.rank(), 34);
}

#[test]
fn matrix_f1024_standard_form() {
    let f1024 = &F2m::generate(1024);

    let mat = Mat::identity(f1024, 19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert_eq!(h, mat);
    assert!(p.is_permutation());

    let mut h = Mat::random(&mut rng, f1024, 13, 31);
    let inv = Mat::invertible_random(&mut rng, f1024, 13);
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
    let f1024 = &F2m::generate(1024);

    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat = Mat::random(&mut rng, f1024, rows, cols);
    assert_eq!(mat, mat.transpose().transpose());
}

#[test]
fn matrix_f1024_rowvec_weight() {
    let f1024 = &F2m::generate(1024);
    let vec = RowVec::zero(f1024, 4);
    assert!(vec.weight() == 0);

    let mut rng = rand::thread_rng();
    let vec = RowVec::random_with_weight(&mut rng, f1024, 35, 13);
    assert_eq!(vec.weight(), 13);
}
