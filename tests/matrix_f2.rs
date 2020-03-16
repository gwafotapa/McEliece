use rand::Rng;

use mceliece::finite_field::*;
use mceliece::matrix::*;

#[test]
fn matrix_f2_new() {
    let f2 = &F2 {};
    let mut mat = Mat::zero(f2, 3, 4);
    mat[(2, 2)] = 1;
    mat[(1, 0)] = 1;
    mat[(2, 3)] = 1;
    let v: Vec<<F2 as Field>::FElt> = vec![0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1];

    assert_eq!(mat.rows(), 3);
    assert_eq!(mat.cols(), 4);
    assert_eq!(mat.data(), v);
}

#[test]
fn matrix_f2_is_permutation() {
    let f2 = &F2 {};
    let mut mat = Mat::zero(f2, 5, 6);
    mat[(0, 4)] = 1;
    mat[(1, 2)] = 1;
    mat[(2, 3)] = 1;
    mat[(3, 5)] = 1;
    mat[(4, 0)] = 1;
    println!("{:?}", mat);
    assert!(!mat.is_permutation());

    let mut mat = Mat::zero(f2, 6, 6);
    mat[(0, 4)] = 1;
    mat[(1, 2)] = 1;
    mat[(3, 3)] = 1;
    mat[(3, 5)] = 1;
    mat[(4, 0)] = 1;
    mat[(5, 1)] = 1;
    println!("{:?}", mat);
    assert!(!mat.is_permutation());

    let mut mat = Mat::zero(f2, 6, 6);
    mat[(0, 4)] = 1;
    mat[(1, 3)] = 1;
    mat[(2, 3)] = 1;
    mat[(3, 5)] = 1;
    mat[(4, 0)] = 1;
    mat[(5, 1)] = 1;
    println!("{:?}", mat);
    assert!(!mat.is_permutation());

    let mut mat = Mat::zero(f2, 6, 6);
    mat[(0, 4)] = 1;
    mat[(1, 2)] = 1;
    mat[(2, 3)] = 1;
    mat[(3, 5)] = 1;
    mat[(4, 0)] = 1;
    mat[(5, 1)] = 1;
    println!("{:?}", mat);
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f2_permutation_random() {
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let mat = Mat::permutation_random(&mut rng, f2, 10);
    println!("{:?}", mat);
    assert!(mat.is_permutation());
}

#[test]
fn matrix_f2_is_invertible() {
    let f2 = &F2 {};
    let id = Mat::identity(f2, 11);
    println!("{:?}", id);
    assert!(id.is_invertible());

    let mut rng = rand::thread_rng();
    let mat = Mat::permutation_random(&mut rng, f2, 20);
    println!("{:?}", mat);
    assert!(mat.is_invertible());
}

#[test]
fn matrix_f2_inverse() {
    let f2 = &F2 {};
    let id = Mat::identity(f2, 11);
    assert_eq!(id.inverse().as_ref(), Some(&id));

    let mut rng = rand::thread_rng();
    let mat = Mat::permutation_random(&mut rng, f2, 11);
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
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let mat = Mat::invertible_random(&mut rng, f2, 15);
    assert!(mat.is_invertible());

    assert_eq!(
        mat.inverse()
            .expect("Cannot inverse invertible matrix")
            .inverse()
            .expect("Cannot inverse invertible matrix"),
        mat
    );

    let prod = &mat * &mat.inverse().expect("Cannot inverse invertible matrix");
    let id = Mat::identity(f2, 15);
    assert_eq!(prod, id);
}

#[test]
fn matrix_f2_add() {
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let a = Mat::random(&mut rng, f2, 11, 11);
    let b = Mat::random(&mut rng, f2, 11, 11);
    let c = Mat::random(&mut rng, f2, 11, 11);
    let z = Mat::zero(f2, 11, 11);
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
fn matrix_f2_mul_wrong_dimensions() {
    let f2 = &F2 {};
    let mut prod: Mat<F2> = Mat::zero(f2, 5, 5);
    let mat1 = Mat::zero(f2, 5, 4);
    let mat2 = Mat::zero(f2, 3, 5);
    prod.prod(&mat1, &mat2);
}

#[test]
fn matrix_f2_mul() {
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let a: Mat<F2> = Mat::random(&mut rng, f2, 10, 8);
    let b: Mat<F2> = Mat::random(&mut rng, f2, 8, 13);
    let c: Mat<F2> = Mat::random(&mut rng, f2, 13, 4);

    // Associativity
    assert_eq!((&a * &b) * &c, &a * (&b * &c));

    // Neutral element
    let i8 = Mat::identity(f2, 8);
    let i10 = Mat::identity(f2, 10);
    assert_eq!(&a * &i8, a);
    assert_eq!(&i10 * &a, a);

    // Zero case
    let z8 = Mat::zero(f2, 8, 8);
    let z10 = Mat::zero(f2, 10, 10);
    let z10_8 = Mat::zero(f2, 10, 8);
    assert_eq!(&a * &z8, z10_8);
    assert_eq!(&z10 * &a, z10_8);

    // Distributivity
    let a: Mat<F2> = Mat::random(&mut rng, f2, 10, 12);
    let b: Mat<F2> = Mat::random(&mut rng, f2, 10, 12);
    let c: Mat<F2> = Mat::random(&mut rng, f2, 12, 9);
    let d: Mat<F2> = Mat::random(&mut rng, f2, 12, 9);

    // Left: (a + b)c = ac + bc
    assert_eq!((&a + &b) * &c, &a * &c + &b * &c);

    // Right: a(b + c) = ab + ac
    assert_eq!(&a * (&c + &d), &a * &c + &a * &d);
}

#[test]
fn matrix_f2_rank() {
    let f2 = &F2 {};
    let mat: Mat<F2> = Mat::zero(f2, 23, 4);
    assert_eq!(mat.rank(), 0);

    let mat: Mat<F2> = Mat::identity(f2, 19);
    assert_eq!(mat.rank(), 19);

    let mut rng = rand::thread_rng();
    let mat: Mat<F2> = Mat::permutation_random(&mut rng, f2, 34);
    assert_eq!(mat.rank(), 34);
}

#[test]
fn matrix_f2_standard_form() {
    let f2 = &F2 {};
    let mat = Mat::identity(f2, 19);
    assert!(mat.is_standard_form());

    let mut rng = rand::thread_rng();
    let (u, h, p) = mat
        .standard_form()
        .expect("Cannot recognize the identity matrix as a standard form");
    assert!(u.is_permutation());
    assert_eq!(h, mat);
    assert!(p.is_permutation());

    let mut h = Mat::random(&mut rng, f2, 13, 31);
    let inv = Mat::invertible_random(&mut rng, f2, 13);
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
fn matrix_f2_transpose() {
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let rows = rng.gen_range(0, 100);
    let cols = rng.gen_range(0, 100);
    let mat = Mat::random(&mut rng, f2, rows, cols);
    assert_eq!(mat, mat.transpose().transpose());
}

#[test]
fn matrix_f2_rowvec_weight() {
    let f2 = &F2 {};
    let vec = RowVec::zero(f2, 4);
    assert!(vec.weight() == 0);

    let mut rng = rand::thread_rng();
    let vec = RowVec::random_with_weight(&mut rng, f2, 35, 13);
    assert_eq!(vec.weight(), 13);
}

#[test]
fn rowvec_f2_save_load() {
    let f2 = &F2 {};
    let mut rng = rand::thread_rng();
    let n = rng.gen_range(10, 1000);
    let vec_save = RowVec::random(&mut rng, f2, n);
    let file_name = "vec_save_load_test";
    vec_save.save_vector(file_name);
    let vec_load = RowVec::load_vector(file_name, f2);
    assert_eq!(vec_save, vec_load);
}
