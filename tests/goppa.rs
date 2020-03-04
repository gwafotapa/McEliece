extern crate log;
extern crate mceliece;

use log::{info, warn};

use mceliece::finite_field::*;
use mceliece::finite_field_1024::*;
use mceliece::finite_field_2::*;
use mceliece::finite_field_8::*;
use mceliece::goppa::*;
use mceliece::matrix::*;
use mceliece::polynomial::*;

use std::sync::Once;

static INIT: Once = Once::new();

/// Setup function that is only run once, even if called multiple times.
fn setup() {
    INIT.call_once(|| {
        env_logger::init();
    });
}

#[test]
fn goppa_f8() {
    setup();

    let mut p = Poly::<F8>::new(3);
    p[2] = F8::one();
    p[1] = F8::one();
    p[0] = F8::one();
    println!("{:?}", p);
    let c: Goppa<F8> = Goppa::new(
        p,
        [F8(0), F8(1), F8(2), F8(4), F8(3), F8(6), F8(7), F8(5)].to_vec(),
    )
    .unwrap();
    // let x = c.parity_check_matrix_x();
    // println!("{:?}", x);
    // let y = c.parity_check_matrix_y();
    // println!("{:?}", y);
    // let z = c.parity_check_matrix_z();
    // println!("{:?}", z);
    let h = c.parity_check_matrix();
    // println!("h: {:?}", h);

    // info!("test");
    let h2 = h.binary_form();
    // info!("parity check matrix in binary form: {:?}", h2);
    // info!("rank: {}", h2.rank());
    // return;
    let g = Goppa::generator_matrix(&h2);
    // println!("g: {:?}", g);

    // let mut m = Mat::zero(h.rows(), g.rows());
    // m.as_prod(&h, &g.transpose());
    let m2 = Mat::prod(&h2, &g.transpose());
    // warn!("h2.g^T: {:?}", m2);
    let z2 = Mat::zero(h2.rows(), g.rows());
    assert_eq!(m2, z2);

    let gg = Mat::from(&g);
    let m = Mat::prod(&h, &gg.transpose());
    // warn!("h.g^T: {:?}", m);
    let z = Mat::zero(h.rows(), gg.rows());
    assert_eq!(m, z);
    // return;
    // let h_bin_form = h.binary_form(3);
    // println!("{:?}", h_bin_form);

    let mut rng = rand::thread_rng();
    let msg = Mat::random(&mut rng, 1, 2);
    // // let mut msg = Mat::zero(1, 2);
    // // msg[(0, 0)] = F2::one();
    // println!("msg: {:?}", msg);
    // // let cdw = Mat::prod(&msg, &g);
    // // println!("cdw: {:?}", cdw);
    // let cdw = c.encode(&msg);
    // println!("cdw: {:?}", cdw);
    // let dcd = c.decode(&cdw);
    // println!("dcd: {:?}", dcd);
    // assert_eq!(cdw, dcd);

    let cdw = c.encode(&msg);
    println!("cdw: {:?}", cdw);
    let err = Mat::weighted_vector_random(&mut rng, 8, 2);
    // // let mut err = Mat::zero(1, 8);
    // // err[(0, 1)] = F2::one();
    println!("err: {:?}", err);
    let rcv = Mat::sum(&cdw, &err);
    println!("rcv: {:?}", rcv);
    let dcd = c.decode(&rcv);
    println!("dcd: {:?}", dcd);
    assert_eq!(cdw, dcd);

    // Problem in the maths: the rcv word which is not a codeword has a null syndrome,
    // that is the syndrome is a multiple of the goppa polynomial.
}

#[test]
fn goppa_f1024() {
    setup();

    let mut rng = rand::thread_rng();
    let c: Goppa<F1024> = Goppa::random(&mut rng, 30, 2);
    let h = c.parity_check_matrix();
    info!("parity check matrix: {:?}", h);
    let g = Goppa::generator_matrix(&h);
    assert_eq!(Mat::prod(&h, &g.transpose()), Mat::zero(h.rows(), g.rows()));

    let h2 = h.binary_form();
    info!("parity check matrix in binary form: {:?}", h2);
    info!("rank: {}", h2.rank());
    let g = Goppa::generator_matrix(&h2);
    // assert_eq!(Mat::prod(&h, &Mat::from(&g.transpose())), Mat::zero(h.rows(), g.rows()));
    assert_eq!(Mat::prod(&h2, &g.transpose()), Mat::zero(h2.rows(), g.rows()));


    
    // info!("generator matrix: {:?}", g2);
    // let g = Mat::from(&g2);
    // assert_eq!(Mat::prod(&h, &g.transpose()), Mat::zero(h.rows(), g.rows()));
    let cdw = Mat::new(1, 30, g.data()[0..30].to_vec());
    info!("codeword: {:?}", cdw);
    let syndrome_cdw = Mat::prod(&h, &Mat::from(&cdw.transpose()));
    let zero = Mat::zero(syndrome_cdw.rows(), syndrome_cdw.cols());
    assert_eq!(syndrome_cdw, zero);
    let err: Mat<F2> = Mat::weighted_vector_random(&mut rng, 30, 2);
    info!("error: {:?}", err);
    let rcv = Mat::sum(&cdw, &err);
    info!("received word: {:?}", rcv);
    let dcd = c.decode(&rcv);
    info!("decoded word: {:?}", dcd);
    assert_eq!(cdw, dcd);
}
