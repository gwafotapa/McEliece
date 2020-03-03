extern crate log;
extern crate mceliece;

use log::info;

use mceliece::finite_field::*;
use mceliece::finite_field_2::*;
use mceliece::finite_field_8::*;
use mceliece::goppa::*;
use mceliece::matrix::*;
use mceliece::polynomial::*;

#[test]
fn goppa_f8() {
    env_logger::init();

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
    println!("h: {:?}", h);

    // info!("test");
    let hb = h.binary_form();
    let g = Goppa::generator_matrix(&hb);
    println!("g: {:?}", g);

    let z = Mat::new(hb.rows(), g.rows());
    // let mut m = Mat::new(h.rows(), g.rows());
    // m.as_prod(&h, &g.transpose());
    let m = Mat::prod(&hb, &g.transpose());
    assert_eq!(m, z);

    // let h_bin_form = h.binary_form(3);
    // println!("{:?}", h_bin_form);

    let mut rng = rand::thread_rng();
    // let msg = Mat::random(&mut rng, 1, 2);
    let mut msg = Mat::new(1, 2);
    msg[(0, 0)] = F2::one();
    println!("msg: {:?}", msg);
    // let cdw = Mat::prod(&msg, &g);
    // println!("cdw: {:?}", cdw);
    let cdw = c.encode(&msg);
    println!("cdw: {:?}", cdw);
    let dcd = c.decode(&cdw);
    println!("dcd: {:?}", dcd);
    assert_eq!(cdw, dcd);

    let cdw = c.encode(&msg);
    println!("cdw: {:?}", cdw);
    let err = Mat::weighted_vector_random(&mut rng, 8, 1);
    // let mut err = Mat::new(1, 8);
    // err[(0, 1)] = F2::one();
    println!("err: {:?}", err);
    let rcv = Mat::sum(&cdw, &err);
    println!("rcv: {:?}", rcv);
    let dcd = c.decode(&rcv);
    println!("dcd: {:?}", dcd);
    assert_eq!(cdw, dcd);

    // Problem in the maths: the rcv word which is not a codeword has a null syndrome,
    // that is the syndrome is a multiple of the goppa polynomial.
}
