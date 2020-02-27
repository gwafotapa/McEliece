extern crate mceliece;

use mceliece::finite_field::*;
use mceliece::finite_field_8::*;
use mceliece::goppa::*;
use mceliece::matrix::*;
use mceliece::polynomial::*;

#[test]
fn goppa_f8() {
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
    let y = c.parity_check_matrix_y();
    println!("{:?}", y);
    let z = c.parity_check_matrix_z();
    println!("{:?}", z);
    let h = c.parity_check_matrix();
    println!("{:?}", h);

    let g = c.generator_matrix();
    println!("{:?}", g);

    let z = Mat::new(h.rows(), g.rows());
    let mut m = Mat::new(h.rows(), g.rows());
    m.as_prod(&h, &g.transpose());
    assert_eq!(m, z);

    let h_bin_form = h.binary_form(3);
    println!("{:?}", h_bin_form);

    // book page 342 add test to confirm matrix h binary form
    // Add test implementing the next code and checking the generator matrix
    // Add decoding test
}
