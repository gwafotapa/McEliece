extern crate mceliece;

use mceliece::finite_field::*;
use mceliece::finite_field_8::*;
use mceliece::polynomial::*;
use mceliece::goppa::*;
use mceliece::matrix::*;

#[test]
fn goppa_f8() {
    let mut p = Poly::<F8>::x_n(2);
    p.set(0, F8::one());
    p.set(1, F8::one());
    println!("{}\n", p.to_str());
    let c: Goppa<F8> = Goppa::new(&p, [F8(0), F8(1), F8(2), F8(3), F8(4), F8(5), F8(6), F8(7)]).unwrap();
    let y = c.parity_check_y();
    y.print();
    println!();
    let z = c.parity_check_z();
    z.print();
    println!();
    let h = c.parity_check();
    h.print();
    println!();

    let g = c.generator();
    g.print();
    println!();

    let z = Mat::new(h.rows(), g.rows());
    let mut m = Mat::new(h.rows(), g.rows());
    m.mul(&h, &g.transpose());
    assert_eq!(m, z);
}
