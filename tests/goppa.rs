use log::{info, warn};

use mceliece::{finite_field::*, goppa::*, matrix::*, polynomial::*};

pub mod common;

const REPEAT: u32 = 10;

#[test]
fn goppa_f8() {
    common::log_setup();
    let f2 = &F2 {};
    let f8 = &F2m::generate(8);
    let p = Poly::support(f8, &[2, 1, 0]);
    let c = Goppa::new(p, [0, 1, 2, 3, 4, 5, 6, 7].to_vec());
    info!("{}", c);

    let x = c.parity_check_x();
    assert_eq!(x, Mat::new(f8, 2, 2, [1, 0, 1, 1].to_vec()));

    let y = c.parity_check_y();
    assert_eq!(
        y,
        Mat::new(
            f8,
            2,
            8,
            [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2, 3, 4, 5, 6, 7].to_vec()
        )
    );

    let h = c.parity_check_matrix(f2);
    info!("h:{}", h);

    let (g, _) = Goppa::<F2m>::generator_from_parity_check(&h);
    info!("g:{}", g);

    let z = &h * &g.transpose();
    assert!(z.is_zero());

    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, 2);
    info!("Message:{}", msg);

    let cdw = c.encode(&msg);
    info!("cdw:{}", cdw);

    let err = RowVec::random_with_weight(&mut rng, f2, 8, 2);
    info!("err:{}", err);

    let rcv = &cdw + err;
    info!("rcv:{}", rcv);

    let dcdw = c.decode(&rcv);
    info!("dcdw:{}", dcdw);

    assert_eq!(cdw, dcdw);
}

#[test]
fn goppa_f1024() {
    common::log_setup();
    let f2 = &F2 {};
    let f1024 = &F2m::generate(1024);
    let mut rng = rand::thread_rng();
    let c = Goppa::random(&mut rng, f1024, 30, 2);
    info!("{}", c);

    let h = c.parity_check_matrix(f2);
    info!("parity check matrix:{}", h);

    let (g, _) = Goppa::<F2m>::generator_from_parity_check(&h);
    assert!((&h * &g.transpose()).is_zero());

    let cdw = RowVec::new(f2, g.data()[0..30].to_vec());
    info!("codeword: {}", cdw);
    let syndrome = Goppa::syndrome_from_xyz(&h, &cdw);

    assert!(syndrome.is_zero());

    let err = RowVec::random_with_weight(&mut rng, f2, 30, 2);
    info!("error:{}", err);

    let rcv = &cdw + &err;
    info!("received word:{}", rcv);

    let dcdw = c.decode(&rcv);
    info!("decoded word:{}", dcdw);

    assert_eq!(cdw, dcdw);
}

#[test]
fn goppa_f256() {
    common::log_setup();
    let f2 = &F2 {};
    let f256 = &F2m::generate(256);
    let mut g = Poly::support(f256, &[22, 17, 15, 12, 5]);
    g[0] = f256.exp(78);
    let mut f256_elts = Vec::with_capacity(256);
    for i in 0..256 {
        f256_elts.push(i);
    }
    let c = Goppa::new(g, f256_elts);
    info!("{}", c);

    let n = f256.order();
    let m = f256.characteristic_exponent() as usize;
    let t = 22;
    let mut rng = rand::thread_rng();
    let msg = RowVec::random(&mut rng, f2, n - m * t);
    info!("msg:{}", msg);

    let cdw = c.encode(&msg);
    info!("cdw:{}", cdw);

    let err = RowVec::random_with_weight(&mut rng, f2, n, t);
    info!("err:{}", err);

    let rcv = &cdw + err;
    info!("rcv:{}", rcv);

    let dcdw = c.decode(&rcv);
    info!("dcdw:{}", dcdw);

    assert_eq!(cdw, dcdw);
}

#[test]
fn goppa_f128() {
    common::log_setup();
    let m = 7;
    let n = 1 << m;
    let t = 10;
    let k = n - m * t;
    let f2 = &F2 {};
    let f128 = &F2m::generate(1 << m);
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f128, n, t);
    info!("{}", goppa);

    let msg = RowVec::random(&mut rng, f2, k);
    info!("msg:{}", msg);

    let cdw = goppa.encode(&msg);
    info!("cdw:{}", cdw);

    let err = RowVec::random_with_weight(&mut rng, f2, n, t);
    info!("err:{}", err);

    let rcv = &cdw + err;
    info!("rcv:{}", rcv);

    let dcdw = goppa.decode(&rcv);
    info!("dcdw:{}", dcdw);

    assert_eq!(cdw, dcdw);
}

#[test]
fn goppa_random() {
    let (m, n, t) = common::goppa_setup();
    info!("F2m=F{}, n={}, t={}", 1 << m, n, t);

    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f2m, n, t);
    info!("{}", goppa);

    let xyz = goppa.parity_check_xyz();
    info!("XYZ:{}", xyz);

    let h = Goppa::parity_check_from_xyz(&xyz, f2);
    info!("Parity-check matrix H:{}", h);

    let (g, _) = Goppa::<F2m>::generator_from_parity_check(&h);
    info!("Generator matrix G:{}", g);

    let k = g.rows();
    info!("Code dimension k = {}", k);

    let msg = RowVec::random(&mut rng, f2, k);
    info!("message:{}", msg);

    let cdw = &msg * &g;
    info!("codeword:{}", cdw);

    let err = RowVec::random_with_weight(&mut rng, f2, n, t);
    info!("error:{}", err);

    let rcv = &cdw + &err;
    info!("received word:{}", rcv);

    let dcdw = goppa.xyz_decode(&xyz, &rcv);
    info!("decoded codeword:{}", dcdw);

    assert_eq!(cdw, dcdw);
}

#[test]
fn goppa_repeat() {
    common::log_setup();
    warn!("Series of {} random encoding - decoding", REPEAT);

    for _r in 0..REPEAT {
        let (m, n, t) = common::goppa_setup();
        warn!(
            "Encode - Decode #{:02}: F2m=F{:<4}, n={:<4}, t={:<4}",
            _r,
            1 << m,
            n,
            t
        );
        goppa_repeated(m, n, t);
    }
}

fn goppa_repeated(m: u32, n: usize, t: usize) {
    let f2 = &F2 {};
    let f2m = &F2m::generate(1 << m);
    let mut rng = rand::thread_rng();
    let goppa = Goppa::random(&mut rng, f2m, n, t);
    info!("{}", goppa);

    let xyz = goppa.parity_check_xyz();
    info!("XYZ:{}", xyz);

    let h = Goppa::parity_check_from_xyz(&xyz, f2);
    info!("Parity-check matrix H:{}", h);

    let (g, _) = Goppa::<F2m>::generator_from_parity_check(&h);
    info!("Generator matrix G:{}", g);

    let k = g.rows();
    info!("Code dimension k = {}", k);

    let msg = RowVec::random(&mut rng, f2, k);
    info!("message:{}", msg);

    let cdw = &msg * &g;
    info!("codeword:{}", cdw);

    let err = RowVec::random_with_weight(&mut rng, f2, n, t);
    info!("error:{}", err);

    let rcv = &cdw + &err;
    info!("received word:{}", rcv);

    let dcdw = goppa.xyz_decode(&xyz, &rcv);
    info!("decoded codeword:{}", dcdw);

    assert_eq!(cdw, dcdw);
}
