use std::sync::Once;
use rand::Rng;

use mceliece::matrix;

static INIT: Once = Once::new();
const GOPPA_N_MIN: u32 = 1;
const GOPPA_N_MAX: u32 = 256 + 1;

const GOPPA_N: u32 = 0;
const GOPPA_T: u32 = 0;

// const GOPPA_N: u32 = 1024;
// const GOPPA_T: u32 = 50;

// const GOPPA_N: u32 = 97;
// const GOPPA_T: u32 = 11;

const REPEAT: u32 = 1;

/// Setup function that is only run once, even if called multiple times.
pub fn log_setup() {
    INIT.call_once(|| {
        env_logger::init();
    });
}

pub fn goppa_setup() -> (u32, u32, u32) {
    log_setup();
    let mut rng = rand::thread_rng();
    let n = match GOPPA_N {
        0 => rng.gen_range(GOPPA_N_MIN, GOPPA_N_MAX),
        value => value,
    };
    let q = n.next_power_of_two();
    let m = q.trailing_zeros();
    let t = match GOPPA_T {
        0 => {
            let t_min = if n == q { 2 } else { 1 };
            let t_max = matrix::div_ceil(n, m) + if n.is_power_of_two() { 1 } else { 0 };
            rng.gen_range(t_min, t_max)
        }
        value => value,
    };
    (m, n, t)
}
