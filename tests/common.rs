use rand::Rng;
use std::sync::Once;

static INIT: Once = Once::new();
const GOPPA_N_MIN: u32 = 2;
const GOPPA_N_MAX: u32 = 256 + 1;

const GOPPA_N: u32 = 0;
const GOPPA_T: u32 = 0;

// const GOPPA_N: u32 = 1024;
// const GOPPA_T: u32 = 50;

// const GOPPA_N: u32 = 97;
// const GOPPA_T: u32 = 11;

/// Setup function that is only run once, even if called multiple times.
pub fn log_setup() {
    INIT.call_once(|| {
        env_logger::init();
    });
}

pub fn goppa_setup() -> (u32, u32, u32) {
    log_setup();
    let mut rng = rand::thread_rng();
    let (mut m, mut n, mut t) = (1, 0, 1);
    while m * t > n {
        n = match GOPPA_N {
            0 => rng.gen_range(GOPPA_N_MIN, GOPPA_N_MAX),
            value => value,
        };
        let q = n.next_power_of_two();
        m = q.trailing_zeros();
        t = match GOPPA_T {
            0 => {
                // if n == q, the degree t of the goppa polynomial must be greater than 1.
                // Otherwise, the set L contains the polynomial root.
                let t_min = if n == q { 2 } else { 1 };
                let t_max = t_min + n / m;
                rng.gen_range(t_min, t_max)
            }
            value => value,
        };
    }
    (m, n, t)
}
