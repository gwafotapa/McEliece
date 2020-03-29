use rand::Rng;
use std::sync::Once;

static INIT: Once = Once::new();
const GOPPA_N_MIN: u32 = 3;
const GOPPA_N_MAX: u32 = 256 + 1;

const GOPPA_N: u32 = 5;
const GOPPA_T: u32 = 0;

// TODO: remove commented code

// const GOPPA_N: u32 = 1024;
// const GOPPA_T: u32 = 50;

// const GOPPA_N: u32 = 97;
// const GOPPA_T: u32 = 11;

pub fn div_ceil(n: u32, d: u32) -> u32 {
    n / d + if n % d == 0 { 0 } else { 1 }
}

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
    while m * t >= n {
        n = match GOPPA_N {
            0 => rng.gen_range(GOPPA_N_MIN, GOPPA_N_MAX),
            value => value,
        };
        let q = n.next_power_of_two();
        m = q.trailing_zeros();
        t = match GOPPA_T {
            0 => {
                let t_min = 1;
                let t_max = div_ceil(n, m);
                rng.gen_range(t_min, t_max)
            }
            value => value,
        };
        if n == q && t == 1 {
            m += 1;
        }
    }
    (m, n, t)
}
