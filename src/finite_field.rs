pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}

pub trait Inv {
    fn inv(self) -> Option<Self>
    where
        Self: Sized;
}

pub trait Exp {
    fn exp(n: usize) -> Self;
}

pub trait Log {
    fn log(self) -> Option<usize>;
}




// struct FFElt(u16);

// struct FF {
//     order: usize,
//     p: usize,
//     m: usize,
//     exp: Vec<FFElt>,
//     log: Vec<usize>,
// }

// impl Add for FFElt {
//     type Output = FFElt;

//     fn add(self, other: FFElt) -> FFElt {
//         FFElt(self.0 ^ other.0)
//     }
// }

// impl Sub for FFElt {
//     type Output = FFElt;

//     fn add(self, other: FFElt) -> FFElt {
//         FFElt(self.0 ^ other.0)
//     }
// }
