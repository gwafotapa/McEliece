pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}

pub trait Inverse {
    // fn is_invertible(&self) -> bool;

    fn inv(&self) -> Option<Self>
    where
        Self: Sized;
}
