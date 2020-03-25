use rand::{rngs::ThreadRng, Rng};
use std::{
    error::Error,
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::{Field, FiniteField, F2};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

pub use perm::Perm;
pub use rowvec::RowVec;

pub fn div_ceil(a: u32, b: u32) -> u32 {
    a / b + if a % b == 0 { 0 } else { 1 }
}

#[derive(Eq, PartialEq)]
pub struct Mat<'a, F: Eq + Field> {
    field: &'a F,
    rows: usize,
    cols: usize,
    data: Vec<F::FElt>,
}

impl<'a, F: Eq + Field> Mat<'a, F> {
    pub fn new(field: &'a F, rows: usize, cols: usize, data: Vec<F::FElt>) -> Self {
        if data.len() != rows * cols {
            panic!("Wrong dimensions");
        }
        Self {
            field,
            rows,
            cols,
            data,
        }
    }

    pub fn zero(field: &'a F, rows: usize, cols: usize) -> Self {
        Self {
            field,
            rows,
            cols,
            data: vec![field.zero(); rows * cols],
        }
    }

    pub fn field(&self) -> &'a F {
        self.field
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn data(&self) -> &Vec<F::FElt> {
        &self.data
    }

    pub fn random(rng: &mut ThreadRng, f: &'a F, n: usize, m: usize) -> Self {
        let mut mat = Self::zero(f, n, m);
        for i in 0..n {
            for j in 0..m {
                mat[(i, j)] = f.random_element(rng);
            }
        }
        mat
    }

    pub fn is_zero(&self) -> bool {
        for i in 0..self.rows {
            for j in 0..self.cols {
                if self[(i, j)] != self.field.zero() {
                    return false;
                }
            }
        }
        true
    }

    pub fn extract_cols(&self, cols: &Vec<usize>) -> Self {
        let mut res = Mat::zero(self.field(), self.rows(), cols.len());
        for j in 0..cols.len() {
            for i in 0..res.rows() {
                res[(i, j)] = self[(i, cols[j])];
            }
        }
        res
    }

    pub fn identity(f: &'a F, n: usize) -> Self {
        let mut id = Self::zero(f, n, n);
        for i in 0..n {
            id[(i, i)] = f.one();
        }
        id
    }

    // pub fn sum(&mut self, mat1: &Self, mat2: &Self) {
    //     let f = self.field;
    //     if f != mat1.field || f != mat2.field {
    //         panic!("Cannot add matrices: fields don't match");
    //     } else if self.rows != mat1.rows
    //         || self.rows != mat2.rows
    //         || self.cols != mat1.cols
    //         || self.cols != mat2.cols
    //     {
    //         panic!("Cannot add matrices: dimensions don't match");
    //     }

    //     for i in 0..self.rows {
    //         for j in 0..self.cols {
    //             self[(i, j)] = f.add(mat1[(i, j)], mat2[(i, j)]);
    //         }
    //     }
    // }

    // pub fn prod(&mut self, mat1: &Self, mat2: &Self) {
    //     let f = self.field;
    //     if f != mat1.field || f != mat2.field {
    //         panic!("Cannot add matrices: fields don't match");
    //     } else if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
    //         panic!("Cannot multiply matrices: dimensions don't match");
    //     }

    //     for i in 0..self.rows {
    //         for j in 0..self.cols {
    //             let mut sum = f.zero();
    //             for k in 0..mat1.cols {
    //                 sum = f.add(sum, f.mul(mat1[(i, k)], mat2[(k, j)]));
    //             }
    //             self[(i, j)] = sum;
    //         }
    //     }
    // }

    pub fn transpose(&self) -> Self {
        let mut t = Self::zero(self.field, self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }


    // TODO: Can I write it better ?
    fn keep_rows(&mut self, rows: &Vec<usize>) {
        let m = self.rows;
        for i in (0..m).rev() {
            if rows.contains(&i) {
                continue;
            }
            self.data.drain(i * self.cols..(i + 1) * self.cols);
            self.rows -= 1;
        }
        // let mut rows_in_descending_order = rows.clone();
        // rows_in_descending_order.sort_by(|a, b| a.cmp(b).reverse());
        // let row_to_keep_iter = rows_in_descending_order.iter();
        // let row_to_keep = row_to_keep_iter.next();
        // for i in (0..self.rows).rev() {
        //     if row_to_keep_iter.next()
        //     self.data.drain(i * self.cols..(i + 1) * self.cols);
        //     self.rows -= 1;
        // }
    }

    // TODO: Can I write it better ?
    fn remove_rows(&mut self, rows: &Vec<usize>) {
        // for i in (0..rows.len()).rev() {
        //     self.data.drain(i * self.cols..(i + 1) * self.cols);
        //     self.rows -= 1;
        // }
        let mut rows_in_descending_order = rows.clone();
        rows_in_descending_order.sort_by(|a, b| a.cmp(b).reverse());
        for i in rows_in_descending_order.iter() {
            self.data.drain(i * self.cols..(i + 1) * self.cols);
            self.rows -= 1;
        }
    }


}

mod f2;
mod gauss;
mod perm;
mod rowvec;
mod traits;
