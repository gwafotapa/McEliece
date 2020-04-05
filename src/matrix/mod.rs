//! Matrices on a field

use std::rc::Rc;

use crate::finite_field::{Field, FieldTrait};

pub use perm::Perm;
pub use rowvec::RowVec;

/// Matrix with coefficients in a field F
#[derive(Eq, PartialEq)]
pub struct Mat<F>
where
    F: FieldTrait,
{
    field: Rc<F>,
    rows: usize,
    cols: usize,
    data: Vec<F::FieldElement>,
}

impl<F> Mat<F>
where
    F: FieldTrait,
{
    /// Creates a new matrix
    ///
    /// The vector data holds the matrix coefficients, row by row.
    ///
    /// # Panics
    ///
    /// Panics if the matrix is empty or if there are not exactly rows * cols coefficients.
    pub fn new(field: Field<F>, rows: usize, cols: usize, data: Vec<F::FieldElement>) -> Self {
        if rows == 0 || cols == 0 || data.is_empty() {
            panic!("Empty matrix");
        } else if data.len() != rows * cols {
            panic!("Dimensions do not match the number of coefficients");
        }
        let field = match field {
            Field::Some(f) => Rc::clone(f),
            Field::Parameters(p) => Rc::new(F::generate(p)),
        };
        Self {
            field,
            rows,
            cols,
            data,
        }
    }

    /// Creates a new matrix whose coefficients are all zero
    ///
    /// # Panics
    ///
    /// Panics if the number of either rows or columns is zero.
    pub fn zero(field: Field<F>, rows: usize, cols: usize) -> Self {
        if rows == 0 || cols == 0 {
            panic!("Empty matrix");
        }
        let field = match field {
            Field::Some(f) => Rc::clone(f),
            Field::Parameters(p) => Rc::new(F::generate(p)),
        };
        let data = vec![field.zero(); rows * cols];
        Self {
            field,
            rows,
            cols,
            data,
        }
    }

    // TODO: Should I return a Rc instead of a &Rc now ?
    pub fn field(&self) -> &Rc<F> {
        &self.field
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn data(&self) -> &Vec<F::FieldElement> {
        &self.data
    }

    pub fn random(field: Field<F>, n: usize, m: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut mat = Self::zero(field, n, m);
        for i in 0..n {
            for j in 0..m {
                mat[(i, j)] = mat.field.random_element(&mut rng);
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

    /// Creates a new matrix by taking the chosen columns in the given order
    pub fn extract_cols(&self, cols: &Vec<usize>) -> Self {
        let mut res = Mat::zero(Field::Some(self.field()), self.rows(), cols.len());
        for j in 0..cols.len() {
            for i in 0..res.rows() {
                res[(i, j)] = self[(i, cols[j])];
            }
        }
        res
    }

    pub fn identity(field: Field<F>, n: usize) -> Self {
        let mut id = Self::zero(field, n, n);
        for i in 0..n {
            id[(i, i)] = id.field.one();
        }
        id
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::zero(Field::Some(self.field()), self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }

    /// Keep only given rows
    ///
    /// # Panics
    ///
    /// Panics if the vector contains an index that does not match any row of the matrix
    pub fn keep_rows(&mut self, rows: &Vec<usize>) {
        let mut rrows = rows.clone();
        rrows.sort_by(|a, b| a.cmp(b).reverse());
        if *rrows.first().unwrap() >= self.rows {
            panic!("invalid row index");
        }
        let mut iter = rrows.into_iter();
        let mut row = iter.next().unwrap_or(self.rows);
        for i in (0..self.rows).rev() {
            if i == row {
                row = iter.next().unwrap_or(self.rows);
                continue;
            }
            self.data.drain(i * self.cols..(i + 1) * self.cols);
            self.rows -= 1;
        }
    }

    /// Remove the given rows from the matrix
    ///
    /// # Panics
    ///
    /// Panics if the vector contains an index that does not match any row of the matrix
    pub fn remove_rows(&mut self, rows: &Vec<usize>) {
        let mut rrows = rows.clone();
        rrows.sort_by(|a, b| a.cmp(b).reverse());
        if *rrows.first().unwrap() >= self.rows {
            panic!("invalid row index");
        }
        for i in rrows.iter() {
            self.data.drain(i * self.cols..(i + 1) * self.cols);
            self.rows -= 1;
        }
    }
}

pub mod f2;
pub mod gauss;
pub mod perm;
pub mod rowvec;
pub mod traits;
