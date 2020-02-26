extern crate rand;

use crate::finite_field::{Inv, One, Zero};
use rand::{distributions, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

// Type T must represent an element from a field, meaning all elements except 0 are inversible.
#[derive(Clone, Eq, PartialEq)]
pub struct Mat<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T> Index<(usize, usize)> for Mat<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for Mat<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}

impl<T> fmt::Debug for Mat<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();

        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                s.push_str(&self[(i, j)].to_string());
                s.push(' ');
            }
            s.push_str(&self[(i, self.cols - 1)].to_string());
            s.push('\n');
        }

        write!(f, "\n{}", s)
    }
}

impl<T> Mat<T>
where
    T: fmt::Display
        + Copy
        + Eq
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Neg<Output = T>
        + Inv
        + AddAssign
        + SubAssign
        + MulAssign,
{
    pub fn new(n: usize, m: usize) -> Mat<T> {
        let mut mat = Mat {
            rows: n,
            cols: m,
            data: Vec::<T>::with_capacity(n * m),
        };
        mat.data.resize(n * m, T::zero());

        mat
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn data(&self) -> Vec<T> {
        self.data.clone()
    }

    // pub fn set(&mut self, row: usize, col: usize, val: T) {
    //     self[(row, col)] = val;
    // }

    // pub fn get(&self, row: usize, col: usize) -> T {
    //     self[(row, col)]
    // }

    // pub fn print(&self) {
    //     for i in 0..self.rows {
    //         for j in 0..self.cols - 1 {
    //             print!("{} ", self[(i, j)]);
    //         }
    //         println!("{}", self[(i, self.cols - 1)]);
    //     }
    
    // //   println!("{:?}", self);
    // }

    pub fn random(rng: &mut rand::rngs::ThreadRng, n: usize, m: usize) -> Mat<T>
    where
        distributions::Standard: distributions::Distribution<T>,
    {
        let mut mat = Mat::new(n, m);
        for i in 0..n {
            for j in 0..m {
                mat[(i, j)] = rng.gen::<T>();
            }
        }

        mat
    }

    pub fn is_permutation(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }

        let n = self.rows;
        let mut cols = Vec::with_capacity(n); // indices of columns containing a '1'
        cols.resize(n, 0);
        for i in 0..n {
            // loop on rows
            let mut has_one = false;
            for j in 0..n {
                let val = self[(i, j)];
                if val == T::zero() {
                    continue;
                } else if val == T::one() {
                    if cols[j] == 1 {
                        return false;
                    }
                    cols[j] = 1;
                    has_one = true;
                } else {
                    return false;
                }
            }
            if !has_one {
                // row contains zeroes only
                return false;
            }
        }

        true
    }

    pub fn permutation_random(rng: &mut rand::rngs::ThreadRng, n: usize) -> Mat<T> {
        let mut mat = Mat::new(n, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..n {
            // Draw a random column index j to set the '1' on this row
            let nbr = rng.gen_range(0, n - i);
            mat[(i, cols[nbr])] = T::one();

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }

        mat
    }

    pub fn identity(n: usize) -> Mat<T> {
        let mut id = Mat::new(n, n);
        for i in 0..n {
            id[(i, i)] = T::one();
        }

        id
    }

    pub fn is_invertible(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }

        self.rank() == self.rows
    }

    pub fn swap_rows(&mut self, row1: usize, row2: usize) {
        if row1 == row2 {
            return;
        }

        let slice = &mut self.data[row1 * self.cols..(row1 + 1) * self.cols].to_vec();
        slice.swap_with_slice(&mut self.data[row2 * self.cols..(row2 + 1) * self.cols]);
        slice.swap_with_slice(&mut self.data[row1 * self.cols..(row1 + 1) * self.cols]);
    }

    pub fn swap_cols(&mut self, col1: usize, col2: usize) {
        if col1 == col2 {
            return;
        }

        for i in 0..self.rows {
            let tmp = self[(i, col1)];
            self[(i, col1)] = self[(i, col2)];
            self[(i, col2)] = tmp;
        }
    }

    pub fn combine_rows(&mut self, row1: usize, lambda: T, row2: usize) {
        for j in 0..self.cols {
            self[(row1, j)] = self[(row1, j)] + lambda * self[(row2, j)];
        }
    }

    // Inv computation via gaussian elimination
    // See https://en.wikipedia.org/wiki/Gaussian_elimination
    pub fn inverse(&self) -> Option<Mat<T>> {
        if self.rows != self.cols {
            return None;
        }

        let n = self.rows;
        let mut mat = self.clone();
        let mut inv: Mat<T> = Mat::identity(n);
        let mut p = 0; // pivot's row and pivot's column

        while p < n {
            // Find pivot
            let mut i = p;
            while i < n && mat[(i, p)] == T::zero() {
                i += 1;
            }
            if i == n {
                return None;
            }

            // Swap rows and mimic operation on matrix 'inv'
            mat.swap_rows(i, p);
            inv.swap_rows(i, p);

            // Normalize pivot's row: L(p) = pivot^-1 * L(p)
            let pivot_inv = mat[(p, p)].inv().unwrap();
            for j in p + 1..n {
                // first p+1 columns are zero
                mat[(p, j)] = pivot_inv * mat[(p, j)];
            }
            mat[(p, p)] = T::one();

            // Mimic normalization on matrix 'inv'
            for j in 0..n {
                inv[(p, j)] = pivot_inv * inv[(p, j)];
            }

            // Adjust all rows below pivot's row: L(k) = L(k) - c(k,p) * L(p)
            // where L(k) is the kth row and c(k,p) is the coefficient [k,p] of our matrix
            for k in p + 1..n {
                if mat[(k, p)] == T::zero() {
                    continue;
                }

                let lambda = mat[(k, p)];

                mat[(k, p)] = T::zero();
                for l in p + 1..n {
                    // first p+1 columns are zero
                    mat[(k, l)] = mat[(k, l)] - lambda * mat[(p, l)];
                }

                // Mimic operation on matrix 'inv'
                for l in 0..n {
                    inv[(k, l)] = inv[(k, l)] - lambda * inv[(p, l)];
                }
            }

            p += 1;
        }

        // Matrix 'mat' is now in triangular form

        for j in (0..n).rev() {
            for i in (0..j).rev() {
                if mat[(i, j)] == T::zero() {
                    continue;
                }

                // Perform the row operation to set c(i, j) to 0:
                // L(i) = L(i) - c(i, j) * L(j)
                // We don't actually need to operate on the original matrix here.
                // Mimic the row operation on matrix 'inv'.
                for l in 0..n {
                    inv[(i, l)] = inv[(i, l)] - mat[(i, j)] * inv[(j, l)];
                }
            }
        }

        Some(inv)
    }

    pub fn invertible_random(rng: &mut rand::rngs::ThreadRng, n: usize) -> Mat<T>
    where
        distributions::Standard: distributions::Distribution<T>,
    {
        let mut mat = Mat::new(n, n);
        let mut i = 0;
        while i < n {
            // Fill line i at random
            for j in 0..n {
                mat[(i, j)] = rng.gen::<T>();
            }

            if mat.rank() == i + 1 {
                i += 1;
            }
        }

        mat
    }

    pub fn row_echelon_form(&mut self) -> usize {
        let n = self.rows;
        let m = self.cols;
        let mut rank = 0;
        let mut row_pivot = 0;
        let mut col_pivot = 0;

        while row_pivot < self.rows && col_pivot < self.cols {
            // Find pivot
            let mut i = row_pivot;
            while i < n && self[(i, col_pivot)] == T::zero() {
                i += 1;
            }
            if i == n {
                col_pivot += 1;
                continue;
            }
            rank += 1;

            // Swap rows
            // if i != row_pivot {
            self.swap_rows(i, row_pivot);

            // for j in 0..m {
            //     let tmp = self.get(row_pivot, j);
            //     self.set(row_pivot, j, self.get(i, j));
            //     self.set(i, j, tmp);
            // }
            // }

            // Normalize pivot's row
            let pivot_inv = self[(row_pivot, col_pivot)].inv().unwrap();
            for j in col_pivot + 1..m {
                self[(row_pivot, j)] = pivot_inv * self[(row_pivot, j)];
            }
            self[(row_pivot, col_pivot)] = T::one();

            // Adjust all rows below pivot's row
            for k in row_pivot + 1..n {
                if self[(k, col_pivot)] == T::zero() {
                    continue;
                }
                for l in col_pivot + 1..m {
                    self[(k, l)] = self[(k, l)] - self[(k, col_pivot)] * self[(row_pivot, l)];
                }
                self[(k, col_pivot)] = T::zero();
            }

            row_pivot += 1;
            col_pivot += 1;
        }

        rank
    }

    pub fn reduced_row_echelon_form(&mut self) -> usize {
        let n = self.rows;
        let m = self.cols;
        let rank = self.row_echelon_form(); // note that all pivots are 1

        for row_pivot in (0..n).rev() {
            // Find the pivot on this row if any
            let mut col_pivot = 0;
            while col_pivot < m && self[(row_pivot, col_pivot)] == T::zero() {
                col_pivot += 1;
            }
            if col_pivot == m {
                continue;
            }

            // Eliminate all non zero elements in the pivot's column
            for i in (0..row_pivot).rev() {
                if self[(i, col_pivot)] == T::zero() {
                    continue;
                }

                for k in col_pivot..m {
                    self[(i, k)] =
                        self[(i, k)] - self[(i, col_pivot)] * self[(row_pivot, col_pivot)];
                }
            }
        }

        rank
    }

    pub fn rank(&self) -> usize {
        let mut mat = self.clone();
        mat.row_echelon_form()
    }

    pub fn sum(mat1: &Mat<T>, mat2: &Mat<T>) -> Mat<T> {
        let mut sum = Mat::new(mat1.rows, mat1.cols);
        sum.as_sum(mat1, mat2);
        sum
    }

    pub fn as_sum(&mut self, mat1: &Mat<T>, mat2: &Mat<T>) {
        if self.rows != mat1.rows
            || self.rows != mat2.rows
            || self.cols != mat1.cols
            || self.cols != mat2.cols
        {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = mat1[(i, j)] + mat2[(i, j)];
            }
        }
    }

    pub fn as_prod(&mut self, mat1: &Mat<T>, mat2: &Mat<T>) {
        if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = T::zero();
                for k in 0..mat1.cols {
                    sum += mat1[(i, k)] * mat2[(k, j)];
                }
                self[(i, j)] = sum;
            }
        }
    }

    // Generates a random row vector of length n and weight t
    pub fn weighted_vector_random(rng: &mut rand::rngs::ThreadRng, n: usize, t: usize) -> Mat<T>
    where
        distributions::Standard: distributions::Distribution<T>,
    {
        let mut vec = Mat::new(1, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..t {
            // Draw a random column index
            let nbr = rng.gen_range(0, n - i);
            let mut elt = rng.gen::<T>();
            while elt == T::zero() {
                elt = rng.gen::<T>();
            }
            vec[(0, cols[nbr])] = elt;

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }

        vec
    }

    pub fn weight(&self) -> Option<usize> {
        if self.rows != 1 {
            return None;
        }

        let mut cnt = 0;
        for j in 0..self.cols {
            if self[(0, j)] != T::zero() {
                cnt += 1;
            }
        }

        Some(cnt)
    }

    // Compute, if possible, (U, H, P) with U invertible, H standard form and P permutation
    // such that self = UHP
    pub fn standard_form(&self) -> Option<(Mat<T>, Mat<T>, Mat<T>)> {
        let m = self.rows;
        let n = self.cols;
        if m > n {
            return None;
        }
        let mut u = Mat::identity(m);
        let mut h = self.clone();
        let mut p = Mat::identity(n);
        let mut col = n; // index of the column to check for a pivot

        // j is the index of the column to "standardize":
        // The first iteration sets a 1 at the last position (m-1) of column n-1.
        // The second iteration sets a 1 at position m-2 of column n-2.
        // ...
        // The last iteration sets a 1 at position 0 of column n-m.
        for j in (n - m..n).rev() {
            // Among the remaining columns, select one with a pivot
            let mut pivot = false;
            let mut row_pivot = 0;
            let mut col_pivot = 0;
            while !pivot && col != 0 {
                col -= 1;

                // Check column 'col' for a pivot
                for row in (0..j + m - n + 1).rev() {
                    if h[(row, col)] != T::zero() {
                        pivot = true;
                        row_pivot = row;
                        col_pivot = col;
                        break;
                    }
                }
            }

            if !pivot {
                return None;
            }

            // Put pivot column in the adequate position and update P
            h.swap_cols(j, col_pivot);
            p.swap_cols(j, col_pivot);

            // Put pivot row in the adequate position and update U
            h.swap_rows(j + m - n, row_pivot);
            u.swap_rows(j + m - n, row_pivot);

            // Pivot is now at (j+m-n, j)

            // Multiply pivot row by pivot^-1 and update U
            let pivot_inv = h[(j + m - n, j)].inv().unwrap();
            for k in 0..n {
                h[(j + m - n, k)] = pivot_inv * h[(j + m - n, k)];
            }
            for k in 0..m {
                u[(j + m - n, k)] = pivot_inv * u[(j + m - n, k)];
            }

            // Nullify the rest of the column and update matrix U accordingly

            // for i in 0..j + m - n {
            //     if h.get(i, j) != T::zero() {
            //         h.combine_rows(i, h.get(i, j), j + m - n);
            //         u.combine_rows(i, h.get(i, j), j + m - n);
            //     }
            // }
            // for i in j + m - n + 1..m {
            //     if h.get(i, j) != T::zero() {
            //         h.combine_rows(i, h.get(i, j), j + m - n);
            //         u.combine_rows(i, h.get(i, j), j + m - n);
            //     }
            // }

            for i in 0..m {
                if h[(i, j)] != T::zero() && i != j + m - n {
                    let lambda = -h[(i, j)];
                    h.combine_rows(i, lambda, j + m - n);
                    u.combine_rows(i, lambda, j + m - n);
                }
            }
        }

        Some((u, h, p))
    }

    pub fn is_standard_form(&self) -> bool {
        let m = self.rows;
        let n = self.cols;
        if m > n {
            return false;
        }
        for i in 0..m {
            for j in n - m..n {
                if n + i == m + j {
                    if self[(i, j)] != T::one() {
                        return false;
                    }
                } else {
                    if self[(i, j)] != T::zero() {
                        return false;
                    }
                }
            }
        }

        true
    }

    pub fn transpose(&self) -> Mat<T> {
        let mut t = Mat::new(self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }

    pub fn prod(a: &Mat<T>, b: &Mat<T>) -> Mat<T> {
        let mut p = Mat::new(a.rows, b.cols);
        p.as_prod(a, b);
        p
    }
}
