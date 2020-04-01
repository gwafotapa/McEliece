// use rand::{rngs::ThreadRng, Rng};
use std::{convert::TryInto, error::Error};

use super::{Mat, F2};
use crate::finite_field::F2FiniteExtension;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

// TODO: clean file from commentaries
impl<'a> Mat<'a, F2> {
    // pub fn add_rows(&mut self, row1: usize, row2: usize) {
    //     let f = self.field;
    //     for j in 0..self.cols {
    //         self[(row1, j)] = f.add(self[(row1, j)], self[(row2, j)]);
    //     }
    // }

    // // Inv computation via gaussian elimination
    // // See https://en.wikipedia.org/wiki/Gaussian_elimination
    // pub fn inverse_f2(&self) -> Option<Self> {
    //     if self.rows != self.cols {
    //         return None;
    //     }

    //     let f = self.field;
    //     let n = self.rows;
    //     let mut mat = self.clone();
    //     let mut inv: Self = Mat::identity(f, n);
    //     let mut p = 0; // pivot's row and pivot's column

    //     while p < n {
    //         // Find pivot
    //         let mut i = p;
    //         while i < n && mat[(i, p)] == f.zero() {
    //             i += 1;
    //         }
    //         if i == n {
    //             return None;
    //         }

    //         // Swap rows and mimic operation on matrix 'inv'
    //         mat.swap_rows(i, p);
    //         inv.swap_rows(i, p);

    //         // Adjust all rows below pivot's row: L(k) = L(k) - c(k,p) * L(p)
    //         // where L(k) is the kth row and c(k,p) is the coefficient [k,p] of our matrix
    //         for k in p + 1..n {
    //             if mat[(k, p)] == f.zero() {
    //                 continue;
    //             }

    //             mat[(k, p)] = f.zero();
    //             for l in p + 1..n {
    //                 // first p+1 columns are zero
    //                 mat[(k, l)] = f.add(mat[(k, l)], mat[(p, l)]);
    //             }

    //             // Mimic operation on matrix 'inv'
    //             for l in 0..n {
    //                 inv[(k, l)] = f.add(inv[(k, l)], inv[(p, l)]);
    //             }
    //         }

    //         p += 1;
    //     }

    //     // Matrix 'mat' is now in triangular form

    //     for j in (0..n).rev() {
    //         for i in (0..j).rev() {
    //             // Perform the row operation to set c(i, j) to 0:
    //             // L(i) = L(i) - c(i, j) * L(j)
    //             // We don't actually need to operate on the original matrix here.
    //             // Mimic the row operation on matrix 'inv'.
    //             if mat[(i, j)] == f.zero() {
    //                 continue;
    //             }

    //             for l in 0..n {
    //                 inv[(i, l)] = f.add(inv[(i, l)], inv[(j, l)]);
    //             }
    //         }
    //     }
    //     Some(inv)
    // }

    // /* Generate a random matrix and put it in standard form with a diagonal of 1.
    //  * Keep track of the applied transformations via a matrix u.
    //  * Return u is as our random invertible matrix. */
    // pub fn invertible_random_f2(rng: &mut ThreadRng, f2: &'a F2, n: usize) -> Self {
    //     let mut mat = Mat::random(rng, f2, n, n);
    //     let mut u = Mat::identity(f2, n);
    //     // Loop on columns
    //     for j in 0..n {
    //         // Find a pivot in column j
    //         let mut i = j;
    //         while i < n && mat[(i, j)] == f2.zero() {
    //             i += 1;
    //         }
    //         // If column j has no pivot, create it
    //         if i == n {
    //             i = rng.gen_range(j, n);
    //             mat[(i, j)] = f2.one();
    //         }
    //         // Put pivot in position (j, j) and mirror operation on matrix u
    //         mat.swap_rows(i, j);
    //         u.swap_rows(i, j);
    //         // Zero coefficients under the pivot and mirror operation on matrix u
    //         for i in j + 1..n {
    //             if mat[(i, j)] == f2.zero() {
    //                 continue;
    //             }
    //             mat[(i, j)] = f2.zero();
    //             for k in j + 1..n {
    //                 mat[(i, k)] = f2.add(mat[(i, k)], mat[(j, k)]);
    //             }
    //             u.add_rows(i, j);
    //         }
    //     }
    //     u
    // }

    // pub fn standard_parity_check_equivalent_f2(&self) -> (Self, Perm) {
    //     let f = self.field;
    //     let m = self.rows;
    //     let n = self.cols;
    //     if m > n {
    //         panic!(
    //             "Matrix must have at most as many rows as columns for parity-check decomposition"
    //         );
    //     }
    //     let mut h = self.clone();
    //     let mut p = Perm::identity(n);
    //     let mut col = n; // index of the column to check for a pivot

    //     // j is the index of the column to "standardize":
    //     // The first iteration sets a 1 at the last position (m-1) of column n-1.
    //     // The second iteration sets a 1 at position m-2 of column n-2.
    //     // ...
    //     // The last iteration sets a 1 at position 0 of column n-m.
    //     for j in (n - m..n).rev() {
    //         // Among the remaining columns, select one with a pivot
    //         let mut pivot = false;
    //         let mut row_pivot = 0;
    //         let mut col_pivot = 0;
    //         while !pivot && col != 0 {
    //             col -= 1;

    //             // Check column 'col' for a pivot
    //             for row in (0..j + m - n + 1).rev() {
    //                 if h[(row, col)] != f.zero() {
    //                     pivot = true;
    //                     row_pivot = row;
    //                     col_pivot = col;
    //                     break;
    //                 }
    //             }
    //         }

    //         if !pivot {
    //             h.remove_rows(&(0..j - (n - m) + 1).collect());
    //             return (h, p);
    //         }

    //         // Put pivot column in the adequate position and update P
    //         h.swap_cols(j, col_pivot);
    //         p.swap(j, col_pivot);

    //         // Put pivot row in the adequate position
    //         h.swap_rows(j + m - n, row_pivot);

    //         // Pivot is now at (j+m-n, j)

    //         // Nullify the rest of the column
    //         for i in 0..m {
    //             if h[(i, j)] != f.zero() && i != j + m - n {
    //                 h.add_rows(i, j + m - n);
    //             }
    //         }
    //     }
    //     (h, p)
    // }

    // pub fn row_echelon_form_f2(&mut self) -> usize {
    //     let f = self.field;
    //     let n = self.rows;
    //     let m = self.cols;
    //     let mut row_pivot = 0;
    //     let mut col_pivot = 0;
    //     let mut rank = 0;

    //     while row_pivot < self.rows && col_pivot < self.cols {
    //         // Find pivot
    //         let mut i = row_pivot;
    //         while i < n && self[(i, col_pivot)] == f.zero() {
    //             i += 1;
    //         }
    //         if i == n {
    //             col_pivot += 1;
    //             continue;
    //         }
    //         rank += 1;
    //         self.swap_rows(i, row_pivot);

    //         // Adjust all rows below pivot's row
    //         for k in row_pivot + 1..n {
    //             if self[(k, col_pivot)] == f.zero() {
    //                 continue;
    //             }
    //             for l in col_pivot + 1..m {
    //                 self[(k, l)] = f.add(self[(k, l)], self[(row_pivot, l)]);
    //             }
    //             self[(k, col_pivot)] = f.zero();
    //         }

    //         row_pivot += 1;
    //         col_pivot += 1;
    //     }
    //     rank
    // }

    // pub fn reduced_row_echelon_form_f2(&mut self) -> usize {
    //     let f = self.field;
    //     let n = self.rows;
    //     let m = self.cols;
    //     let rank = self.row_echelon_form_f2(); // note that all pivots are 1

    //     for row_pivot in (0..n).rev() {
    //         // Find the pivot on this row if any
    //         let mut col_pivot = 0;
    //         while col_pivot < m && self[(row_pivot, col_pivot)] == f.zero() {
    //             col_pivot += 1;
    //         }
    //         if col_pivot == m {
    //             continue;
    //         }

    //         // Eliminate all non zero elements in the pivot's column
    //         for i in (0..row_pivot).rev() {
    //             if self[(i, col_pivot)] == f.zero() {
    //                 continue;
    //             }

    //             for k in col_pivot..m {
    //                 self[(i, k)] = f.add(self[(i, k)], self[(row_pivot, col_pivot)]);
    //             }
    //         }
    //     }
    //     rank
    // }

    // pub fn rank_f2(&self) -> usize {
    //     let mut mat = self.clone();
    //     mat.row_echelon_form_f2()
    // }

    // // Compute, if possible, (U, S, P) with U invertible, S standard form and P permutation
    // // such that S = U * self * P
    // pub fn standard_form_f2(&self) -> Option<(Self, Self, Perm)> {
    //     let f = self.field;
    //     let m = self.rows;
    //     let n = self.cols;
    //     if m > n {
    //         return None;
    //     }
    //     let mut u = Mat::identity(f, m);
    //     let mut h = self.clone();
    //     let mut p = Perm::identity(n);
    //     let mut col = n; // index of the column to check for a pivot

    //     // j is the index of the column to "standardize":
    //     // The first iteration sets a 1 at the last position (m-1) of column n-1.
    //     // The second iteration sets a 1 at position m-2 of column n-2.
    //     // ...
    //     // The last iteration sets a 1 at position 0 of column n-m.
    //     for j in (n - m..n).rev() {
    //         // Among the remaining columns, select one with a pivot
    //         let mut pivot = false;
    //         let mut row_pivot = 0;
    //         let mut col_pivot = 0;
    //         while !pivot && col != 0 {
    //             col -= 1;

    //             // Check column 'col' for a pivot
    //             for row in (0..j + m - n + 1).rev() {
    //                 if h[(row, col)] != f.zero() {
    //                     pivot = true;
    //                     row_pivot = row;
    //                     col_pivot = col;
    //                     break;
    //                 }
    //             }
    //         }

    //         if !pivot {
    //             return None;
    //         }

    //         // Put pivot column in the adequate position and update P
    //         h.swap_cols(j, col_pivot);
    //         p.swap(j, col_pivot);

    //         // Put pivot row in the adequate position and update U
    //         h.swap_rows(j + m - n, row_pivot);
    //         u.swap_rows(j + m - n, row_pivot);

    //         // Pivot is now at (j+m-n, j)

    //         // Nullify the rest of the column and update matrix U accordingly
    //         for i in 0..m {
    //             if h[(i, j)] != f.zero() && i != j + m - n {
    //                 h.add_rows(i, j + m - n);
    //                 u.add_rows(i, j + m - n);
    //             }
    //         }
    //     }
    //     Some((u, h, p))
    // }

    // pub fn max_set_of_independant_rows_f2(&mut self) -> Vec<usize> {
    //     let f = self.field;
    //     let n = self.rows;
    //     let m = self.cols;
    //     let mut row_pivot = 0;
    //     let mut col_pivot = 0;
    //     let mut p = Perm::identity(n);
    //     let mut rank = 0;

    //     while row_pivot < self.rows && col_pivot < self.cols {
    //         // Find pivot
    //         let mut i = row_pivot;
    //         while i < n && self[(i, col_pivot)] == f.zero() {
    //             i += 1;
    //         }
    //         if i == n {
    //             col_pivot += 1;
    //             continue;
    //         }
    //         self.swap_rows(i, row_pivot);
    //         p.swap(i, row_pivot);
    //         rank += 1;

    //         // Adjust all rows below pivot's row
    //         for k in row_pivot + 1..n {
    //             if self[(k, col_pivot)] == f.zero() {
    //                 continue;
    //             }
    //             for l in col_pivot + 1..m {
    //                 self[(k, l)] = f.add(self[(k, l)], self[(row_pivot, l)]);
    //             }
    //             self[(k, col_pivot)] = f.zero();
    //         }

    //         row_pivot += 1;
    //         col_pivot += 1;
    //     }
    //     p.data()[0..rank].to_vec()
    // }

    // pub fn remove_redundant_rows_f2(&mut self) {
    //     let mut tmp = self.clone();
    //     let max_set_of_indep_rows = tmp.max_set_of_independant_rows_f2();
    //     self.keep_rows(&max_set_of_indep_rows);
    // }

    /// Encodes the matrix in bytes
    ///
    /// We start by encoding numbers of rows and columns on four bytes each.
    /// The matrix data follows.
    pub fn to_bytes(&self) -> Vec<u8> {
        let len = 4 + 4 + div_ceil(self.rows * self.cols, 8);
        let mut vec = Vec::with_capacity(len);
        vec.extend_from_slice(&(self.rows as u32).to_be_bytes());
        vec.extend_from_slice(&(self.cols as u32).to_be_bytes());
        let mut byte = 0;
        let mut shift = 7;
        for i in 0..self.rows {
            for j in 0..self.cols {
                byte |= (self[(i, j)] as u8) << shift;
                if shift == 0 {
                    vec.push(byte);
                    byte = 0;
                    shift = 7;
                } else {
                    shift -= 1;
                }
            }
        }
        if (self.rows * self.cols) % 8 != 0 {
            vec.push(byte);
        }
        vec
    }

    /// Decodes bytes encoded with [to_bytes()](struct.Mat.html#method.to_bytes) to a matrix
    pub fn from_bytes(vec: &Vec<u8>, f2: &'a F2) -> Result<Self> {
        let rows = u32::from_be_bytes(vec[0..4].try_into()?) as usize;
        let cols = u32::from_be_bytes(vec[4..8].try_into()?) as usize;
        let mut mat = Mat::zero(f2, rows, cols);
        let mut k = 8;
        let mut shift = 7;
        for i in 0..rows {
            for j in 0..cols {
                mat[(i, j)] = ((vec[k] >> shift) & 1).into();
                if shift == 0 {
                    k += 1;
                    shift = 7;
                } else {
                    shift -= 1;
                }
            }
        }
        Ok(mat)
    }
}

// impl<'a, F: CharacteristicTwo + Eq> Mat<'a, F> {
//     pub fn from<'b>(f: &'a F, mat_f2: &Mat<'b, F2>) -> Self {
//         let mut mat_f2m = Mat::zero(f, mat_f2.rows(), mat_f2.cols());
//         for i in 0..mat_f2.rows() {
//             for j in 0..mat_f2.cols() {
//                 mat_f2m[(i, j)] = f.from(mat_f2.field(), mat_f2[(i, j)]);
//             }
//         }
//         mat_f2m
//     }
// }

impl<'a, F: Eq + F2FiniteExtension> Mat<'a, F> {
    /// Takes a t * n matrix on F<sub>2<sup>m</sup></sub>
    /// and outputs a mt * n matrix on F<sub>2</sub>
    /// by decomposing each coefficient on the canonical basis
    pub fn binary<'b>(&self, f2: &'b F2) -> Mat<'b, F2> {
        let f = self.field();
        let m = f.characteristic_exponent() as usize;
        let mut bin = Mat::zero(f2, m * self.rows, self.cols);
        for j in 0..self.cols {
            for i in 0..self.rows {
                let elt_as_u32 = f.elt_to_u32(self[(i, j)]);
                for k in 0..m {
                    bin[(m * i + k, j)] = (elt_as_u32 >> k) & 1;
                }
            }
        }
        bin
    }
}

// TODO: Should I rewrite this in each module ?
fn div_ceil(a: usize, b: usize) -> usize {
    a / b + if a % b == 0 { 0 } else { 1 }
}
