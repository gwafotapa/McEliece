use rand::{rngs::ThreadRng, Rng};
use std::{
    fmt::{Debug, Display, Formatter, Result},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::finite_field::{CharacteristicTwo, F2FiniteExtension, Field, FiniteField, F2};

pub use rowvec::RowVec;

#[derive(Eq, PartialEq)]
pub struct Mat<'a, F: Eq + Field> {
    field: &'a F,
    rows: usize,
    cols: usize,
    data: Vec<F::FElt>,
}

impl<'a, F: Eq + Field> Clone for Mat<'a, F> {
    fn clone(&self) -> Self {
        Mat {
            field: self.field, // field is shared
            rows: self.rows,
            cols: self.cols,
            data: self.data.clone(),
        }
    }
}

impl<'a, F: Eq + Field> Add for Mat<'a, F> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        &self + &other
    }
}

impl<'a, F: Eq + Field> Add<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<'a, F: Eq + Field> Add<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn add(self, other: Mat<'a, F>) -> Self::Output {
        self + &other
    }
}

impl<'a, F: Eq + Field> Add for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn add(self, other: Self) -> Self::Output {
        let f = self.field;
        if f != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        let mut sum = Mat::zero(f, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            sum.data[i] = f.add(self.data[i], other.data[i]);
        }
        sum
    }
}

impl<'a, F: Eq + Field> AddAssign<Mat<'a, F>> for Mat<'a, F> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<'a, F: Eq + Field> AddAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn add_assign(&mut self, other: &Self) {
        let f = self.field;
        if f != other.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = f.add(self.data[i], other.data[i]);
        }
    }
}

impl<'a, F: Eq + Field> Sub for Mat<'a, F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<'a, F: Eq + Field> Sub<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        &self - other
    }
}

impl<'a, F: Eq + Field> Sub<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn sub(self, other: Mat<'a, F>) -> Self::Output {
        self - &other
    }
}

impl<'a, F: Eq + Field> Sub for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn sub(self, other: Self) -> Self::Output {
        let f = self.field;
        if f != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        let mut diff = Mat::zero(f, self.rows, self.cols);
        for i in 0..self.rows * self.cols {
            diff.data[i] = f.sub(self.data[i], other.data[i]);
        }
        diff
    }
}

impl<'a, F: Eq + Field> SubAssign<Mat<'a, F>> for Mat<'a, F> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<'a, F: Eq + Field> SubAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn sub_assign(&mut self, other: &Self) {
        let f = self.field;
        if f != other.field {
            panic!("Cannot substract matrices: fields don't match");
        } else if self.rows != other.rows || self.cols != other.cols {
            panic!("Cannot substract matrices: dimensions don't match");
        }

        for i in 0..self.rows * self.cols {
            self.data[i] = f.sub(self.data[i], other.data[i]);
        }
    }
}

impl<'a, F: Eq + Field> Mul for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<'a, F: Eq + Field> Mul<&Mat<'a, F>> for Mat<'a, F> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<'a, F: Eq + Field> Mul<Mat<'a, F>> for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, other: Mat<'a, F>) -> Self::Output {
        self * &other
    }
}

impl<'a, F: Eq + Field> Mul for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn mul(self, other: Self) -> Self::Output {
        let f = self.field;
        if f != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let mut prod = Mat::zero(f, self.rows, other.cols);
        for i in 0..prod.rows {
            for j in 0..prod.cols {
                for k in 0..self.cols {
                    prod[(i, j)] = f.add(prod[(i, j)], f.mul(self[(i, k)], other[(k, j)]));
                }
            }
        }
        prod
    }
}

impl<'a, F: Eq + Field> MulAssign<Mat<'a, F>> for Mat<'a, F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<'a, F: Eq + Field> MulAssign<&Mat<'a, F>> for Mat<'a, F> {
    fn mul_assign(&mut self, other: &Self) {
        let f = self.field;
        if f != other.field {
            panic!("Cannot multiply matrices: fields don't match");
        } else if self.cols != other.rows || self.cols != other.cols {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        let tmp = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = f.zero();
                for k in 0..tmp.cols {
                    self[(i, j)] = f.add(self[(i, j)], f.mul(tmp[(i, k)], other[(k, j)]));
                }
            }
        }
    }
}

impl<'a, F: Eq + Field> Neg for Mat<'a, F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, F: Eq + Field> Neg for &Mat<'a, F> {
    type Output = Mat<'a, F>;

    fn neg(self) -> Self::Output {
        let f = self.field;
        let mut opp = Mat::zero(f, self.rows, self.cols);

        for i in 0..self.rows * self.cols {
            opp.data[i] = f.neg(self.data[i]);
        }
        opp
    }
}

impl<'a, F: Eq + Field> Index<(usize, usize)> for Mat<'a, F> {
    type Output = F::FElt;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.cols + index.1]
    }
}

impl<'a, F: Eq + Field> IndexMut<(usize, usize)> for Mat<'a, F> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.cols + index.1]
    }
}

// impl<'a> Debug for Mat<'a, F2> {
//     fn fmt(&self, f: &mut Formatter<'_>) -> Result {
//         let k = self.field;
//         write!(f, "\n")?;
//         for i in 0..self.rows {
//             for j in 0..self.cols - 1 {
//                 write!(f, "{} ", self[(i, j)].to_string())?;
//             }
//             write!(f, "{}\n", self[(i, self.cols - 1)].to_string())?;
//         }
//         Ok(())
//     }
// }

impl<'a, F: Eq + F2FiniteExtension> Debug for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let k = self.field;
        let chr = k.characteristic_exponent() as usize;
        let width = chr / 4 + if chr % 4 == 0 { 0 } else { 1 };
        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(
                    f,
                    "{:>w$x} ",
                    k.elt_to_u32(self[(i, j)]),
                    w = width
                )?;
            }
            write!(
                f,
                "{:>w$}\n",
                k.elt_to_u32(self[(i, self.cols - 1)]),
                w = width
            )?;
        }
        Ok(())
    }
}

impl<'a, F: Eq + FiniteField> Display for Mat<'a, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let k = self.field;

        // number of digits of order
        let digits = ((32 - k.order().leading_zeros()) / 3 + 1) as usize;

        write!(f, "\n")?;
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                write!(
                    f,
                    "{:>width$} ",
                    k.elt_to_str(self[(i, j)]),
                    width = if k.order() == 2 { 1 } else { 2 + digits }
                )?;
            }
            write!(
                f,
                "{:>width$}\n",
                k.elt_to_str(self[(i, self.cols - 1)]),
                width = if k.order() == 2 { 1 } else { 2 + digits }
            )?;
        }
        Ok(())
    }
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

    pub fn data(&self) -> Vec<F::FElt> {
        self.data.clone()
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

    pub fn is_permutation(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }

        let f = self.field;
        let n = self.rows;
        let mut cols = Vec::with_capacity(n); // indices of columns containing a '1'
        cols.resize(n, 0);

        for i in 0..n {
            // loop on rows
            let mut has_one = false;
            for j in 0..n {
                if self[(i, j)] == f.zero() {
                    continue;
                } else if self[(i, j)] == f.one() {
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

    pub fn permutation_random(rng: &mut ThreadRng, f: &'a F, n: usize) -> Self {
        let mut mat = Self::zero(f, n, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..n {
            // Draw a random column index j to set the '1' on this row
            let nbr = rng.gen_range(0, n - i);
            mat[(i, cols[nbr])] = f.one();

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }
        mat
    }

    pub fn identity(f: &'a F, n: usize) -> Self {
        let mut id = Self::zero(f, n, n);
        for i in 0..n {
            id[(i, i)] = f.one();
        }
        id
    }

    pub fn sum(&mut self, mat1: &Self, mat2: &Self) {
        let f = self.field;
        if f != mat1.field || f != mat2.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != mat1.rows
            || self.rows != mat2.rows
            || self.cols != mat1.cols
            || self.cols != mat2.cols
        {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] = f.add(mat1[(i, j)], mat2[(i, j)]);
            }
        }
    }

    pub fn prod(&mut self, mat1: &Self, mat2: &Self) {
        let f = self.field;
        if f != mat1.field || f != mat2.field {
            panic!("Cannot add matrices: fields don't match");
        } else if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = f.zero();
                for k in 0..mat1.cols {
                    sum = f.add(sum, f.mul(mat1[(i, k)], mat2[(k, j)]));
                }
                self[(i, j)] = sum;
            }
        }
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::zero(self.field, self.cols, self.rows);
        for i in 0..t.rows {
            for j in 0..t.cols {
                t[(i, j)] = self[(j, i)];
            }
        }
        t
    }
}

impl<'a, F: CharacteristicTwo + Eq> Mat<'a, F> {
    pub fn from<'b>(f: &'a F, mat_f2: &Mat<'b, F2>) -> Self
// where
    //     F2: Eq + F2FiniteExtension,
    {
        let mut mat_f2m = Mat::zero(f, mat_f2.rows(), mat_f2.cols());
        for i in 0..mat_f2.rows() {
            for j in 0..mat_f2.cols() {
                mat_f2m[(i, j)] = f.from(mat_f2.field(), mat_f2[(i, j)]);
            }
        }
        mat_f2m
    }
}

impl<'a, F: Eq + F2FiniteExtension> Mat<'a, F> {
    pub fn binary_form<'b>(&self, f2: &'b F2) -> Mat<'b, F2>
// where
    //     F2: Eq + F2FiniteExtension,
    {
        let f = self.field();
        let m = f.characteristic_exponent();
        let mut bin = Mat::zero(f2, m as usize * self.rows, self.cols);
        for j in 0..self.cols {
            for i in 0..self.rows {
                for k in 0..m as usize {
                    bin[(m as usize * i + k, j)] = match (f.elt_to_u32(self[(i, j)]) >> k) & 1 {
                        0 => f2.zero(),
                        1 => f2.one(),
                        _ => panic!("Unexpected value"),
                    }
                }
            }
        }
        bin
    }
}

impl<'a> Mat<'a, F2> {
    pub fn to_hex_string(&self) -> String { // TODO: Should return type be String or str or &str ??
        if self.rows > 255 || self.cols > 255 {
            panic!("Cannot convert matrix to hex string: dimensions not supported");
        }
        let len = 2 * (4 + 1) + 2 * (self.rows * self.cols / 8 + 1);
        let mut s = String::with_capacity(len);
        s.push_str(format!("{:x}#{:x}#", self.rows, self.cols).as_str());
        let mut byte = 0;
        let mut shift = 7;
        for i in 0..self.rows {
            for j in 0..self.cols {
                // let shift = (self.cols * i + j) % 8;
                byte |= self[(i, j)] << shift;
                if shift == 0 {
                    s.push_str(format!("{:02x}", byte).as_str());
                    byte = 0;
                    shift = 7;
                } else {
                    shift -= 1;
                }
            }
        }
        s
    }

    // TODO: Shouldn't I create a f2 field and return it along the matrix ?
    pub fn from_hex_string(s: &str, f2: &'a F2) -> Self {
        // TODO: how to deal with badly formed input ?
        let v: Vec<&str> = s.split('#').collect();
        let rows = usize::from_str_radix(v[0], 16).unwrap();
        let cols = usize::from_str_radix(v[1], 16).unwrap();
        let data = hex::decode(v[2]).expect("Hex decoding failed"); // TODO: use a result and ?
                                                                 // if data.len() < 3 {
                                                                 //     panic!("Cannot convert hex string to matrix: string is too short");
                                                                 // }
        // let rows = data[0];
        // let cols = data[1];
        // if data.len() != 2 + rows * cols {
        //     panic!("Cannot convert hex string to matrix: wrong matrix dimensions");
        // }
        // println!("rows: {}\ncols: {}\n", rows, cols);
        let mut mat = Mat::zero(f2, rows, cols); // TODO: into() or as ?
        let mut d = 0;
        let mut shift = 7;
        for i in 0..rows {
            for j in 0..cols {
                mat[(i, j)] = ((data[d] >> shift) & 1).into();
                if shift == 0 {
                    d += 1;
                    shift = 7;
                } else {
                    shift -= 1;
                }
            }
        }
        // for i in 0..data.len() {
        //     for j in 0..8 {
        //         let bit_index = 8 * i + j;
        //         let row = bit_index / cols;
        //         let col = bit_index % cols;
        //         mat[(row, col)] = ((data[i] >> 7 - j) & 1).into();
        //     }
        // }
        // println!("{}\n", mat);
        mat
    }
}

mod gauss;
mod rowvec;
