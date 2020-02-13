extern crate rand;

use rand::Rng;
use std::char;
use std::fmt;

#[derive(Eq, PartialEq)]
pub struct Mat {
    rows: usize,
    cols: usize,
    data: Vec<u8>,
}

impl Clone for Mat {
    fn clone(&self) -> Mat {
        let mut cln = Mat::new(self.rows, self.cols);
        cln.data = self.data.clone();

        cln
    }
}

impl fmt::Debug for Mat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\n{}", self.to_str())
    }
}

impl Mat {
    pub fn new(n: usize, m: usize) -> Mat {
        let mut mat = Mat {
            rows: n,
            cols: m,
            data: Vec::with_capacity(n * m),
        };
        mat.data.resize(n * m, 0);

        mat
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn set(&mut self, row: usize, col: usize, val: u8) {
        self.data[row * self.cols + col] = val;
    }

    pub fn get(&self, row: usize, col: usize) -> u8 {
        self.data[row * self.cols + col]
    }

    pub fn print(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                print!("{} ", self.get(i, j));
            }
            println!("{}", self.get(i, self.cols - 1));
        }
    }

    pub fn to_str(&self) -> String {
        let mut s = String::new();

        for i in 0..self.rows {
            for j in 0..self.cols - 1 {
                match char::from_digit(self.get(i, j) as u32, 10) {
                    Some(chr) => s.push(chr),
                    None => panic!("char::from_digit returned None"),
                }
                s.push(' ');
            }
            match char::from_digit(self.get(i, self.cols - 1) as u32, 10) {
                Some(chr) => s.push(chr),
                None => panic!("char::from_digit returned None"),
            }
            s.push('\n');
        }

        s
    }

    pub fn random(rng: &mut rand::rngs::ThreadRng, n: usize, m: usize) -> Mat {
        let mut mat = Mat::new(n, m);
        for i in 0..n {
            for j in 0..m {
                mat.set(i, j, rng.gen_range(0, 2));
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
            let mut one = false; // true when '1' is found
            let mut j = 0; // column iterator
            while !one && j != n {
                if self.get(i, j) == 1 {
                    if cols[j] == 1 {
                        return false;
                    }
                    cols[j] = 1;
                    one = true;
                }
                j += 1;
            }
            if !one {
                // row contains zeroes only
                return false;
            }
        }

        true
    }

    pub fn permutation_random(rng: &mut rand::rngs::ThreadRng, n: usize) -> Mat {
        let mut mat = Mat::new(n, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..n {
            // Draw a random column index j to set the '1' on this row
            let nbr = rng.gen_range(0, n - i);
            mat.set(i, cols[nbr], 1);

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }

        mat
    }

    pub fn identity(n: usize) -> Mat {
        let mut id = Mat::new(n, n);
        for i in 0..n {
            id.set(i, i, 1);
        }

        id
    }

    pub fn is_invertible(&self) -> bool {
        if self.rows != self.cols {
            return false;
        }

        let n = self.rows;
        let mut mat = self.clone();
        for j in 0..n {
            // Find a pivot in column j
            let mut row_pivot = n;
            for i in j..n {
                if mat.get(i, j) == 1 {
                    row_pivot = i;
                }
            }
            if row_pivot == n {
                return false;
            }

            // Swap row j with the pivot's row
            if row_pivot != j {
                for k in 0..n {
                    let tmp = mat.get(j, k);
                    mat.set(j, k, mat.get(row_pivot, k));
                    mat.set(row_pivot, k, tmp);
                }
            }

            // Add line j to remaining lines as needed
            for i in 0..n {
                if mat.get(i, j) == 1 && i != j {
                    for k in 0..n {
                        mat.set(i, k, mat.get(i, k) ^ mat.get(j, k));
                    }
                }
            }
        }

        true
    }

    // Inverse computation via gaussian elimination
    // See https://en.wikipedia.org/wiki/Gaussian_elimination
    pub fn inverse(&self) -> Option<Mat> {
        if self.rows != self.cols {
            return None;
        }

        let n = self.rows;
        let mut mat = self.clone();
        let mut inv = Mat::identity(n);
        for j in 0..n {
            // Find a pivot in column j
            let mut row_pivot = n;
            for i in j..n {
                if mat.get(i, j) == 1 {
                    row_pivot = i;
                }
            }
            if row_pivot == n {
                return None;
            }

            // Swap row j with the pivot's row
            if row_pivot != j {
                for k in 0..n {
                    let tmp = mat.get(j, k);
                    mat.set(j, k, mat.get(row_pivot, k));
                    mat.set(row_pivot, k, tmp);

                    // Mimic operation on identity matrix
                    let tmp = inv.get(j, k);
                    inv.set(j, k, inv.get(row_pivot, k));
                    inv.set(row_pivot, k, tmp);
                }
            }

            // Add line j to remaining lines as needed
            for i in 0..n {
                if mat.get(i, j) == 1 && i != j {
                    for k in 0..n {
                        mat.set(i, k, mat.get(i, k) ^ mat.get(j, k));

                        // Mimic operation on identity matrix
                        inv.set(i, k, inv.get(i, k) ^ inv.get(j, k));
                    }
                }
            }
        }

        Some(inv)
    }

    pub fn invertible_random(rng: &mut rand::rngs::ThreadRng, n: usize) -> Mat {
        let mut mat = Mat::new(n, n);
        let mut cln = Mat::new(n, n);
        let mut i = 0;
        while i < n {
            // Fill line i at random
            for j in 0..n {
                mat.set(i, j, rng.gen_range(0, 2));
            }

            // Copy matrix (remember lines after ith are null)
            for k in 0..i + 1 {
                for j in 0..n {
                    cln.set(k, j, mat.get(k, j));
                }
            }

            // Check line i is independant from lines above
            if cln.gauss() == i + 1 {
                i += 1;
            }
        }

        mat
    }

    // Applies gaussian elimination and returns matrix rank
    pub fn gauss(&mut self) -> usize {
        let mut rank = 0;
        let mut j = 0;
        loop {
            // Find a pivot in column j
            let mut row_pivot = self.rows;
            for i in j..self.rows {
                if self.get(i, j) == 1 {
                    row_pivot = i;
                }
            }
            if row_pivot == self.rows {
                return rank;
            }

            // Swap row j with the pivot's row
            if row_pivot != j {
                for k in 0..self.cols {
                    let tmp = self.get(j, k);
                    self.set(j, k, self.get(row_pivot, k));
                    self.set(row_pivot, k, tmp);
                }
            }

            // Add line j to remaining lines as needed
            for i in 0..self.rows {
                if self.get(i, j) == 1 && i != j {
                    for k in 0..self.cols {
                        self.set(i, k, self.get(i, k) ^ self.get(j, k));
                    }
                }
            }

            j += 1;
            rank += 1;
            if rank == self.rows || j == self.cols {
                return rank;
            }
        }
    }

    pub fn add(&mut self, mat1: &Mat, mat2: &Mat) {
        if self.rows != mat1.rows
            || self.rows != mat2.rows
            || self.cols != mat1.cols
            || self.cols != mat2.cols
        {
            panic!("Cannot add matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                self.set(i, j, mat1.get(i, j) ^ mat2.get(i, j));
            }
        }
    }

    pub fn mul(&mut self, mat1: &Mat, mat2: &Mat) {
        if self.rows != mat1.rows || self.cols != mat2.cols || mat1.cols != mat2.rows {
            panic!("Cannot multiply matrices: dimensions don't match");
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = 0;
                for k in 0..mat1.cols {
                    sum ^= mat1.get(i, k) & mat2.get(k, j);
                }
                self.set(i, j, sum);
            }
        }
    }

    // Generates a random row vector of length n and weight t
    pub fn weighted_vector_random(rng: &mut rand::rngs::ThreadRng, n: usize, t: usize) -> Mat {
        let mut vec = Mat::new(1, n);
        let mut cols = Vec::with_capacity(n); // remaining column indices
        for i in 0..n {
            cols.push(i);
        }

        for i in 0..t {
            // Draw a random column index
            let nbr = rng.gen_range(0, n - i);
            vec.set(0, cols[nbr], 1);

            // Remove the index from the list by putting it at the end
            cols.swap(nbr, n - 1 - i);
        }

        vec
    }

    pub fn weight(&self) -> Option<u8> {
        if self.rows != 1 {
            return None;
        }

        let mut cnt = 0;
        for j in 0..self.cols {
            cnt += self.get(0, j);
        }

        Some(cnt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let mut mat = Mat::new(3, 4);
        mat.set(2, 2, 1);
        mat.set(1, 0, 1);
        mat.set(2, 3, 1);

        let v: Vec<u8> = vec![0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1];

        assert_eq!(3, mat.rows());
        assert_eq!(4, mat.cols());
        assert_eq!(&mat.data, &v);
    }

    #[test]
    fn test_is_permutation() {
        let mut mat = Mat::new(5, 6);
        mat.set(0, 4, 1);
        mat.set(1, 2, 1);
        mat.set(2, 3, 1);
        mat.set(3, 5, 1);
        mat.set(4, 0, 1);
        println!();
        mat.print();
        assert!(!mat.is_permutation());

        let mut mat = Mat::new(6, 6);
        mat.set(0, 4, 1);
        mat.set(1, 2, 1);
        mat.set(3, 3, 1);
        mat.set(3, 5, 1);
        mat.set(4, 0, 1);
        mat.set(5, 1, 1);
        println!();
        mat.print();
        assert!(!mat.is_permutation());

        let mut mat = Mat::new(6, 6);
        mat.set(0, 4, 1);
        mat.set(1, 3, 1);
        mat.set(2, 3, 1);
        mat.set(3, 5, 1);
        mat.set(4, 0, 1);
        mat.set(5, 1, 1);
        println!();
        mat.print();
        assert!(!mat.is_permutation());

        let mut mat = Mat::new(6, 6);
        mat.set(0, 4, 1);
        mat.set(1, 2, 1);
        mat.set(2, 3, 1);
        mat.set(3, 5, 1);
        mat.set(4, 0, 1);
        mat.set(5, 1, 1);
        println!();
        mat.print();
        assert!(mat.is_permutation());
    }

    #[test]
    fn test_permutation_random() {
        let mut rng = rand::thread_rng();
        let mat = Mat::permutation_random(&mut rng, 10);
        println!();
        mat.print();
        assert!(mat.is_permutation());
    }

    #[test]
    fn test_is_invertible() {
        let id = Mat::identity(11);
        println!();
        id.print();
        assert!(id.is_invertible());

        let mut rng = rand::thread_rng();
        let mat = Mat::permutation_random(&mut rng, 20);
        println!();
        mat.print();
        assert!(mat.is_invertible());
    }

    #[test]
    fn test_inverse() {
        let id = Mat::identity(11);
        assert_eq!(id.inverse().as_ref(), Some(&id));

        let mut rng = rand::thread_rng();
        let mat = Mat::permutation_random(&mut rng, 11);
        assert_eq!(
            mat.inverse()
                .expect("Singular permutation matrix")
                .inverse()
                .expect("Singular inverse matrix"),
            mat
        );
    }

    #[test]
    fn test_invertible_random() {
        let mut rng = rand::thread_rng();
        let mat = Mat::invertible_random(&mut rng, 15);
        assert!(mat.is_invertible());
        assert_eq!(
            mat.inverse()
                .expect("Singular permutation matrix")
                .inverse()
                .expect("Singular inverse matrix"),
            mat
        );
        let mut prod = Mat::new(15, 15);
        prod.mul(&mat, &mat.inverse().expect("Singular inverse matrix"));
        let id = Mat::identity(15);
        assert_eq!(prod, id);
    }

    // Check 'a + a = 0' and 'a + 0 = 0 + a = a'
    #[test]
    fn test_add_neutral_element() {
        let mut rng = rand::thread_rng();
        let mat = Mat::random(&mut rng, 11, 11);
        let mut sum = Mat::new(11, 11);
        sum.add(&mat, &mat);
        let zero = Mat::new(11, 11);
        assert_eq!(sum, zero);

        sum.add(&mat, &zero);
        assert_eq!(sum, mat);

        sum.add(&zero, &mat);
        assert_eq!(sum, mat);
    }

    #[test]
    #[should_panic]
    fn test_mul_wrong_dimensions() {
        let mut prod = Mat::new(5, 5);
        let mat1 = Mat::new(5, 4);
        let mat2 = Mat::new(3, 5);
        prod.mul(&mat1, &mat2);
    }

    // Check (ab)c = a(bc)
    #[test]
    fn test_mul_associativity() {
        let mut rng = rand::thread_rng();

        let a = Mat::random(&mut rng, 10, 8);
        let b = Mat::random(&mut rng, 8, 13);
        let c = Mat::random(&mut rng, 13, 4);

        let mut ab = Mat::new(10, 13);
        ab.mul(&a, &b);
        let mut bc = Mat::new(8, 4);
        bc.mul(&b, &c);

        let mut abc1 = Mat::new(10, 4);
        abc1.mul(&ab, &c);
        let mut abc2 = Mat::new(10, 4);
        abc2.mul(&a, &bc);

        assert_eq!(abc1, abc2);
    }

    #[test]
    fn test_weighted_vector_random() {
        let mat = Mat::new(3, 4);
        assert!(mat.weight() == None);

        let mut rng = rand::thread_rng();
        let vec = Mat::weighted_vector_random(&mut rng, 35, 13);
        assert!(vec.weight().expect("Cannot compute vector's weight") == 13);
    }
}
