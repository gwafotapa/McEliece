pub struct Mat {
    rows: usize,
    cols: usize,
    data: Vec<Vec<u8>>,
}

impl Mat {
    fn new(n: usize, m: usize) -> Mat {
        let mut mat = Mat {
            rows: n,
            cols: m,
            data: Vec::with_capacity(n),
        };
        for _i in 0..mat.rows {
            let mut col = Vec::with_capacity(m);
            col.resize(m, 0);
            mat.data.push(col);
        }
        mat
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn set(&mut self, row: usize, col: usize, val: u8) {
        self.data[row][col] = val;
    }

    pub fn get(&self, row: usize, col: usize) -> u8 {
        self.data[row][col]
    }

    pub fn print(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{}", self.get(i, j));
            }
            println!();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let mut mat = Mat::zero(3, 4);
        mat.set(2, 2, 1);
        mat.set(1, 0, 1);
        mat.set(2, 3, 1);

        let m: Vec<Vec<u8>> = vec![vec![0, 0, 0, 0], vec![1, 0, 0, 0], vec![0, 0, 1, 1]];

        assert_eq!(3, mat.rows());
        assert_eq!(4, mat.cols());
        assert_eq!(&m, &mat.data);
    }
}
