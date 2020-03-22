use rand::{rngs::ThreadRng, Rng};
use std::ops::{Index, IndexMut};

#[derive(Debug, Eq, PartialEq)]
pub struct Perm(Vec<usize>);

impl Index<usize> for Perm {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for Perm {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Perm {
    pub fn new(vec: Vec<usize>) -> Result<Self, &'static str> {
        let n = vec.len();
        let mut list = vec![0; n];
        for i in 0..n {
            if vec[i] >= n {
                return Err("Not a permutation");
            }
            if list[vec[i]] == 0 {
                list[vec[i]] = 1;
            } else {
                return Err("Not a permutation");
            }
        }
        Ok(Perm(vec))
    }

    pub fn random(rng: &mut ThreadRng, n: usize) -> Self {
        let mut cols = Vec::with_capacity(n);
        for i in 0..n {
            cols.push(i);
        }

        let mut vec = vec![0; n];
        for i in 0..n {
            let index = rng.gen_range(0, cols.len());
            vec[i] = cols[index];
            cols.swap_remove(index);
        }
        Perm(vec)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn data(&self) -> &Vec<usize> {
        &self.0
    }

    pub fn inverse(&self) -> Self {
        let mut inv = vec![0; self.len()];
        for i in 0..self.len() {
            inv[self[i]] = i;
        }
        Perm(inv)
    }

    pub fn is_permutation(&self) -> bool {
        let n = self.len();
        let mut list = vec![false; n];
        for i in 0..n {
            match list.get(self[i]) {
                Some(false) => {
                    list[self[i]] = true;
                }
                _ => return false,
            }
        }
        true
    }
}
