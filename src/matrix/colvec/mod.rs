use rand::Rng;
use std::rc::Rc;

use super::Mat;
use crate::finite_field::{Field, FieldTrait};

#[derive(Eq, PartialEq)]
pub struct ColVec<F>(pub Mat<F>)
where
    F: FieldTrait;

impl<F> ColVec<F>
where
    F: FieldTrait,
{
    pub fn new(field: Field<F>, data: Vec<F::FieldElement>) -> Self {
        ColVec(Mat::new(field, data.len(), 1, data))
    }

    pub fn zero(field: Field<F>, rows: usize) -> Self {
        ColVec(Mat::zero(field, rows, 1))
    }

    pub fn field(&self) -> &Rc<F> {
        self.0.field()
    }

    pub fn rows(&self) -> usize {
        self.0.rows()
    }

    pub fn cols(&self) -> usize {
        1
    }

    pub fn data(&self) -> &Vec<F::FieldElement> {
        self.0.data()
    }

    pub fn weight(&self) -> usize {
        let mut weight = 0;
        for i in 0..self.rows() {
            if self[i] != self.0.field.zero() {
                weight += 1;
            }
        }
        weight
    }

    pub fn random(field: Field<F>, n: usize) -> Self {
        ColVec(Mat::random(field, n, 1))
    }

    pub fn random_with_weight(field: Field<F>, n: usize, w: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut vec = ColVec::zero(field, n);
        let mut rows = Vec::with_capacity(n);
        for i in 0..n {
            rows.push(i);
        }

        for _i in 0..w {
            let mut elt;
            loop {
                elt = vec.field().random_element(&mut rng);
                if elt != vec.field().zero() {
                    break;
                }
            }
            let index = rng.gen_range(0, rows.len());
            vec[rows[index]] = elt;
            rows.swap_remove(index);
        }
        vec
    }

    pub fn transpose(&self) -> Mat<F> {
        self.0.transpose()
    }

    pub fn is_zero(&self) -> bool {
        for i in 0..self.rows() {
            if self[i] != self.field().zero() {
                return false;
            }
        }
        true
    }

    // pub fn extract_rows(&self, perm: &Vec<usize>) -> Self {
    //     ColVec(self.0.extract_rows(perm))
    // }
}

// pub mod io;
pub mod traits;
