use crate::{affine::AffineBody, mass::Vec6};

pub struct NewtonSolver {
    pub max_iter: usize,
    pub tol: f32,
}

impl NewtonSolver {
    pub fn solve(&self, ab: &AffineBody, q0: Vec6) -> Vec6 {
        let mut q = q0;
        for _ in 0..self.max_iter {
            let grad = ab.grad(&q);
            let hess = ab.hess(&q);

            let dir = hess.try_inverse().expect("hess not invertible!") * (-grad);

            if dir.magnitude() < self.tol {
                break;
            }
            q += dir;
        }

        q
    }
}
