use crate::{affine::AffineBody, mass::Vec6};

pub struct NewtonSolver {
    pub max_iter: usize,
    pub tol: f32,
    pub c: f32,    // Armijo 条件中的常数
    pub beta: f32, // 步长缩放因子
}

impl NewtonSolver {
    /// This converges!
    pub fn solve_damped(&self, ab: &AffineBody, q0: Vec6) -> Vec6 {
        let mut q = q0;
        let mut alpha = 1.0; // 初始步长

        for _ in 0..self.max_iter {
            let grad = ab.grad(&q);
            let hess = ab.hess(&q);

            let dir = hess.try_inverse().expect("hess not invertible!") * (-grad);

            if dir.magnitude() < self.tol {
                break;
            }

            // Armijo 条件
            while ab.potential(&(q + alpha * dir))
                > ab.potential(&q) + self.c * alpha * grad.dot(&dir)
            {
                alpha *= self.beta;
            }

            q += alpha * dir;
            alpha = 1.0; // 重置步长
        }

        q
    }

    /// This does not converge!
    pub fn solve(&self, ab: &AffineBody, q0: Vec6) -> Vec6 {
        let mut q = q0;
        for _ in 0..self.max_iter {
            let grad = ab.grad(&q);
            let hess = ab.hess(&q);

            let dir = hess.try_inverse().expect("hess not invertible!") * (-grad);

            dbg!(&dir);
            if dir.magnitude() < self.tol {
                break;
            }
            q += dir;
        }

        q
    }
}
