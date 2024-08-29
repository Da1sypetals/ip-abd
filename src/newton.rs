use crate::{
    affine::{AffineBody, AffineDof},
    hess::subhess,
    mass::Vec6,
};

pub struct NewtonSolver {
    pub max_iter: usize,
    pub tol: f32,
    pub c: f32,    // Armijo 条件中的常数
    pub beta: f32, // 步长缩放因子
}

impl NewtonSolver {
    /*
        pub fn solve_affine(&self, ab: &AffineBody, q0: Vec6) -> Vec6 {
            let mut q = q0;

            let mut cnt = 0;
            for _ in 0..self.max_iter {
                let grad = ab.orth.grad(&q);
                let grad4 = glm::vec4(grad.z, grad.w, grad.a, grad.b);
                let hess4 = subhess(q.z, q.w, q.a, q.b);
                println!("{}", hess4);

                let dir4 = hess4.try_inverse().expect("hess not invertible!") * (-grad4);
                let dir = Vec6::new(0.0, 0.0, dir4.x, dir4.y, dir4.z, dir4.w);

                if dir.magnitude() < self.tol {
                    break;
                }

                // Armijo 条件
                let mut alpha = 1f32;
                while ab.potential(&(q + alpha * dir))
                    > ab.potential(&q) + self.c * alpha * grad.dot(&dir)
                {
                    alpha *= self.beta;
                }
                println!("alpha = {}", alpha);

                q += alpha * dir;

                cnt += 1;
            }
            println!("Newton solver iterated {} times", cnt);

            q
        }
    */

    /// This converges!
    pub fn solve_damped(&self, ab: &AffineBody, q0: Vec6) -> Vec6 {
        let mut q = q0;

        let mut cnt = 0;
        for _ in 0..self.max_iter {
            let grad = ab.grad(&q);
            let hess = ab.hess(&q);

            let dir = hess.try_inverse().expect("hess not invertible!") * (-grad);
            if dir.magnitude() < self.tol {
                break;
            }

            let mut alpha = 1f32;
            // // inversion aware line search
            // while !((q + alpha * dir).all_eig_of_a_is_positive()) {
            //     alpha *= self.beta;
            // }

            // armijo condition
            while ab.potential(&(q + alpha * dir))
                > ab.potential(&q) + self.c * alpha * grad.dot(&dir)
            {
                alpha *= self.beta;
            }
            println!("alpha = {}", alpha);

            q += alpha * dir;

            cnt += 1;
        }
        println!("Newton solver iterated {} times", cnt);

        q
    }

    /// This does not converge!
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
