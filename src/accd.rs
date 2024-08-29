use crate::{affine::AffineBody, bound::Boundary, mass::Vec6};

#[derive(Clone, Debug)]
pub struct CcdPair {
    pub p: glm::Vec2,
    pub e: (glm::Vec2, glm::Vec2),
}

impl CcdPair {
    pub fn distance(&self) -> f32 {
        let p = self.p;
        let (a, b) = (&self.e.0, &self.e.1);

        let ab = b - a;
        let ap = p - a;
        let bp = p - b;

        // Check if the point is closest to the start of the edge
        if ab.dot(&ap) <= 0.0 {
            return ap.norm();
        }

        // Check if the point is closest to the end of the edge
        if ab.dot(&bp) >= 0.0 {
            return bp.norm();
        }

        // The point is closest to the edge itself
        let ap_proj_on_ab = ap.dot(&ab) / ab.norm_squared();
        let proj_point = a + ap_proj_on_ab * ab;
        (p - proj_point).norm()
    }
}

#[derive(Debug)]
pub struct CcdDir {
    pub p: glm::Vec2,
    pub e: (glm::Vec2, glm::Vec2),
}

pub struct AccdMassive {
    pub s: f32,
    pub t_c: f32,
    pub max_iter: u32,
}

impl AccdMassive {
    pub fn new(s: f32, max_iter: u32) -> Self {
        Self {
            s,
            t_c: 1f32,
            max_iter,
        }
    }

    pub fn toi(&self, ab: &AffineBody, q: &Vec6, qdir: &Vec6) -> f32 {
        let mut t = 1f32;
        let accd = Accd {
            s: self.s,
            t_c: self.t_c,
            max_iter: self.max_iter,
        };

        for i in 0..ab.poly.n {
            let x = ab.pos(q, i);
            let xdir = ab.posdir(q, qdir, i);
            for (e1, e2) in Boundary::edges() {
                let cpos = CcdPair { p: x, e: (e1, e2) };
                let cdir = CcdDir {
                    p: xdir,
                    e: (glm::Vec2::zeros(), glm::Vec2::zeros()),
                };
                t = t.min(accd.toi(&cpos, &cdir));
            }
        }
        t
    }
}

pub struct Accd {
    pub s: f32,
    pub t_c: f32,
    max_iter: u32,
}

impl Accd {
    pub fn new(s: f32, max_iter: u32) -> Self {
        Self {
            s,
            t_c: 1f32,
            max_iter,
        }
    }

    /// See: https://ipc-sim.github.io/C-IPC/
    /// - currently `xi` is not supported
    pub fn toi(
        &self,
        cpos: &CcdPair, // x
        cdir: &CcdDir,  // p
    ) -> f32 {
        let pbar = (cdir.p + cdir.e.0 + cdir.e.1) / 3f32;

        let mut x = cpos.clone();
        let p = CcdDir {
            p: cdir.p - pbar,
            e: (cdir.e.0 - pbar, cdir.e.1 - pbar),
        };

        let l_p = p.p.magnitude() + p.e.0.magnitude().max(p.e.1.magnitude());
        if l_p == 0f32 {
            return 1f32;
        }

        // let dsqr = pair.distance().powi(2);
        let d = x.distance();
        let g = self.s * d;

        let mut t = 0f32;
        let mut t_l = (1f32 - self.s) * d / l_p;

        for _ in 0..self.max_iter {
            x.p += t_l * p.p;
            x.e.0 += t_l * p.e.0;
            x.e.1 += t_l * p.e.1;

            let d = x.distance();
            if t > 0f32 && d < g {
                return t;
            }

            t += t_l;
            if t > self.t_c {
                return self.t_c;
            }
            t_l = 0.9f32 * d / l_p;
        }

        t
    }
}
