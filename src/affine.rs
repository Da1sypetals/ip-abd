use crate::{
    dt, g,
    mass::{Mat6x6, Vec6},
    poly::Polygon,
};

pub struct AffineBody {
    // @ init
    pub poly: Polygon,
    pub q: Vec6,
    pub dq: Vec6,
    pub qtilde: Vec6,
}

impl AffineBody {
    pub fn new(poly: Polygon) -> Self {
        let q0 = Vec6::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0);
        Self {
            poly,
            q: q0,
            dq: Vec6::zeros(),
            qtilde: Vec6::zeros(),
        }
    }

    pub fn prepare(&mut self) {
        let grav = glm::vec2(0f32, g);
        let g_q = self.poly.uniform_accel(grav);
        self.qtilde = self.q + dt * self.dq + dt * dt * self.poly.mass_inv() * g_q;
    }

    pub fn potential(&self, q: &Vec6) -> f32 {
        // inertial
        let qdiff = q - self.qtilde;
        0.5f32 * (qdiff.transpose() * self.poly.mass_matrix() * qdiff).x
    }

    pub fn grad(&self, q: &Vec6) -> Vec6 {
        // mass matrix is symmetric
        self.poly.mass_matrix() * (q - self.qtilde)
    }

    pub fn hess(&self, q: &Vec6) -> Mat6x6 {
        self.poly.mass_matrix()
    }

    pub fn post(&mut self, new_q: &Vec6) {
        self.dq = (new_q - self.q) / dt;
        self.q = new_q.clone();
    }
}
