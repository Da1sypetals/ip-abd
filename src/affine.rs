use crate::{
    contact::ContactIp,
    dt, g,
    grad::_grad,
    grad_contact::PointContactArg,
    hess::_hess,
    mass::{Mat6x6, Vec6},
    newton::collect_contacts,
    poly::Polygon,
};

pub trait AffineDof {
    fn t(&self) -> (f32, f32);
    fn a(&self) -> (f32, f32, f32, f32);
    fn all_eig_of_a_is_positive(&self) -> bool;
}

impl AffineDof for Vec6 {
    fn t(&self) -> (f32, f32) {
        (self.x, self.y)
    }

    fn a(&self) -> (f32, f32, f32, f32) {
        (self.z, self.w, self.a, self.b)
    }

    fn all_eig_of_a_is_positive(&self) -> bool {
        // sequential major of a
        let (a11, a12, a21, a22) = self.a();

        if a11 <= 0f32 {
            // 1st sequential major
            return false;
        }
        if a11 * a22 - a21 * a12 <= 0f32 {
            // 2nd sequential major
            return false;
        }
        true
    }
}

pub struct AffineBody {
    // @ init
    pub poly: Polygon,
    pub q: Vec6,
    pub dq: Vec6,
    pub qtilde: Vec6,
    pub orth: OrthIp,
    pub contact_ip: ContactIp,
}

impl AffineBody {
    pub fn new(poly: Polygon, kappa: f32) -> Self {
        let scale = 1f32;
        let q0 = Vec6::new(0.0, 0.0, scale, 0.0, 0.0, scale);
        Self {
            poly,
            q: q0,
            dq: Vec6::zeros(),
            qtilde: Vec6::zeros(),
            orth: OrthIp { kappa },
            contact_ip: ContactIp,
        }
    }

    pub fn prepare(&mut self) {
        let grav = glm::vec2(0f32, g);
        let g_q = self.poly.uniform_accel(grav);
        self.qtilde = self.q + dt * self.dq + dt * dt * self.poly.mass_inv() * g_q;
    }

    pub fn potential(&self, q: &Vec6) -> f32 {
        // inertial
        let con = &collect_contacts(self, q);
        let qdiff = q - self.qtilde;
        let inertial = 0.5f32 * (qdiff.transpose() * self.poly.mass_matrix() * qdiff).x;
        inertial + self.orth.potential(q) + self.contact_ip.potential(q, con)
    }

    pub fn grad(&self, q: &Vec6) -> Vec6 {
        // mass matrix is symmetric
        let con = &collect_contacts(self, q);
        let inertial = self.poly.mass_matrix() * (q - self.qtilde);
        inertial + self.orth.grad(q) + self.contact_ip.grad(q, con)
    }

    pub fn hess(&self, q: &Vec6) -> Mat6x6 {
        let con = &collect_contacts(self, q);
        let inertial = self.poly.mass_matrix();
        inertial + self.orth.hess(q) + self.contact_ip.hess(q, con)
    }

    pub fn post(&mut self, new_q: &Vec6) {
        self.dq = (new_q - self.q) / dt;
        self.q = new_q.clone();
    }

    pub fn pos(&self, q: &Vec6, i: usize) -> glm::Vec2 {
        let (a11, a12, a21, a22) = q.a();
        let (tx, ty) = q.t();
        let p = self.poly.pos_init(i);
        glm::vec2(a11 * p.x + a12 * p.y + tx, a21 * p.x + a22 * p.y + ty)
    }
    pub fn posdir(&self, q: &Vec6, qdir: &Vec6, i: usize) -> glm::Vec2 {
        let xcur = self.pos(q, i);
        let xnext = self.pos(&(q + qdir), i);
        xnext - xcur
    }
}

pub struct OrthIp {
    kappa: f32,
}

impl OrthIp {
    pub fn potential(&self, q: &Vec6) -> f32 {
        let (a11, a12, a21, a22) = q.a();
        self.kappa
            * (2f32 * (a11 * a12 + a21 * a22).powi(2)
                + (a11.powi(2) + a21.powi(2) - 1f32).powi(2)
                + (a12.powi(2) + a22.powi(2) - 1f32).powi(2))
    }

    pub fn grad(&self, q: &Vec6) -> Vec6 {
        let (a11, a12, a21, a22) = q.a();
        self.kappa * _grad(a11, a12, a21, a22)
    }

    pub fn hess(&self, q: &Vec6) -> Mat6x6 {
        let (a11, a12, a21, a22) = q.a();
        self.kappa * _hess(a11, a12, a21, a22)
    }
}
