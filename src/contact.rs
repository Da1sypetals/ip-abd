use crate::{
    dhat,
    grad_contact::PointContactArg,
    mass::{Mat6x6, Vec6},
};

pub struct ContactIp;

impl ContactIp {
    pub fn potential(&self, q: &Vec6, con: &Vec<PointContactArg>) -> f32 {
        let mut res = 0f32;
        for c in con {
            let d = c.distance();
            res += -(d - dhat) * (d - dhat) * ((d / dhat).ln());
        }
        res
    }

    pub fn grad(&self, q: &Vec6, con: &Vec<PointContactArg>) -> Vec6 {
        let mut res = Vec6::zeros();
        for c in con {
            let d = c.distance();
            let diff1 = (d - dhat) * (-2f32 * d * (d / dhat).ln() - d + dhat) / d;
            res += c.grad() * diff1;
        }
        res
    }

    pub fn hess(&self, q: &Vec6, con: &Vec<PointContactArg>) -> Mat6x6 {
        let mut res = Mat6x6::zeros();

        for c in con {
            let d = c.distance();
            let diff1 = (d - dhat) * (-2f32 * d * (d / dhat).ln() - d + dhat) / d;
            let diff2 = -2f32 * (d / dhat).ln() - 3f32 + 2f32 * dhat / d + (dhat * dhat) / (d * d);
            let d_grad = c.grad();
            let d_hess = c.hess();

            res += diff2 * &d_grad * d_grad.transpose() + diff1 * d_hess;
        }

        res
    }
}
