use crate::generated::grad::*;
use crate::generated::hess::case1::hess_case1;
use crate::generated::hess::case2::hess_case2;
use crate::generated::hess::case3::hess_case3;
use crate::mass::Mat6x6;
use crate::{affine::AffineDof, mass::Vec6};

pub struct PointContactArg {
    // edge ends
    pub u: glm::Vec2,
    pub v: glm::Vec2,

    // rest point position
    pub p: glm::Vec2,

    // q: abd configuration vector
    pub q: Vec6,
}

impl PointContactArg {
    pub fn x(&self) -> glm::Vec2 {
        let (a11, a12, a21, a22) = self.q.a();
        let (tx, ty) = self.q.t();
        glm::vec2(
            a11 * self.p.x + a12 * self.p.y + tx,
            a21 * self.p.x + a22 * self.p.y + ty,
        )
    }
    pub fn distance(&self) -> f32 {
        let uv = self.v - self.u;
        let ux = self.x() - self.u;

        let uv_len_sq = uv.dot(&uv);
        if uv_len_sq == 0.0 {
            // u and v are the same point, so return the distance from p to u
            return ux.norm();
        }

        let t = ux.dot(&uv) / uv_len_sq;

        if t < 0.0 {
            // p is closer to u
            return ux.norm();
        } else if t > 1.0 {
            // p is closer to v
            return (self.x() - self.v).norm();
        } else {
            // p is closest to the line segment uv
            let projection = self.u + t * uv;
            return (self.x() - projection).norm();
        }
    }

    pub fn grad(&self) -> Vec6 {
        let uv = self.v - self.u;
        let ux = self.x() - self.u;

        let uv_len_sq = uv.dot(&uv);
        let t = ux.dot(&uv) / uv_len_sq;

        let (a11, a12, a21, a22) = self.q.a();
        let (tx, ty) = self.q.t();
        let (px, py) = (self.p.x, self.p.y);
        let (ux, uy, vx, vy) = (self.u.x, self.u.y, self.v.x, self.v.y);

        if t < 0.0 {
            // p is closer to u
            Vec6::new(
                grad_case1_0(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
                grad_case1_1(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
                grad_case1_2(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
                grad_case1_3(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
                grad_case1_4(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
                grad_case1_5(a11, a12, a21, a22, px, py, tx, ty, ux, uy),
            )
        } else if t > 1.0 {
            // p is closer to v
            Vec6::new(
                grad_case2_0(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
                grad_case2_1(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
                grad_case2_2(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
                grad_case2_3(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
                grad_case2_4(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
                grad_case2_5(a11, a12, a21, a22, px, py, tx, ty, vx, vy),
            )
        } else {
            // p is closest to the line segment uv
            Vec6::new(
                grad_case3_0(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
                grad_case3_1(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
                grad_case3_2(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
                grad_case3_3(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
                grad_case3_4(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
                grad_case3_5(a11, a12, a21, a22, px, py, tx, ty, ux, uy, vx, vy),
            )
        }
    }

    pub fn hess(&self) -> Mat6x6 {
        let uv = self.v - self.u;
        let ux = self.x() - self.u;

        let uv_len_sq = uv.dot(&uv);
        let t = ux.dot(&uv) / uv_len_sq;

        let (a11, a12, a21, a22) = self.q.a();
        let (tx, ty) = self.q.t();
        let (px, py) = (self.p.x, self.p.y);
        let (ux, uy, vx, vy) = (self.u.x, self.u.y, self.v.x, self.v.y);

        if t < 0.0 {
            // p is closer to u
            hess_case1(a11, a12, a21, a22, px, py, tx, ty, ux, uy)
        } else if t > 1.0 {
            // p is closer to v
            hess_case2(a11, a12, a21, a22, px, py, tx, ty, ux, uy)
        } else {
            // p is closest to the line segment uv
            hess_case3()
        }
    }
}
