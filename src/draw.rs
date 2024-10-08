use crate::affine::AffineBody;
use macroquad::prelude::*;

pub fn draw_point(p: &glm::Vec2) {
    let _p = (p + glm::vec2(1f32, 1f32)) / 2f32;
    let p_unnorm = glm::vec2(_p[0] * screen_width(), _p[1] * screen_height());
    // dbg!(&p_unnorm);
    draw_circle(p_unnorm[0], p_unnorm[1], 5.0, YELLOW);
}

pub fn draw_link(p1: &glm::Vec2, p2: &glm::Vec2) {
    let _p1 = (p1 + glm::vec2(1f32, 1f32)) / 2f32;
    let _p2 = (p2 + glm::vec2(1f32, 1f32)) / 2f32;
    let p1_unnorm = glm::vec2(_p1[0] * screen_width(), _p1[1] * screen_height());
    let p2_unnorm = glm::vec2(_p2[0] * screen_width(), _p2[1] * screen_height());
    draw_line(
        p1_unnorm.x,
        p1_unnorm.y,
        p2_unnorm.x,
        p2_unnorm.y,
        1.,
        WHITE,
    );
}

impl AffineBody {
    pub fn draw(&self) {
        for i in 0..self.poly.n {
            let x = self.pos(&self.q, i);
            draw_point(&x);
        }
        for i in 0..self.poly.n - 1 {
            let x1 = self.pos(&self.q, i);
            let x2 = self.pos(&self.q, i + 1);
            draw_link(&x1, &x2);
        }
        let x1 = self.pos(&self.q, self.poly.n - 1);
        let x2 = self.pos(&self.q, 0);
        draw_link(&x1, &x2);
    }
}
