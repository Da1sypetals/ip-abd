use std::io::{self, Write};

use affine::AffineBody;
use newton::NewtonSolver;
use poly::Polygon;

mod affine;
mod mass;
mod newton;
mod poly;
extern crate nalgebra_glm as glm;

const dt: f32 = 0.01;
const g: f32 = 9.8;

fn pause() {
    println!("Press Enter to continue...");
    io::stdout().flush().unwrap();
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
}

fn main() {
    let nodes = vec![
        glm::Vec2::new(1.0, 0.0),
        glm::Vec2::new(0.0, 1.0),
        glm::Vec2::new(1.0, 1.0),
    ];

    let polygon = Polygon::new(nodes, 3f32);
    dbg!(polygon.mass_matrix());

    let mut ab = AffineBody::new(polygon);

    let n_step = 100;
    let newton = NewtonSolver {
        max_iter: 100,
        tol: 0.001,
    };

    for i in 0..n_step {
        ab.prepare();
        let q0 = ab.qtilde;
        let q = newton.solve(&ab, q0);
        ab.post(&q);

        println!("step {}", i);
        dbg!(&ab.q);
        pause();
    }
}
