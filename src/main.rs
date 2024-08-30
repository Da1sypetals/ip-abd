use std::io::{self, Write};

use affine::AffineBody;
use macroquad::prelude::*;
use newton::NewtonSolver;
use poly::Polygon;

mod accd;
mod affine;
mod bound;
mod contact;
mod draw;
mod generated;
mod grad;
mod grad_contact;
mod hess;
mod mass;
mod newton;
mod poly;
extern crate nalgebra_glm as glm;

fn window_conf() -> Conf {
    Conf {
        window_title: String::from("Springs"),
        window_width: 1000,
        window_height: 1000,
        ..Default::default()
    }
}

const dt: f32 = 0.01;
const g: f32 = 9.8;
const dhat: f32 = 0.001;

fn pause() {
    println!("Press Enter to continue...");
    io::stdout().flush().unwrap();
    let mut input = String::new();

    io::stdin().read_line(&mut input).unwrap();
}

#[macroquad::main(window_conf)]
async fn main() {
    let polygon = Polygon::from_file("a.poly").expect("Failure when reading polygon from file!");
    // dbg!(polygon.mass_matrix());

    let kappa: f32 = 1.0;

    let mut ab = AffineBody::new(polygon, kappa);

    let newton = NewtonSolver {
        max_iter: 500,
        tol: 0.000001,
        c: 0.001,
        beta: 0.5,
    };

    for i in 0.. {
        clear_background(BLACK);

        ab.prepare();

        let q0 = ab.q;
        let q = newton.solve_damped(&ab, q0);
        ab.post(&q);

        println!("step {}", i);
        dbg!(&ab.q);
        dbg!(ab.orth.potential(&q));

        ab.draw();

        pause();
        next_frame().await
    }
}
