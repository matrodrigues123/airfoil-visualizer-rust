use std::f64::consts::PI;

use nalgebra::{Matrix6, Vector6};

pub struct ParsecParameters {
    r_le_up: f64,   // radius of the leading edge in the upper surface
    r_le_lo: f64,   // radius of the leading edge in the lower surface
    x_up: f64,     // x coordinate of the maximum position of the upper surface
    z_up: f64,     // z coordinate of the maximum position of the upper surface
    x_lo: f64,     // x coordinate of the minimum position of the lower surface
    z_lo: f64,     // z coordinate of the minimum position of the lower surface
    z_xxup: f64,   // curvature in the position of maximum thickness of the upper surface
    z_xxlo: f64,   // curvature in the position of maximum thickness of the lower surface
    z_te: f64,     // z coordinate of the trailing edge
    delta_z_te: f64,  // separation between the upper and lower surface in the trailing edge 
    alpha_te: f64,   // angle of the inclination of the mean chamber line in the tailing edge
    beta_te: f64,    // angle of separation of the upper and lower surface in the tailing edge
}

impl ParsecParameters {
    pub fn new(
        r_le_up: f64,   
        r_le_lo: f64,   
        x_up: f64,     
        z_up: f64,     
        x_lo: f64,     
        z_lo: f64,     
        z_xxup: f64,   
        z_xxlo: f64,   
        z_te: f64,     
        delta_z_te: f64,  
        alpha_te: f64,   
        beta_te: f64,   
    ) -> Self {
        ParsecParameters {
            r_le_up,   
            r_le_lo,   
            x_up,     
            z_up,     
            x_lo,     
            z_lo,     
            z_xxup,   
            z_xxlo,   
            z_te,     
            delta_z_te,  
            alpha_te,   
            beta_te,   
        }
    }
}

pub enum Parameters {
    Parsec(ParsecParameters),
    BSpline,
}

struct Coordinate {
    x: f64,
    y: f64,
}

pub struct Airfoil {
    coordinates: Vec<Coordinate>,
}

impl Airfoil {
    pub fn new(params: Parameters, num_points: usize) -> Self {
        match params {
            Parameters::Parsec(p) => {
                Airfoil {
                    coordinates: generate_coordinates(&p, num_points)
                }
            },
            Parameters::BSpline => todo!(),
        }
    }
}

fn coef_up(p: &ParsecParameters) -> Vector6<f64> {

    let a = Matrix6::new(
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        p.x_up.powi(1/2), p.x_up.powi(3/2), p.x_up.powi(5/2), p.x_up.powi(7/2), p.x_up.powi(9/2), p.x_up.powi(11/2),
        0.5, 1.5, 2.5, 3.5, 4.5, 5.5,
        0.5*p.x_up.powi(-1/2), 1.5*p.x_up.powi(1/2), 2.5*p.x_up.powi(3/2), 3.5*p.x_up.powi(5/2), 4.5*p.x_up.powi(7/2), 5.5*p.x_up.powi(9/2),
        -0.25*p.x_up.powi(-3/2), 0.75*p.x_up.powi(-1/2), 3.75*p.x_up.powi(1/2), 8.75*p.x_up.powi(3/2), 15.75*p.x_up.powi(5/2), 24.75*p.x_up.powi(7/2),
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0
    );

    let b = Vector6::new(
        p.z_te + 0.5*p.delta_z_te, 
        p.z_up, 
        (p.alpha_te - 0.5*p.beta_te).tan(), 
        0.0,
        p.z_xxup, 
        p.r_le_up.sqrt());

    // Solve the linear system: a * coef_up = b
    let decomp = a.lu();
    decomp.solve(&b).expect("Linear solution failed")
}

fn coef_low(p: &ParsecParameters) -> Vector6<f64> {

    let a = Matrix6::new(
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        p.x_up.powi(1/2), p.x_up.powi(3/2), p.x_up.powi(5/2), p.x_up.powi(7/2), p.x_up.powi(9/2), p.x_up.powi(11/2),
        0.5, 1.5, 2.5, 3.5, 4.5, 5.5,
        0.5*p.x_up.powi(-1/2), 1.5*p.x_up.powi(1/2), 2.5*p.x_up.powi(3/2), 3.5*p.x_up.powi(5/2), 4.5*p.x_up.powi(7/2), 5.5*p.x_up.powi(9/2),
        -0.25*p.x_up.powi(-3/2), 0.75*p.x_up.powi(-1/2), 3.75*p.x_up.powi(1/2), 8.75*p.x_up.powi(3/2), 15.75*p.x_up.powi(5/2), 24.75*p.x_up.powi(7/2),
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0
    );

    let b = Vector6::new(
        p.z_te - 0.5*p.delta_z_te, 
        p.z_lo, 
        (p.alpha_te - 0.5*p.beta_te).tan(), 
        0.0,
        p.z_xxlo, 
        -p.r_le_lo.sqrt());

    // Solve the linear system: a * coef_low = b
    let decomp = a.lu();
    decomp.solve(&b).expect("Linear solution failed")
}

fn generate_coordinates(p: &ParsecParameters, num_points: usize) -> Vec<Coordinate> {
    // cosine distribution to generate x coordinates
    let mut x: Vec<f64> = Vec::new();
    for i in 1..=num_points {
        x.push(0.5 * (1.0 - f64::cos(PI*(i as f64 - 1.0)/(num_points as f64 - 1.0))));
    }

    let mut coords: Vec<Coordinate> = Vec::new();

    let a_coef: Vector6<f64> = coef_up(p);
    for x_coord in x.iter() {
        let mut z_coord_up = 0.0;
        for n in 1..=6 {
            z_coord_up += a_coef[n] * x_coord.powf(n as f64 - 0.5);
        }
        coords.push(
            Coordinate { x: *x_coord, y: z_coord_up }
        );
    }

    let b_coef: Vector6<f64> = coef_low(p);
    for x_coord in x.iter().rev() {
        let mut z_coord_low = 0.0;
        for n in 1..=6 {
            z_coord_low += b_coef[n] * x_coord.powf(n as f64 - 0.5);
        }
        coords.push(
            Coordinate { x: *x_coord, y: z_coord_low }
        );
    }

    coords
}
