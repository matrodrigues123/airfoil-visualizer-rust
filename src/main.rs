use std::f64::consts::PI;

use airfoil_visualizer::parametrizations::parsec::{ParsecParameters, Airfoil, Parameters};

fn main() {
    let parsec_params = Parameters::Parsec(ParsecParameters::new(0.02,
         0.005, 
         0.43, 
         0.12, 
         0.23, 
         -0.018, 
         -0.8, 
         0.4, 
         -0.001, 
         0.0, 
         -5.0*PI/180.0, 
         25.0*PI/180.0));

    let airfoil = Airfoil::new(parsec_params, 35);

    
    
}
