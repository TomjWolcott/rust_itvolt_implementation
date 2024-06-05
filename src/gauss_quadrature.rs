use itertools::Itertools;
use rgsl::legendre::polynomials::legendre_Pl;
use rgsl::numerical_differentiation::{deriv_central};

#[allow(dead_code)]
pub enum QuadratureType {
    Trapezoidal,
    // Simpson,
    // Gauss,
    Lobatto,
}

#[allow(dead_code)]
pub fn get_quad_points(
    ty: QuadratureType, 
    from: f64, 
    to: f64, 
    num_points: usize
) -> Vec<f64> {
    match ty {
        QuadratureType::Trapezoidal => {
            // let delta = (to - from) / num_points as f64;

            (0..num_points).map(|i| (
                from + (i as f64 / num_points as f64) * (to - from)
                // if i == 0 || i == num_points - 1 { delta / 2.0 } else { delta }
            )).collect()
        },
        QuadratureType::Lobatto => {
            todo!()
        }
    }
}

pub fn gauss_lobatto_quad(num_steps: usize, from: f64, to: f64) -> Vec<f64> {
    let legendre_poly = |x: f64| legendre_Pl(num_steps as i32 - 1, x.clamp(-1.0, 1.0));

    let legendre_poly_deriv = |x|  {
        deriv_central(legendre_poly, x, 1e-4).unwrap().0
    };

    let intervals = 3 * num_steps;

    let mut quad_points: Vec<f64> = (0..=intervals)
        .map(|i| 2.0 * (i as f64) / intervals as f64 - 1.0)
        .tuple_windows()
        .filter(|(x1, x2)| legendre_poly_deriv(*x1) * legendre_poly_deriv(*x2) <= 0.0)
        .map(|(x1, x2)| {
            let mut x = (x1 + x2) / 2.0;
            let mut prev_x = x + 2.0;

            while (x - prev_x).abs() > 1e-10 {
                prev_x = x;
                x = x - legendre_poly_deriv(x) / deriv_central(legendre_poly_deriv, x, 1e-2).unwrap().0
            }

            x
        }).map(|x| (0.5 * x + 0.5) * (to - from) + from)
        .collect();

    quad_points.insert(0, from);
    quad_points.push(to);

    quad_points
}

#[test]
fn lobatto_points() {
    let quad_points = gauss_lobatto_quad(5, -1.0, 1.0);

    println!("{:?}", quad_points);
}