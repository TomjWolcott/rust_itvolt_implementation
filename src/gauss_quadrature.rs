use itertools::Itertools;
use rgsl::legendre::polynomials::legendre_Pl;
use rgsl::numerical_differentiation::{deriv_central};

pub fn gauss_lobatto_quad_points(num_steps: usize, from: f64, to: f64) -> Vec<f64> {
    // let legendre_poly = |x: f64| legendre_Pl(num_steps as i32 - 1, x.clamp(-1.0, 1.0));

    let legendre_poly_deriv = |x: f64| {
        (0..=(num_steps / 2 - 1)).map(|i| {
            let n = (num_steps % 2 + 2 * i) as i32;

            (2.0 * n as f64 + 1.0) * legendre_Pl(n, x.clamp(-1.0, 1.0))
        }).sum()
    };

    // let legendre_poly_deriv = |x|  {
    //     deriv_central(legendre_poly, x, 1e-8).unwrap().0
    // };

    let intervals = 100 * num_steps;

    let mut quad_points: Vec<f64> = (0..=intervals)
        .map(|i| 2.0 * (i as f64) / intervals as f64 - 1.0001)
        .tuple_windows()
        .filter(|(x1, x2)| legendre_poly_deriv(*x1) * legendre_poly_deriv(*x2) < 0.0)
        .map(|(x1, x2)| {
            let mut x = (x1 + x2) / 2.0;

            for _ in 0..100 {
                let prev_x = x;
                x = x - legendre_poly_deriv(x) / deriv_central(legendre_poly_deriv, x, 1e-8).unwrap().0;
                if (x - prev_x).abs() > 1e-11 { break; }
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
    let quad_points = gauss_lobatto_quad_points(3, -1.0, 1.0);

    println!("{:?}", quad_points);
}