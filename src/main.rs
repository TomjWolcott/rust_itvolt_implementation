mod complex_wrapper;
use std::ops::{Add, Mul};

use complex_wrapper::*;

mod gauss_quadrature;
use gauss_quadrature::*;

use gauss_quad::*;


fn main() {}

// fn get_quad_points(num: usize, from: f64, to: f64) -> Vec<f64> {
//     // Right now, I'm just doing evenly spaced points
//     (0..num)
//         .map(|i| (i as f64 / num as f64) * (to - from) + from)
//         .collect()
// }]

#[test]
fn test_gauss_quad() {
    let gauss = GaussLegendre::init(3);

    println!("{:?}", gauss);
}

#[test]
fn make_lagrange_weights() {
    let poly1 = lagrange(&vec![1., 2., 3.11, 4.32], 2);
    
    println!("{} -- {} -- {} -- {} -- {}", poly1(1.), poly1(2.), poly1(3.11), poly1(4.32), poly1(5.17));

    println!("{}", integral(|x| x.powi(2), -1.0, 1.0, 0.01))
}

fn lagrange(points: &Vec<f64>, i: usize) -> Box<dyn Fn(f64) -> f64> {
    let points = points.clone();
    Box::new(move |x| points.iter()
        .enumerate()
        .fold(1.0, |acc, (j, point)| 
            if i == j { acc } else { acc * (x - point) / (points[i] - point) }
        )
    )
}

fn integral<T>(f: impl Fn(f64) -> T, from: f64, to: f64, delta: f64) -> T 
    where T: Add<T, Output = T> + Mul<f64, Output = T> + Default
{
    let mut x = from;
    let mut result: T = T::default();

    while x < to {
        result = result + (f(x) + f(x+delta)) * 0.5 * delta;
        x += delta;
    }

    result
}

fn weights(t_0: f64, quad_points: &Vec<f64>) -> Vec<Vec<f64>> {
    quad_points.iter().map(|&t_i| {
        quad_points.iter().enumerate().map(|(j, _)| {
            integral(lagrange(&quad_points, j), t_0, t_i, 0.001)
        }).collect()
    }).collect()
}

fn quad_points(from: f64, to: f64, amount: usize) -> Vec<f64> {
    // Evenly spaced
    // (0..=amount).map(|i| (i as f64 / amount as f64) * (to - from) + from).collect()

    // Lobatto
    gauss_lobatto_quad(amount, from, to)

    // Legendre
    // let mut quad_points = GaussLegendre::init(amount).nodes;
    // quad_points.insert(0, 1.0);
    // quad_points.push(-1.0);
    //
    // for point in quad_points.iter_mut() {
    //     *point = (1.0 - *point) * (to - from) / 2.0 + from
    // }
    //
    // quad_points
}

#[allow(non_snake_case)]
#[test]
fn approx_simple_ode() {
    let num_quad_points = 20;
    let delta_time = 1.0;
    let total_time = 25.0;
    let max_k = 200;
    let mut max_err = 0.0;

    let mut f: Vec<Complex> = Vec::new();

    for j in 0..((total_time / delta_time) as usize) {
        let tau_start = j as f64 * delta_time;
        let tau_end = (j + 1) as f64 * delta_time;
        let t_a = (tau_start + tau_end) / 2.0;
        let last_f = f.last().cloned().unwrap_or(ONE);
        // let last_f = E.pow(-I * tau_start.powi(2) / 2.0);
        let g = |t: f64| E.pow(-I * t_a * (t - tau_start)) * last_f;
        let kernel = |t_i: f64, t_j: f64| -I * E.pow(-I * t_a * (t_i - t_j)) * (t_j - t_a);
        let quad_points = quad_points(tau_start, tau_end, num_quad_points);

        let weights = weights(tau_start, &quad_points);

        if j == 0 {
            f = quad_points.iter().map(|&t_i| g(t_i)).collect();
        }

        iterative_volterra_integration(&quad_points, &weights, &mut f, g, kernel, max_k);

        // println!("\n\nINTERVAL: {tau_start:.4} to {tau_end:.4}");

        for (&f_i, &t_i) in f.iter().zip(quad_points.iter()) {
            let correct = E.pow(-I * t_i.powi(2) / 2.0);
            // println!("f({t_i:.4}) = {f_i} vs. {correct} DIFF: {:.4e} ", (correct - f_i).len());

            max_err = (f_i - correct).len().max(max_err);
        }
    }

    println!("MAX ERROR: {:.4e}", max_err)
}

fn iterative_volterra_integration(
    quad_points: &Vec<f64>,
    weights: &Vec<Vec<f64>>,
    f: &mut Vec<Complex>, 
    g: impl Fn(f64) -> Complex, 
    kernel: impl Fn(f64, f64) -> Complex,
    steps: usize
) {
    for _k in 0..steps {
        let mut max_diff: f64 = 0.0;
        let mut new_f = Vec::new();

        for (i, &t_i) in quad_points.iter().enumerate() {
            let mut new_f_i = g(t_i);

            for (j, &t_j) in quad_points.iter().enumerate() {
                new_f_i += kernel(t_i, t_j) * weights[i][j] * f[j]
            }

            new_f.push(new_f_i);

            max_diff = max_diff.max((f[i] - new_f_i).len());
        }

        *f = new_f;

        if max_diff < 1e-8 {
            // println!("FINISHED at k = {k}: max_diff = {max_diff}");
            break;
        }
    }
}