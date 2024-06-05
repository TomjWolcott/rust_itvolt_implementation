mod complex_wrapper;
use std::ops::{Add, AddAssign, Div, Mul, Sub};

use complex_wrapper::*;

mod gauss_quadrature;
use gauss_quadrature::*;

use rgsl::integration::qk61;


fn main() {

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

fn weights(t_0: f64, quad_points: &Vec<f64>) -> Vec<Vec<f64>> {
    quad_points.iter().map(|&t_i| {
        (0..quad_points.len()).map(|j| {
            qk61(lagrange(&quad_points, j), t_0, t_i).0
        }).collect()
    }).collect()
}

#[test]
fn test_paper_inputs() {
    // Inputs from the ITVOLT Paper
    let inputs = [
        (false, ItvoltType::Jacobi, 0.1, 6),
        (false, ItvoltType::Jacobi, 0.1, 12),
        (false, ItvoltType::Jacobi, 0.1, 18),
        (false, ItvoltType::Jacobi, 1.0, 15),
        (false, ItvoltType::Jacobi, 1.0, 30),
        (false, ItvoltType::Jacobi, 1.0, 45),
        (false, ItvoltType::GaussSeidel, 0.1, 6),
        (false, ItvoltType::GaussSeidel, 0.1, 12),
        (false, ItvoltType::GaussSeidel, 0.1, 18),
        (false, ItvoltType::GaussSeidel, 1.0, 15),
        (false, ItvoltType::GaussSeidel, 1.0, 30),
        (false, ItvoltType::GaussSeidel, 1.0, 45),
        (true, ItvoltType::Jacobi, 0.1, 3),
        (true, ItvoltType::Jacobi, 0.1, 6),
        (true, ItvoltType::Jacobi, 0.1, 12),
        (true, ItvoltType::Jacobi, 1.0, 5),
        (true, ItvoltType::Jacobi, 1.0, 15),
        (true, ItvoltType::Jacobi, 1.0, 30),
        (true, ItvoltType::GaussSeidel, 0.1, 3),
        (true, ItvoltType::GaussSeidel, 0.1, 6),
        (true, ItvoltType::GaussSeidel, 0.1, 12),
        (true, ItvoltType::GaussSeidel, 1.0, 5),
        (true, ItvoltType::GaussSeidel, 1.0, 15),
        (true, ItvoltType::GaussSeidel, 1.0, 30),
    ];

    for (use_midpoint, ty, delta_time, num_quad_points) in inputs.iter() {
        let (max_err, max_k) = simple_one_dimensional_ode(*use_midpoint, *delta_time, *num_quad_points, *ty);
        println!(
            "t_α: {}, Δτ: {delta_time}, n: {num_quad_points}, type: {ty:?}\n    ε_sol: {max_err:.6e}\n    k_max: {max_k}\n",
            if *use_midpoint { "(τ_j + τ_{j+1}) / 2" } else { "0" },
        );
    }
}

#[derive(Debug, Clone, Copy)]
pub enum ItvoltType {
    Jacobi,
    GaussSeidel
}

fn simple_one_dimensional_ode(
    use_midpoint: bool,        // t_α = 0 or (τ_j + τ_{j+1}) / 2
    delta_time: f64,           // Δτ
    num_quad_points: usize,    // n
    ty: ItvoltType
) -> (f64, usize) {
    let total_time = 25.0;
    let max_iterations = 5 * num_quad_points;
    let mut quad_points = gauss_lobatto_quad_points(num_quad_points, 0.0, delta_time); // t_i
    let weights = weights(0.0, &quad_points); // w_{i,j}
    let mut max_err = 0.0; // ε_{sol}
    let mut max_k = 0;    // k_{max}

    let mut f: Vec<Complex> = Vec::new();

    for j in 0..((total_time / delta_time) as usize) {
        let tau_start = j as f64 * delta_time;             // τ_j
        let tau_end = (j + 1) as f64 * delta_time;         // τ_{j+1}
        let t_a = if use_midpoint { (tau_start + tau_end) / 2.0 } else { 0.0 };       // t_α
        let last_f = f.last().cloned().unwrap_or(ONE); // f(τ_j) or f(0) = 1
        let g = |t: f64| E.pow(-I * t_a * (t - tau_start)) * last_f; // g(t)
        let kernel = |t_i: f64, t_j: f64| -I * E.pow(-I * t_a * (t_i - t_j)) * (t_j - t_a); // K(t_i, t_j)

        f = quad_points.iter().map(|&t_i| g(t_i)).collect();

        max_k = max_k.max(match ty {
            ItvoltType::Jacobi =>
                iterative_volterra_propagator_jacobi(&quad_points, &weights, &mut f, g, kernel, max_iterations),
            ItvoltType::GaussSeidel => {
                iterative_volterra_propagator_gauss_seidel(&quad_points, &weights, &mut f, g, kernel, max_iterations)
            }
        });

        for (&f_i, &t_i) in f.iter().zip(quad_points.iter()) {
            let correct = E.pow(-I * t_i.powi(2) / 2.0);
            max_err = (f_i - correct).magnitude().max(max_err);
        }

        quad_points.iter_mut().for_each(|t| *t += delta_time);
    }

    (max_err, max_k)
}

pub trait NumericalError {
    fn error(&self, other: &Self) -> f64;
}

impl NumericalError for Complex {
    fn error(&self, other: &Self) -> f64 {
        (*self - *other).magnitude()
    }
}

impl NumericalError for f64 {
    fn error(&self, other: &Self) -> f64 {
        (self - other).abs()
    }
}

fn iterative_volterra_propagator_jacobi<T, S>(
    quad_points: &Vec<f64>,
    weights: &Vec<Vec<f64>>,
    f: &mut Vec<T>,
    g: impl Fn(f64) -> T,
    kernel: impl Fn(f64, f64) -> S,
    max_iterations: usize
) -> usize where
    T: Add<T, Output = T> + AddAssign<T> + NumericalError + Copy,
    S: Mul<f64, Output = S> + Mul<T, Output = T> + Copy
{
    for k in 0..max_iterations {
        let mut max_diff: f64 = 0.0;
        let mut new_f = Vec::new(); // f^{(k+1)}

        for (i, &t_i) in quad_points.iter().enumerate() {
            let mut new_f_i = g(t_i);

            for (j, &t_j) in quad_points.iter().enumerate() {
                new_f_i += kernel(t_i, t_j) * weights[i][j] * f[j]
            }

            new_f.push(new_f_i);

            max_diff = max_diff.max(f[i].error(&new_f_i));
        }

        *f = new_f;

        if max_diff < 1e-10 {
            return k + 1;
        }
    }

    max_iterations
}

trait Identity {
    fn identity() -> Self;
}

impl Identity for Complex {
    fn identity() -> Self {
        ONE
    }
}

impl Identity for f64 {
    fn identity() -> Self {
        1.0
    }
}

fn iterative_volterra_propagator_gauss_seidel<T, S>(
    quad_points: &Vec<f64>,
    weights: &Vec<Vec<f64>>,
    f: &mut Vec<T>,
    g: impl Fn(f64) -> T,
    kernel: impl Fn(f64, f64) -> S,
    max_iterations: usize
) -> usize where
    T: Add<T, Output = T> + AddAssign<T> + Div<S, Output = T> + NumericalError + Copy,
    S: Mul<f64, Output = S> + Mul<T, Output = T> + Sub<S, Output = S> + Copy + Identity
{
    for k in 0..max_iterations {
        let mut max_diff: f64 = 0.0;
        let mut new_f = Vec::new(); // f^{(k+1)}

        for (i, &t_i) in quad_points.iter().enumerate() {
            let mut new_f_i = g(t_i);

            for j in 0..i {
                let t_j = quad_points[j];
                new_f_i += kernel(t_i, t_j) * weights[i][j] * new_f[j]
            }

            for j in i+1..quad_points.len() {
                let t_j = quad_points[j];
                new_f_i += kernel(t_i, t_j) * weights[i][j] * f[j]
            }

            new_f_i = new_f_i / (S::identity() - kernel(t_i, t_i) * weights[i][i]);

            new_f.push(new_f_i);

            max_diff = max_diff.max(f[i].error(&new_f_i));
        }

        *f = new_f;

        if max_diff < 1e-10 {
            return k + 1;
        }
    }

    max_iterations
}