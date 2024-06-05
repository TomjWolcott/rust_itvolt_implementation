# Rust ITVOLT implementation
This is a small Rust implementation of ITVOLT to solve this simple 1D ODE numerically:
$$\left[ i \frac{\partial}{\partial t}-t \right ]\psi(t) = 0, \space \space \space \psi(0)=1$$

To run the code, first install [rust](https://www.rust-lang.org/learn/get-started) and then follow [these instructions](https://docs.rs/GSL/latest/rgsl/index.html) to download GSL.  Next clone this repository and run the paper inputs test using `cargo test test_paper_inputs -- --nocapture`.  You can enable optimizations on this code by using `cargo test test_paper_inputs --release -- --nocapture`.

## Output of `test_paper_inputs`:
```
t_α: 0, Δτ: 0.1, n: 6, type: Jacobi
    ε_sol: 9.658089e-5
    k_max: 26

t_α: 0, Δτ: 0.1, n: 12, type: Jacobi
    ε_sol: 1.024754e-10
    k_max: 20

t_α: 0, Δτ: 0.1, n: 18, type: Jacobi
    ε_sol: 1.018039e-10
    k_max: 20

t_α: 0, Δτ: 1, n: 15, type: Jacobi
    ε_sol: 4.459584e72
    k_max: 75

t_α: 0, Δτ: 1, n: 30, type: Jacobi
    ε_sol: 5.077035e-7
    k_max: 150

t_α: 0, Δτ: 1, n: 45, type: Jacobi
    ε_sol: 7.875728e-7
    k_max: 225

t_α: 0, Δτ: 0.1, n: 6, type: GaussSeidel
    ε_sol: 9.658072e-5
    k_max: 11

t_α: 0, Δτ: 0.1, n: 12, type: GaussSeidel
    ε_sol: 5.002168e-11
    k_max: 8

t_α: 0, Δτ: 0.1, n: 18, type: GaussSeidel
    ε_sol: 3.487181e-11
    k_max: 7

t_α: 0, Δτ: 1, n: 15, type: GaussSeidel
    ε_sol: 3.247329e-1
    k_max: 75

t_α: 0, Δτ: 1, n: 30, type: GaussSeidel
    ε_sol: 6.278638e-10
    k_max: 22

t_α: 0, Δτ: 1, n: 45, type: GaussSeidel
    ε_sol: 2.206529e-11
    k_max: 17

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 3, type: Jacobi
    ε_sol: 7.812500e-7
    k_max: 2

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 6, type: Jacobi
    ε_sol: 4.058449e-11
    k_max: 4

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 12, type: Jacobi
    ε_sol: 4.057995e-11
    k_max: 4

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 5, type: Jacobi
    ε_sol: 4.660516e-5
    k_max: 8

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 15, type: Jacobi
    ε_sol: 1.480701e-12
    k_max: 7

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 30, type: Jacobi
    ε_sol: 1.448017e-12
    k_max: 7

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 3, type: GaussSeidel
    ε_sol: 7.812500e-7
    k_max: 2

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 6, type: GaussSeidel
    ε_sol: 4.058430e-11
    k_max: 3

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 0.1, n: 12, type: GaussSeidel
    ε_sol: 4.057763e-11
    k_max: 3

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 5, type: GaussSeidel
    ε_sol: 4.660516e-5
    k_max: 5

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 15, type: GaussSeidel
    ε_sol: 3.069663e-14
    k_max: 5

t_α: (τ_j + τ_{j+1}) / 2, Δτ: 1, n: 30, type: GaussSeidel
    ε_sol: 5.240547e-14
    k_max: 4
```