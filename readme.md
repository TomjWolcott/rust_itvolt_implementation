# Rust ITVOLT implementation
This is a small Rust implementation of ITVOLT to solve this simple 1D ODE numerically:
$$\left[ i \frac{\partial}{\partial t}-t \right ]\psi(t) = 0, \space \space \space \psi(0)=1$$

To run the code, first install [rust](https://www.rust-lang.org/learn/get-started) and then follow [these instructions](https://docs.rs/GSL/latest/rgsl/index.html) to download GSL.  Next clone this repository and run the paper inputs test using `cargo test test_paper_inputs -- --nocapture`.  You can enable optimizations on this code by using `cargo test test_paper_inputs --release -- --nocapture`.
