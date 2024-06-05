use rgsl::types::ComplexF64;
use std::{iter::Sum, ops::*};

#[derive(Clone, Copy)]
pub struct Complex(ComplexF64);

pub const I: Complex = Complex(ComplexF64 { dat: [0.0, 1.0] });
pub const E: Complex = Complex(ComplexF64 {
    dat: [std::f64::consts::E, 0.0],
});
pub const ONE: Complex = Complex(ComplexF64 { dat: [1.0, 0.0] });

impl Neg for Complex {
    type Output = Self;

    fn neg(self) -> Self {
        -1.0 * self
    }
}

impl Mul<f64> for Complex {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self(self.0.mul_real(rhs))
    }
}

impl Mul<Complex> for f64 {
    type Output = Complex;

    fn mul(self, rhs: Complex) -> Complex {
        rhs * self
    }
}

impl Div<f64> for Complex {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self(self.0.div_real(rhs))
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0.add(&rhs.0))
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl Sum for Complex {
    fn sum<I: Iterator<Item = Complex>>(iter: I) -> Self {
        iter.fold(Complex(ComplexF64 { dat: [0.0, 0.0] }), |a, b| a + b)
    }
}

impl Sub for Complex {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0.sub(&rhs.0))
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self(self.0.mul(&rhs.0))
    }
}

impl Div for Complex {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self(self.0.div(&rhs.0))
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({:.4} + {:.4}i)", self.0.real(), self.0.imaginary())
    }
}

impl Complex {
    pub fn pow(self, n: Self) -> Self {
        Self(self.0.pow(&n.0))
    }

    pub fn len(self) -> f64 {
        self.0.mul(&self.0.conjugate()).real().sqrt()
    }

    pub fn exp(self) -> Self {
        Self(self.0.exp())
    }

    pub fn sin(self) -> Self {
        Self(self.0.sin())
    }

    pub fn cos(self) -> Self {
        Self(self.0.cos())
    }

    pub fn tan(self) -> Self {
        Self(self.0.tan())
    }

    pub fn sinh(self) -> Self {
        Self(self.0.sinh())
    }

    pub fn cosh(self) -> Self {
        Self(self.0.cosh())
    }

    pub fn tanh(self) -> Self {
        Self(self.0.tanh())
    }

    pub fn sqrt(self) -> Self {
        Self(self.0.sqrt())
    }

    pub fn log10(self) -> Self {
        Self(self.0.log10())
    }

    pub fn arg(self) -> f64 {
        self.0.arg()
    }

    pub fn real(self) -> f64 {
        self.0.real()
    }

    pub fn imag(self) -> f64 {
        self.0.imaginary()
    }
}

impl Deref for Complex {
    type Target = ComplexF64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Complex {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
