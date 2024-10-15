// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use ff::Field;
use rand_core::{CryptoRng, RngCore};

/// A trait representing a polynomial with operations for managing coefficients and evaluation.
///
/// This trait provides methods for:
/// - Creating random polynomials
/// - Setting coefficients
/// - Evaluating the polynomial at a given point
/// - Retrieving coefficients of the polynomial
///
/// The type of the polynomial's coefficients is specified by the associated `Scalar` type.
pub trait Polynomial {
    /// The type of the scalar used for the polynomial's coefficients.
    type Scalar;

    /// Creates a random polynomial of a given degree.
    ///
    /// This method requires the `random` feature to be enabled.
    ///
    /// # Parameters
    ///
    /// - `rng`: A mutable reference to a random number generator.
    /// - `degree`: The degree of the polynomial to generate. This will result in a polynomial with
    ///   `degree + 1` coefficients.
    ///
    /// # Returns
    ///
    /// A polynomial represented as a vector of random coefficients.
    #[cfg(feature = "random")]
    fn random<R: RngCore + CryptoRng>(rng: &mut R, degree: usize) -> Self;

    /// Sets the coefficient of the polynomial at the specified index.
    ///
    /// # Parameters
    ///
    /// - `index`: The index of the coefficient to set. The index corresponds to the power of `x`
    ///   in the polynomial term.
    /// - `scalar`: The scalar value to set at the given index.
    ///
    /// # Panics
    ///
    /// Panics if `index` is out of bounds of the polynomial's coefficient vector.
    fn set_coefficient(&mut self, index: usize, scalar: &Self::Scalar);

    /// Evaluates the polynomial at a given point.
    ///
    /// This method uses Horner's method for evaluation.
    ///
    /// # Parameters
    ///
    /// - `point`: The point at which to evaluate the polynomial.
    ///
    /// # Returns
    ///
    /// The result of evaluating the polynomial at the given point.
    fn evaluate(&self, point: &Self::Scalar) -> Self::Scalar;

    /// Retrieves the coefficients of the polynomial.
    ///
    /// # Returns
    ///
    /// A vector containing the coefficients of the polynomial. The coefficients are ordered
    /// from the constant term (index 0) to the highest degree term.
    fn coefficients(&self) -> Vec<Self::Scalar>;
}

impl Polynomial for Vec<blstrs::Scalar> {
    type Scalar = blstrs::Scalar;

    #[cfg(feature = "random")]
    fn random<R: RngCore + CryptoRng>(rng: &mut R, degree: usize) -> Self {
        let mut res: Vec<blstrs::Scalar> = (0..=degree)
            .map(|_| crate::scalar::Scalar::random(rng))
            .collect();
        while blstrs::Scalar::ZERO == res[degree] {
            res[degree] = crate::scalar::Scalar::random(rng);
        }
        res
    }

    fn set_coefficient(&mut self, index: usize, scalar: &Self::Scalar) {
        self[index] = *scalar
    }

    fn evaluate(&self, point: &Self::Scalar) -> Self::Scalar {
        let mut result = crate::scalar::Scalar::zero();
        for coef in self.iter().rev() {
            result = result * point + coef;
        }
        result
    }

    fn coefficients(&self) -> Vec<Self::Scalar> {
        self.clone()
    }
}
