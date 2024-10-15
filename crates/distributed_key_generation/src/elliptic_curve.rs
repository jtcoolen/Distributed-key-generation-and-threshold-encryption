// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use group::Group;

/// A trait representing an elliptic curve group, where the elements of the curve can be
/// combined using group operations such as addition and scalar multiplication.
pub trait EllipticCurve {
    /// The associated scalar type used for scalar multiplication. This can be any arbitrary-sized integer,
    /// but depending on the implementation, it may be fixed-sized integer whose size matches the
    /// bit size of the order of the elliptic curve.
    type S;

    /// Returns the generator point of the elliptic curve. The generator is a fixed
    /// point on the curve used as a base for scalar multiplication.
    fn generator() -> Self;

    /// Adds another point on the elliptic curve to `self`, modifying `self` in place.
    ///
    /// # Arguments
    /// * `other` - A reference to another point on the elliptic curve to be added to `self`.
    fn add_assign(&mut self, other: &Self);

    /// Multiplies `self` by a scalar, modifying `self` in place.
    ///
    /// # Arguments
    /// * `scalar` - A reference to a scalar value used for multiplication.
    fn mul_assign(&mut self, scalar: &Self::S);

    /// Converts the elliptic curve point into a byte array representation.
    ///
    /// # Returns
    /// A `Vec<u8>` containing the byte representation of the point.
    fn to_bytes(&self) -> Vec<u8>;

    /// Compares `self` with another point on the elliptic curve for equality.
    ///
    /// # Arguments
    /// * `other` - A reference to another point on the elliptic curve.
    ///
    /// # Returns
    /// True if and only if `self` and `other` are equal.
    fn equals(&self, other: &Self) -> bool;

    /// Returns the neutral element (also known as the zero or identity element) of the elliptic curve.
    /// This is the point that acts as the additive identity in the group, meaning that adding this
    /// point to any other point on the curve results in the original point.
    fn zero() -> Self;

    /// Computes a multi-exponentiation, i.e., the sum of multiple points each multiplied by their
    /// corresponding scalar.
    ///
    /// # Arguments
    /// * `points` - A slice of points on the elliptic curve.
    /// * `scalars` - A slice of scalar values corresponding to the points.
    ///
    /// # Returns   
    /// A single elliptic curve point representing the sum of the pointwise products of the points and scalars.
    fn multiexp(points: &[Self], scalars: &[Self::S]) -> Self
    where
        Self: Sized;
}

/// Implementation of the `EllipticCurve` trait for the `blstrs::G1Projective` type, which represents
/// a point on the G1 curve in the BLS12-381 pairing-friendly elliptic curve.
impl EllipticCurve for blstrs::G1Projective {
    type S = blstrs::Scalar;

    fn generator() -> Self {
        group::Group::generator()
    }

    fn add_assign(&mut self, other: &Self) {
        *self += other
    }

    fn mul_assign(&mut self, scalar: &Self::S) {
        *self *= scalar;
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.to_compressed().to_vec()
    }

    fn equals(&self, other: &Self) -> bool {
        self.eq(other)
    }

    fn zero() -> Self {
        blstrs::G1Projective::identity()
    }

    fn multiexp(points: &[Self], scalars: &[Self::S]) -> Self {
        blstrs::G1Projective::multi_exp(points, scalars)
    }
}
