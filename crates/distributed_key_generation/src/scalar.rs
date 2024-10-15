// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use crate::z;
use blst::{blst_fr, blst_fr_from_scalar};
use ff::Field;
use ff::PrimeField;
use rand::RngCore;

const DOMAIN_SEPARATOR_TAG: &[u8; 35] = b"TEZOS_NIDKG_V1_BLS12_381_SCALAR_DST";

/// A trait defining operations for scalar field elements.
///
/// This trait provides methods for arithmetic operations, conversions,
/// and hashing of scalar field elements.
pub trait Scalar {
    /// Generates a random scalar field element.
    ///
    /// This method requires the `random` feature to be enabled.
    ///
    /// # Parameters
    ///
    /// - `rng`: A mutable reference to a random number generator.
    ///
    /// # Returns
    ///
    /// A random scalar field element.
    #[cfg(feature = "random")]
    fn random<R: RngCore>(rng: &mut R) -> Self;

    /// Adds another scalar field element to `self` in place.
    ///
    /// # Parameters
    ///
    /// - `rhs`: The scalar field element to add.
    fn add_assign(&mut self, rhs: &Self);

    /// Multiplies `self` by another scalar field element in place.
    ///
    /// # Parameters
    ///
    /// - `rhs`: The scalar field element to multiply by.
    fn mul_assign(&mut self, rhs: &Self);

    /// Computes `self` raised to the power of `exp`.
    ///
    /// # Parameters
    ///
    /// - `exp`: The exponent to raise `self` to.
    ///
    /// # Returns
    ///
    /// The result of raising `self` to the power of `exp`.
    fn pow(&self, exp: u64) -> Self;

    /// Returns the modulus of the scalar field as a `z::Z` type.
    ///
    /// # Type Parameters
    ///
    /// - `Z`: A type implementing the `z::Z` trait.
    ///
    /// # Returns
    ///
    /// The modulus of the scalar field represented as a `z::Z` type.
    fn modulus_as_z<Z: z::Z>() -> Z;

    /// Converts the scalar field element to a `z::Z` type.
    ///
    /// # Type Parameters
    ///
    /// - `Z`: A type implementing the `z::Z` trait.
    ///
    /// # Returns
    ///
    /// The scalar field element represented as a `z::Z` type.
    fn to_z<Z: z::Z>(&self) -> Z;

    /// Converts a byte vector into a scalar field element.
    ///
    /// # Parameters
    ///
    /// - `b`: A byte vector representing the scalar field element.
    ///
    /// # Returns
    ///
    /// The scalar field element corresponding to the byte vector.
    fn from_bytes_be(_: Vec<u8>) -> Self;

    /// Computes the square of `self`.
    ///
    /// # Returns
    ///
    /// The square of `self`.
    fn sqr(&self) -> Self;

    /// Returns the zero element of the scalar field.
    ///
    /// # Returns
    ///
    /// The zero scalar field element.
    fn zero() -> Self;

    /// Creates a scalar field element from a `u64` value.
    ///
    /// # Parameters
    ///
    /// - `n`: The `u64` value to convert to a scalar field element.
    ///
    /// # Returns
    ///
    /// The scalar field element corresponding to the `u64` value.
    fn from(n: u64) -> Self;

    /// Computes a scalar field element from a hash of the input byte vector.
    ///
    /// The hash is expected to be performed with an appropriate Domain Separator Tag.
    ///
    /// # Parameters
    ///
    /// - `b`: A byte vector to hash.
    ///
    /// # Returns
    ///
    /// The scalar field element derived from the hash of the input byte vector.
    fn hash(b: Vec<u8>) -> Self;
}

impl crate::scalar::Scalar for blstrs::Scalar {
    #[cfg(feature = "random")]
    fn random<R: RngCore + rand_core::RngCore>(rng: &mut R) -> Self {
        <blstrs::Scalar as Field>::random(rng)
    }

    fn add_assign(&mut self, rhs: &Self) {
        *self += rhs
    }

    fn mul_assign(&mut self, rhs: &Self) {
        *self *= rhs
    }

    fn pow(&self, exp: u64) -> Self {
        <blstrs::Scalar as Field>::pow(self, [exp])
    }

    fn modulus_as_z<Z: z::Z>() -> Z {
        Z::from_string(&blstrs::Scalar::MODULUS[2..], 16)
    }

    fn to_z<Z: z::Z>(&self) -> Z {
        Z::from_bytes_be(self.to_bytes_be().to_vec(), true)
    }

    fn from_bytes_be(b: Vec<u8>) -> Self {
        assert!(b.len() <= 32, "Input byte array is too long");
        let mut padded = [0u8; 32];
        padded[32 - b.len()..].copy_from_slice(&b);
        Self::from_bytes_be(&padded).unwrap()
    }

    fn sqr(&self) -> Self {
        self.square()
    }

    fn zero() -> Self {
        blstrs::Scalar::ZERO.clone()
    }

    fn from(n: u64) -> Self {
        <blstrs::Scalar as From<u64>>::from(n)
    }

    fn hash(b: Vec<u8>) -> Self {
        let hash = blst::blst_scalar::hash_to(&b, DOMAIN_SEPARATOR_TAG).unwrap();
        let mut ret = blst::blst_fr::default();
        unsafe {
            blst_fr_from_scalar(&mut ret, &hash);
        }
        <blstrs::Scalar as From<blst_fr>>::from(ret)
    }
}
