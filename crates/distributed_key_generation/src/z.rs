// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

//! Signed big integer trait.

use std::cmp::Ordering;

#[cfg(feature = "random")]
use rand_core::CryptoRng;
use rug::integer::{Order, UnsignedPrimitive};

#[derive(Debug)]
pub struct EuclideanDivResult<Z> {
    pub quotient: Z,
    pub remainder: Z,
}

/// A trait defining operations for arbitrary-precision integers.
///
/// This trait provides a wide range of methods for arithmetic, modular arithmetic,
/// and number theory operations on arbitrary-precision integers.
///
/// Invalid values for the arguments of the functions below can trigger exceptions.
pub trait Z {
    fn zero() -> Self;


     fn from_digits<T>(digits: &[T], order: Order) -> Self
    where
        T: UnsignedPrimitive;

    /// Returns the default value of the integer type, typically zero.
    fn default() -> Self;

    fn to_string(&self) -> String;

    fn from(n: u64) -> Self
    where
        Self: Sized;

    /// Constructs an integer from a string representation.
    ///
    /// # Arguments
    ///
    /// * `s` - A string slice representing the number.
    /// * `base` - The base of the number system (e.g., 10 for decimal, 16 for hexadecimal).
    ///
    /// # Returns
    ///
    /// An integer representation of the string.
    fn from_string(s: &str, base: u64) -> Self;

    fn from_i64(i: i64) -> Self;

    /// Constructs an integer from a vector of bytes in big-endian order.
    ///
    /// # Arguments
    ///
    /// * `b` - A vector of bytes representing the integer.
    /// * `positive` - A boolean indicating if the integer is positive.
    ///
    /// # Returns
    ///
    /// An integer representation of the byte vector.
    fn from_bytes_be(b: Vec<u8>, positive: bool) -> Self;

    /// Converts the integer to a vector of bytes in big-endian order.
    ///
    /// # Returns
    ///
    /// A tuple containing a vector of bytes and a boolean indicating if the integer is positive.
    fn to_bytes_be(&self) -> (Vec<u8>, bool);

    /// Compares the absolute values of two integers for equality.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to compare against.
    ///
    /// # Returns
    ///
    /// `true` if the absolute values of `self` and `rhs` are equal, otherwise `false`.
    fn eq_abs(&self, rhs: &Self) -> bool;

    /// Checks if the absolute value of the integer is strictly less than that of another integer.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to compare against.
    ///
    /// # Returns
    ///
    /// `true` if the absolute value of `self` is strictly less than `rhs`,
    /// otherwise `false`.
    fn strictly_less_than_abs(&self, rhs: &Self) -> bool;

    /// Checks if the integer is strictly less than another integer.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to compare against.
    ///
    /// # Returns
    ///
    /// `true` if `self` is strictly less than `rhs`, otherwise `false`.
    fn strictly_less_than(&self, rhs: &Self) -> bool;

    /// Adds another integer to this integer.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to add.
    ///
    /// # Returns
    ///
    /// The result of the addition.
    fn add(&self, rhs: &Self) -> Self;

    /// Subtracts another integer from this integer.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to subtract.
    ///
    /// # Returns
    ///
    /// The result of the subtraction.
    fn sub(&self, rhs: &Self) -> Self;

    /// Multiplies this integer by another integer.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The integer to multiply with.
    ///
    /// # Returns
    ///
    /// The result of the multiplication.
    fn mul(&self, rhs: &Self) -> Self;

    /// Squares this integer.
    ///
    /// # Returns
    ///
    /// The result of squaring the integer.
    fn sqr(&self) -> Self;

    /// Negates this integer.
    ///
    /// # Returns
    ///
    /// The result of negating the integer.
    fn neg(&self) -> Self;

    /// Performs an exact division of the integer by 2.
    fn divide_by_2_exact(&mut self);

    /// Performs an exact division of the integer by 4.
    fn divide_by_4_exact(&mut self);

    /// Checks if the integer is odd.
    ///
    /// # Returns
    ///
    /// `true` if the integer is odd, otherwise `false`.
    fn is_odd(&self) -> bool;

    /// Performs Euclidean division with rounding towards +infinity.
    ///
    /// The remainder gets the opposite sign of the denominator.
    ///
    /// # Arguments
    ///
    /// * `other` - The divisor integer.
    ///
    /// # Returns
    ///
    /// A result containing the quotient and remainder.
    fn euclidean_div_ceil(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized;

    /// Performs division and returns both quotient and remainder.
    ///
    /// # Arguments
    ///
    /// * `other` - The divisor integer.
    ///
    /// # Returns
    ///
    /// A result containing the quotient and remainder.
    fn div_rem(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized;

    /// Negates the integer in place.
    fn oppose(&mut self);

    /// Checks if the integer is positive.
    ///
    /// # Returns
    ///
    /// `true` if the integer is positive, otherwise `false`.
    fn is_positive(&self) -> bool;

    /// Computes the n-th root of the integer and truncates it.
    ///
    /// # Arguments
    ///
    /// * `n` - The degree of the root to compute.
    ///
    /// # Returns
    ///
    /// The n-th root of the integer, truncated.
    fn root(&self, n: u32) -> Self;

    /// Computes the greatest common divisor (GCD) of this integer and another integer.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to compute the GCD with.
    ///
    /// # Returns
    ///
    /// The GCD of the two integers.
    fn gcd(&self, other: &Self) -> Self;

    /// Computes the extended GCD of this integer and another integer.
    ///
    /// This method returns a tuple `(gcd, x, y)` such that `gcd = x * self + y * other`.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to compute the extended GCD with.
    ///
    /// # Returns
    ///
    /// A tuple containing the GCD and the coefficients x and y: `(gcd, x, y)`.
    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized;

    /// Divides this integer by another integer and returns the result.
    ///
    /// # Arguments
    ///
    /// * `other` - The divisor integer.
    ///
    /// # Returns
    ///
    /// The result of the division.
    fn divide_exact(&self, other: &Self) -> Self;

    /// Checks if this integer divides another integer exactly.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to check divisibility.
    ///
    /// # Returns
    ///
    /// `true` if `self` divides `other` exactly, otherwise `false`.
    fn divides(&self, other: &Self) -> bool;

    /// Adds another integer modulo a given modulo.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to add.
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// The result of the addition modulo `modulo`.
    fn add_mod(&self, other: &Self, modulo: &Self) -> Self;

    /// Subtracts another integer modulo a given modulo.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to subtract.
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// The result of the subtraction modulo `modulo`.
    fn sub_mod(&self, other: &Self, modulo: &Self) -> Self;

    /// Takes the modulo of this integer.
    ///
    /// # Arguments
    ///
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// The result of taking `self` modulo `modulo`.
    fn take_mod(&self, modulo: &Self) -> Self;

    /// Multiplies this integer by another integer modulo a given modulo.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to multiply.
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// The result of the multiplication modulo `modulo`.
    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self;

    /// Sets the sign of the integer.
    ///
    /// # Arguments
    ///
    /// * `positive` - A boolean indicating if the integer should be positive.
    fn set_sign(&mut self, positive: bool);

    /// Solves a congruence equation `self * x ≡ other (mod modulo)` for `x`.
    ///
    /// # Arguments
    ///
    /// * `other` - The right-hand side of the congruence.
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// A tuple `(x, factor)` where `x` is a solution to the congruence and `factor` is the factor
    /// by which solutions differ.
    ///
    /// # Panics
    ///
    /// Panics if there is no solution to the congruence.
    fn solve_congruence(&self, other: &Self, modulo: &Self) -> (Self, Self)
    where
        Self: Sized,
    {
        let (gcd, s, _) = self.extended_gcd(modulo);
        let r = other.take_mod(&gcd);
        // a solution exists iff other is divisible by the GCD of self and modulo
        if !r.eq_abs(&Self::zero()) {
            panic!("no solution");
        };
        // The solutions are of the form
        // for i=1 to gcd-1: (other/gcd)*s+i*(modulo/gcd)
        (
            other.divide_exact(&gcd).mul_mod(&s, modulo),
            modulo.divide_exact(&gcd),
        )
    }

    /// Divides this integer by another integer and rounds towards negative infinity.
    ///
    /// # Arguments
    ///
    /// * `other` - The divisor integer.
    ///
    /// # Returns
    ///
    /// The result of the division rounded towards negative infinity.
    fn div_floor(&self, other: &Self) -> Self;

    /// Returns the bit size of the integer.
    ///
    /// # Returns
    ///
    /// The number of bits needed to represent the integer.
    fn bit_size(&self) -> u64;

    /// Gets the bit value at a specific index.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the bit to retrieve.
    ///
    /// # Returns
    ///
    /// `true` if the bit is set, otherwise `false`.
    fn get_bit(&self, index: u64) -> bool;

    /// Computes the square root of the integer.
    ///
    /// # Returns
    ///
    /// The integer square root of the integer.
    fn sqrt(&self) -> Self;

    /// Finds the next prime number greater than the integer.
    ///
    /// # Returns
    ///
    /// The next prime number greater than the integer.
    fn next_prime(&self) -> Self;

    /// Computes the Kronecker symbol of two integers.
    /// See [https://mathworld.wolfram.com/KroneckerSymbol.html].
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to compute the Kronecker symbol with.
    ///
    /// # Returns
    ///
    /// The Kronecker symbol of the two integers.
    fn kronecker(&self, other: &Self) -> i32;

    /// Computes the modular inverse of this integer with respect to a given modulus.
    ///
    /// # Arguments
    ///
    /// * `modulo` - The modulus.
    ///
    /// # Returns
    ///
    /// An `Option` containing the modular inverse if it exists, otherwise `None`.
    fn invert_mod(&self, modulo: &Self) -> Option<Self>
    where
        Self: Sized;

    /// Compares this integer with another integer.
    ///
    /// # Arguments
    ///
    /// * `other` - The integer to compare against.
    ///
    /// # Returns
    ///
    /// An `Ordering` indicating the relative order of the two integers.
    fn compare(&self, other: &Self) -> Ordering;

    /// Computes the square root modulo a prime number.
    ///
    /// # Arguments
    ///
    /// * `prime` - The prime modulus.
    ///
    /// # Returns
    ///
    /// An `Option` containing the square root modulo `prime` if it exists, otherwise `None`.
    fn sqrt_mod_prime(&self, prime: &Self) -> Option<Self>
    where
        Self: Sized;

    /// Computes the absolute value of the integer.
    ///
    /// # Returns
    ///
    /// The absolute value of the integer.
    fn abs(&self) -> Self;

    /// Computes the remainder of the division by another integer and modifies the integer in place.
    ///
    /// # Arguments
    ///
    /// * `other` - The divisor integer.
    fn rem_trunc_inplace(&mut self, other: &Self);

    /// Divides the integer by a power of 2 and modifies it in place.
    ///
    /// # Arguments
    ///
    /// * `other` - The power of 2 to divide by.
    fn divide_pow2(&mut self, other: usize);

    /// Divides the integer by 2 and modifies it in place.
    fn divide_by_2(&mut self);

    #[cfg(feature = "gmp")]
    /// Computes the ceiling of the absolute value of the logarithm of `self` squared.
    ///
    /// # Returns
    ///
    /// The ceiling of the absolute value of the logarithm of `self` squared.
    fn ceil_abslog_square(&self) -> Self;

    /// Performs a left bit shift on the integer.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of positions to shift.
    ///
    /// # Returns
    ///
    /// The result of the left bit shift.
    fn shl(&self, n: u32) -> Self;

    /// Increments `self` by one.
    fn incr(&self) -> Self;
}

#[cfg(feature = "random")]
pub trait Randomizable: Ord {
    fn random<R: CryptoRng + rand_core::RngCore>(rng: &mut R) -> Self;

    fn sample_bits<R: rand_core::CryptoRng + rand_core::RngCore>(nbits: u32, rng: &mut R) -> Self;

    fn sample_range<R: CryptoRng + rand_core::RngCore>(rng: &mut R, min: &Self, max: &Self) -> Self
    where
        Self: Sized;
}
