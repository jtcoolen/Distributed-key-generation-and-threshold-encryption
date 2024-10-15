// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

//! 4096-bit signed big integer (crypto-bigint backend).

use std::cmp::Ordering::Less;
use std::cmp::{Ordering, PartialEq};
use std::ops::{Add, Rem, Sub};

use crypto_bigint::{Encoding, Integer, NonZero, U4096};
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::z::{EuclideanDivResult, Z};

/// Signed integer with range (-2^4096, 2^4096)
/// Does not check for overflows!
#[derive(Debug, Clone)]
pub struct Bignum {
    pub positive: bool,
    pub uint: U4096,
}

impl PartialEq for Bignum {
    fn eq(&self, other: &Self) -> bool {
        (self.uint.eq(&other.uint)) && (self.positive == other.positive)
    }
}

impl Serialize for Bignum {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Use a Serde struct serializer to serialize the Bignum struct
        let mut state = serializer.serialize_struct("Bignum", 2)?;
        state.serialize_field("positive", &self.positive)?;

        // Serialize the uint field as a byte array
        let uint_bytes = self.uint.to_be_bytes();
        state.serialize_field("uint", &uint_bytes[..])?;

        state.end()
    }
}

impl PartialEq for EuclideanDivResult<Bignum> {
    fn eq(&self, other: &Self) -> bool {
        self.quotient.eq_abs(&other.quotient)
            && self.remainder.eq_abs(&other.remainder)
            && self.quotient.positive == other.quotient.positive
            && self.remainder.positive == other.remainder.positive
    }
}

impl Rem for &Bignum {
    type Output = Bignum;

    fn rem(self, rhs: Self) -> Self::Output {
        Bignum {
            positive: self.positive,
            uint: self.uint.rem_vartime(&NonZero::new(rhs.uint).unwrap()),
        }
    }
}

impl crate::z::Z for Bignum {
    fn zero() -> Self {
        Bignum {
            positive: true,
            uint: U4096::ZERO,
        }
    }

    fn default() -> Self {
        Bignum {
            positive: true,
            uint: U4096::ZERO,
        }
    }

    fn from(n: u64) -> Self {
        Bignum {
            positive: true,
            uint: <U4096 as From<_>>::from(n),
        }
    }

    fn from_string(_s: &str, _base: u64) -> Self {
        todo!()
    }

    fn from_i64(_i: i64) -> Self {
        todo!()
    }

    fn from_bytes_be(b: Vec<u8>, positive: bool) -> Self {
        let mut padded = [0u8; 512];
        let length = b.len();

        if length > 512 {
            panic!("Input vector is too large");
        }

        padded[512 - length..].copy_from_slice(&b);

        Bignum {
            positive,
            uint: U4096::from_be_bytes(padded),
        }
    }

    fn to_bytes_be(&self) -> (Vec<u8>, bool) {
        let b: Vec<u8> = U4096::to_be_bytes(&self.uint).into();
        (b, self.positive)
    }

    fn eq_abs(&self, rhs: &Self) -> bool {
        U4096::eq(&self.uint, &rhs.uint)
    }

    fn strictly_less_than_abs(&self, rhs: &Self) -> bool {
        U4096::cmp(&self.uint, &rhs.uint) == Less
    }

    fn strictly_less_than(&self, rhs: &Self) -> bool {
        let both_positive = self.positive && rhs.positive;
        let both_negative = !self.positive && !rhs.positive;
        let self_negative_rhs_positive = !self.positive && rhs.positive;
        let less_than_abs = self.strictly_less_than_abs(rhs);

        let result_if_both_positive = both_positive && less_than_abs;
        let result_if_both_negative = both_negative && !less_than_abs;
        let result_if_self_negative_rhs_positive = self_negative_rhs_positive;

        result_if_both_positive || result_if_both_negative || result_if_self_negative_rhs_positive
    }

    fn add(&self, rhs: &Self) -> Self {
        let mut res = Self::zero();
        match (
            self.positive,
            rhs.positive,
            self.strictly_less_than_abs(rhs),
        ) {
            (true, true, _) => {
                res.uint = <U4096 as Add>::add(self.uint, rhs.uint);
                res.positive = true;
            }
            (false, false, _) => {
                res.uint = <U4096 as Add>::add(self.uint, rhs.uint);
                res.positive = false;
            }
            (true, false, true) => {
                res.uint = <U4096 as Sub>::sub(rhs.uint, self.uint);
                res.positive = false;
            }
            (true, false, false) => {
                res.uint = <U4096 as Sub>::sub(self.uint, rhs.uint);
                res.positive = true;
            }
            (false, true, true) => {
                res.uint = <U4096 as Sub>::sub(rhs.uint, self.uint);
                res.positive = true;
            }
            (false, true, false) => {
                res.uint = <U4096 as Sub>::sub(self.uint, rhs.uint);
                res.positive = false;
            }
        }
        res
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    fn mul(&self, rhs: &Self) -> Self {
        let mut res = Self::zero();
        let (lo, _hi) = U4096::split_mul(&self.uint, &rhs.uint);
        res.uint = lo;
        res.positive = self.positive == rhs.positive;
        res
    }

    fn sqr(&self) -> Self {
        Bignum {
            positive: true,
            uint: self.uint.square_wide().0,
        }
    }

    fn neg(&self) -> Self {
        Bignum {
            positive: !self.positive,
            uint: self.uint,
        }
    }

    fn divide_by_2_exact(&mut self) {
        self.uint = self.uint.shr_vartime(1u32);
    }

    fn divide_by_4_exact(&mut self) {
        self.uint = self.uint.shr_vartime(2u32);
    }

    fn is_odd(&self) -> bool {
        crypto_bigint::Uint::is_odd(&self.uint).into()
    }

    fn euclidean_div_ceil(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized,
    {
        let divisor = NonZero::new(other.uint).unwrap();
        let (q, r) = U4096::div_rem_vartime(&self.uint, &divisor);

        let a_positive = self.positive;
        let b_positive = other.positive;
        let ab_positive = a_positive && b_positive;

        let mut quotient = Bignum {
            positive: ab_positive,
            uint: q,
        };
        let mut remainder = Bignum {
            positive: a_positive,
            uint: r,
        };

        if remainder.uint != U4096::ZERO
            && (self.compare(&Self::zero()).is_gt() && other.compare(&Self::zero()).is_gt()
                || (self.compare(&Self::zero()).is_lt() && other.compare(&Self::zero()).is_lt()))
        {
            quotient.uint = quotient.uint.add(&U4096::from(1u32));
            remainder = remainder.sub(other)
        }

        EuclideanDivResult {
            quotient,
            remainder,
        }
    }

    fn oppose(&mut self) {
        self.positive = !self.positive
    }

    fn is_positive(&self) -> bool {
        self.positive
    }

    fn root(&self, _n: u32) -> Self {
        todo!()
    }

    fn gcd(&self, other: &Self) -> Self {
        Bignum {
            positive: !(self.uint == U4096::ZERO && other.uint == U4096::ZERO),
            uint: U4096::gcd(&self.uint, &other.uint),
        }
    }

    fn div_rem(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized,
    {
        let divisor = NonZero::new(other.uint).unwrap();
        let (q, r) = U4096::div_rem_vartime(&self.uint, &divisor);
        EuclideanDivResult {
            quotient: Bignum {
                positive: true,
                uint: q,
            },
            remainder: Bignum {
                positive: true,
                uint: r,
            },
        }
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized,
    {
        let zero = <Bignum as crate::z::Z>::from(0);
        let one = <Bignum as crate::z::Z>::from(1);

        let mut s = zero.clone();
        let mut old_s = one;
        let mut r = other.clone();
        let mut old_r = self.clone();

        while !r.eq_abs(&zero) {
            let EuclideanDivResult { quotient, .. } = old_r.div_rem(&r);

            // Update old_r and r using `std::mem::swap` to avoid unnecessary cloning
            old_r = old_r.sub(&quotient.mul(&r));
            old_r = std::mem::replace(&mut r, old_r);

            // Update old_s and s in the same manner
            old_s = old_s.sub(&quotient.mul(&s));
            old_s = std::mem::replace(&mut s, old_s);
        }

        let bezout_t = if !other.eq_abs(&zero) {
            old_r.sub(&old_s.mul(self)).div_rem(other).quotient
        } else {
            zero
        };

        (old_r, old_s, bezout_t)
    }

    fn divide_exact(&self, other: &Self) -> Self {
        Bignum {
            positive: self.positive == other.positive,
            uint: U4096::wrapping_div_vartime(&self.uint, &NonZero::new(other.uint).unwrap()),
        }
    }

    fn divides(&self, _other: &Self) -> bool {
        todo!()
    }

    fn add_mod(&self, other: &Self, modulo: &Self) -> Self {
        assert!(modulo.is_positive());
        Self::add(self, other).take_mod(modulo)
    }

    fn sub_mod(&self, other: &Self, modulo: &Self) -> Self {
        assert!(modulo.is_positive());
        let mut res = &((self % modulo).sub(&(other % modulo))) % modulo;

        if !res.is_positive() {
            res = res.add(modulo)
        }
        res.take_mod(modulo)
    }

    fn take_mod(&self, modulo: &Self) -> Self {
        assert!(modulo.is_positive());
        Bignum {
            positive: true,
            uint: U4096::rem_vartime(&self.uint, &NonZero::new(modulo.uint).unwrap()),
        }
    }

    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self {
        assert!(modulo.is_positive());
        let mut res = &Self::mul(self, other) % modulo;
        if !res.is_positive() {
            res = res.add(modulo)
        }
        res.take_mod(modulo)
    }

    fn set_sign(&mut self, _positive: bool) {
        todo!()
    }

    fn div_floor(&self, _other: &Self) -> Self {
        todo!()
    }

    fn bit_size(&self) -> u64 {
        U4096::bits(&self.uint).into()
    }

    fn get_bit(&self, index: u64) -> bool {
        U4096::bit(&self.uint, index as u32).into()
    }

    fn sqrt(&self) -> Self {
        Bignum {
            positive: self.positive,
            uint: U4096::checked_sqrt(&self.uint).unwrap(),
        }
    }

    fn next_prime(&self) -> Self {
        todo!()
    }

    fn kronecker(&self, _other: &Self) -> i32 {
        todo!()
    }

    fn invert_mod(&self, modulo: &Self) -> Option<Self>
    where
        Self: Sized,
    {
        assert!(modulo.positive);
        let res: Option<U4096> = U4096::inv_mod(&self.uint, &modulo.uint).into();
        res.map(|uint| Bignum {
            positive: self.positive,
            uint,
        })
    }

    fn compare(&self, other: &Self) -> Ordering {
        match (self.positive, other.positive) {
            (true, true) => U4096::cmp(&self.uint, &other.uint),
            (false, false) => U4096::cmp(&other.uint, &self.uint),
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
        }
    }

    fn sqrt_mod_prime(&self, _prime: &Self) -> Option<Self>
    where
        Self: Sized,
    {
        todo!()
    }

    fn abs(&self) -> Self {
        todo!()
    }

    fn rem_trunc_inplace(&mut self, other: &Self) {
        let other_non_zero = NonZero::new(other.uint).expect("Non-zero value expected");
        self.uint %= other_non_zero;
    }

    fn divide_pow2(&mut self, other: usize) {
        self.uint = self.uint.shr_vartime(other as u32);
    }

    fn divide_by_2(&mut self) {
        let mut s = self.clone();
        s.uint >>= 1u32;
        // rounding toward -infinity
        if self.is_odd() && !self.positive {
            s.uint = s.uint.add(U4096::from(1u32));
        }
        *self = s;
    }

    #[cfg(feature = "gmp")]
    fn ceil_abslog_square(&self) -> Self {
        todo!()
    }

    fn shl(&self, n: u32) -> Self {
        let mut s = self.clone();
        s.uint <<= n;
        s
    }

    fn incr(&self) -> Self {
        let mut s = self.clone();
        s.uint = s.uint + U4096::one();
        s
    }
}
