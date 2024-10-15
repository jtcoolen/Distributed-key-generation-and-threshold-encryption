// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

//! Arbitrary precision signed big integer (GMP backend).

use core::cmp::Ordering;
use std::cmp::Ordering::Less;
use std::ops::{AddAssign, ShrAssign};

#[cfg(feature = "random")]
use rand_core::{CryptoRng, RngCore};
use rug::integer::{IntegerExt64, Order};
use rug::ops::NegAssign;
use rug::ops::{DivRoundingAssign, RemRoundingAssign};
use rug::Complete;
pub use rug::Integer;

use crate::z::{EuclideanDivResult, Z};

#[cfg(feature = "random")]
impl crate::z::Randomizable for rug::Integer {
    fn random<R: rand_core::CryptoRng + rand_core::RngCore>(rng: &mut R) -> Self {
        let mut limbs = [0u64; 64];

        let mut dest = [0u8; 8 * 64];
        rng.fill_bytes(&mut dest);
        for (i, chunk) in dest.chunks_exact(8).enumerate() {
            limbs[i] = u64::from_le_bytes(chunk.try_into().expect("Chunk size is not 8"));
        }
        Integer::from_digits(&limbs, Order::Msf)
    }

    fn sample_bits<R: rand_core::CryptoRng + rand_core::RngCore>(nbits: u32, rng: &mut R) -> Self {
        // Calculate the number of bytes needed
        let nbytes = nbits / 8 + 1;

        // Create a buffer to hold the random bytes
        let mut buffer = vec![0u8; nbytes as usize];

        // Fill the buffer with random bytes
        rng.fill_bytes(&mut buffer);

        // Create a rug::Integer from the random bytes
        let mut integer = Integer::from_digits(&buffer, rug::integer::Order::Lsf);

        // Ensure the integer has exactly nbits by masking the most significant bits
        if nbits % 8 != 0 {
            let mask = (1u8 << (nbits % 8)) - 1;
            let last_byte_index = buffer.len() - 1;
            buffer[last_byte_index] &= mask;
            integer = Integer::from_digits(&buffer, rug::integer::Order::Lsf);
        }

        integer
    }

    fn sample_range<R: CryptoRng + RngCore>(rng: &mut R, min: &Self, max: &Self) -> Self
    where
        Self: Sized,
    {
        assert!(min <= max, "min should be less than or equal to max");
        let rd = Self::random(rng);
        rd.modulo_ref(&((max - min).complete() + (&<Self as From<_>>::from(1))))
            .complete()
            .add(min)
    }
}

impl crate::z::Z for rug::Integer {
    fn zero() -> Self {
        <Self as From<u64>>::from(0u64)
    }

    fn default() -> Self {
        Integer::new()
    }

    fn from(n: u64) -> Self {
        <Self as From<u64>>::from(n)
    }

    fn from_bytes_be(b: Vec<u8>, positive: bool) -> Self {
        let res = Integer::from_digits(&b, rug::integer::Order::Msf);
        if !positive {
            res.neg()
        } else {
            res
        }
    }

    fn to_bytes_be(&self) -> (Vec<u8>, bool) {
        (
            Integer::to_digits::<u8>(self, rug::integer::Order::Msf),
            self.is_positive(),
        )
    }

    fn eq_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Equal
    }

    fn strictly_less_than_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Less
    }

    fn strictly_less_than(&self, rhs: &Self) -> bool {
        self < rhs
    }

    fn add(&self, rhs: &Self) -> Self {
        (self + rhs).complete()
    }

    fn sub(&self, rhs: &Self) -> Self {
        (self - rhs).complete()
    }

    fn mul(&self, rhs: &Self) -> Self {
        (self * rhs).complete()
    }

    fn sqr(&self) -> Self {
        self.square_ref().complete()
    }

    fn neg(&self) -> Self {
        (-self).complete()
    }

    fn divide_by_2_exact(&mut self) {
        self.div_exact_u_mut(2)
    }

    fn divide_by_2(&mut self) {
        rug::Integer::shr_assign(self, 1u32)
    }

    fn divide_by_4_exact(&mut self) {
        self.div_exact_u_mut(4)
    }

    fn is_odd(&self) -> bool {
        self.is_odd()
    }

    fn euclidean_div_ceil(&self, other: &Self) -> crate::z::EuclideanDivResult<Self>
    where
        Self: Sized,
    {
        let (q, r) = self.div_rem_ceil_ref(other).complete();
        EuclideanDivResult {
            quotient: q,
            remainder: r,
        }
    }

    fn oppose(&mut self) {
        self.neg_assign()
    }

    fn is_positive(&self) -> bool {
        self.is_positive()
    }

    fn root(&self, n: u32) -> Self {
        self.root_ref(n).complete()
    }

    fn gcd(&self, other: &Self) -> Self {
        self.gcd_ref(other).complete()
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized,
    {
        self.extended_gcd_ref(other).complete()
    }

    fn divide_exact(&self, other: &Self) -> Self {
        self.div_exact_ref(other).complete()
    }

    fn divide_pow2(&mut self, other: usize) {
        Integer::shr_assign(self, other)
    }

    fn divides(&self, other: &Self) -> bool {
        other.is_divisible(self)
    }

    fn add_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self + other).complete().modulo(modulo)
    }

    fn sub_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self - other).complete().modulo_ref(modulo).complete()
    }

    fn take_mod(&self, modulo: &Self) -> Self {
        self.modulo_ref(modulo).complete()
    }

    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self * other).complete().modulo_ref(modulo).complete()
    }

    fn set_sign(&mut self, positive: bool) {
        match (self.is_positive(), positive) {
            (true, true) => (),
            (true, false) => self.oppose(),
            (false, true) => self.oppose(),
            (false, false) => (),
        }
    }

    fn div_floor(&self, other: &Self) -> Self {
        let mut x = self.clone();
        x.div_floor_assign(other);
        x
    }

    fn bit_size(&self) -> u64 {
        Integer::significant_bits_64(self)
    }

    fn get_bit(&self, index: u64) -> bool {
        Integer::get_bit_64(self, index)
    }

    fn sqrt(&self) -> Self {
        Integer::sqrt_ref(self).complete()
    }

    fn next_prime(&self) -> Self {
        Integer::next_prime_ref(self).complete()
    }

    fn kronecker(&self, other: &Self) -> i32 {
        Integer::kronecker(self, other)
    }

    fn invert_mod(&self, modulo: &Self) -> Option<Self>
    where
        Self: Sized,
    {
        self.invert_ref(modulo).map(|b| b.into())
    }

    fn compare(&self, other: &Self) -> Ordering {
        Integer::cmp(self, other)
    }

    fn sqrt_mod_prime(&self, prime: &Self) -> Option<Self> {
        // Shanks-Tonelli
        if self == &Integer::ZERO {
            return Some(Integer::ZERO);
        }

        // Check if a is a quadratic residue modulo p
        if self.legendre(prime) != 1 {
            return None;
        }

        if prime.take_mod(&4u32.into()) == 3u32 {
            let exp: Integer = (prime.add(&1.into())).div_exact_u(4);
            return Some(self.pow_mod_ref(&exp, prime).unwrap().into());
        }

        // Find n such that n is a quadratic non-residue modulo p
        let mut n: Integer = 2u32.into();
        while n.compare(prime) == Less {
            if n.legendre(prime) == -1 {
                break;
            }
            n += 1;
        }

        // p - 1 = 2^s * q
        let mut q = (prime - Integer::ONE).complete();
        let mut s: i32 = 0;
        while q.is_even() {
            q >>= 1;
            s += 1;
        }

        let inv_prime = self.invert_mod(prime).unwrap();

        let r =
            <Integer as From<_>>::from(self.pow_mod_ref(&((q.clone() + 1) >> 1), prime).unwrap());
        let y = r
            .sqr()
            .take_mod(prime)
            .mul_mod(&inv_prime, prime)
            .take_mod(prime);
        let b: Integer = n.pow_mod_ref(&q, prime).unwrap().into();
        let mut j = Integer::ZERO.clone();

        // Calculate the power j of b such that b^(2*j)*r^2/a = 1 mod p.
        for k in 0..=(s - 2) {
            let exp = <Integer as From<_>>::from(Integer::ONE << ((s - 2 - k) as u32));
            let b_pow =
                <Integer as From<_>>::from(b.pow_mod_ref(&(j.clone() << 1), prime).unwrap());
            let b_pow = <Integer as From<_>>::from(
                b_pow.mul_mod(&y, prime).pow_mod_ref(&exp, prime).unwrap(),
            );
            if b_pow != 1 {
                j.add_assign(<Integer as From<_>>::from(Integer::ONE << k));
            }
        }

        // b^(2*j)*r^2/a = 1 mod p => (b^j*r)^2 = a mod p.
        Some(
            <Integer as From<_>>::from(b.pow_mod_ref(&j, prime).unwrap())
                .mul(&r)
                .take_mod(prime),
        )
    }

    fn abs(&self) -> Self {
        Integer::abs_ref(self).complete()
    }

    fn from_i64(i: i64) -> Self {
        <Integer as From<_>>::from(i)
    }

    fn from_string(s: &str, base: u64) -> Self {
        Integer::from_str_radix(s, base as i32).unwrap()
    }

    fn rem_trunc_inplace(&mut self, other: &Self) {
        self.rem_trunc_assign(other)
    }

    fn div_rem(&self, _other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized,
    {
        todo!()
    }

    #[cfg(feature = "gmp")]
    fn ceil_abslog_square(&self) -> Self {
        use rug::{ops::CompleteRound, Float};
        // Determine the bit size of `n`
        let nbits = self.significant_bits();

        // Initialize precision for floating-point calculations
        let precision = nbits + 10; // Add some extra bits for safety
        let mut nf = Float::with_val(precision, self);
        nf.abs_mut();

        // Variables used in the computation
        let mut acc = Float::with_val(precision, 0);
        let mut z = Float::with_val(precision, 0);
        let mut pow = Float::with_val(precision, 0);
        let mut h = Float::with_val(precision, 1);
        let mut t = Float::with_val(precision, 0);
        let mut tmp = Float::with_val(precision, 0);

        // Compute `m` such that log(|n|)^2 = 4^m * log(nf)^2 with 0 < nf < 2
        let mut m = 0;
        while nf.cmp_abs(&Float::with_val(precision, 2)).unwrap().is_ge() {
            nf.sqrt_mut();
            m += 1;
        }

        // Calculate z = 1 - nf, with -1 < z < 1
        z = Float::with_val(precision, 1);
        z -= &nf;

        pow = z.clone().square();

        for i in 1..1_000_000 {
            // Max iterations: 1e6
            tmp = Float::with_val(precision, i + 1);
            tmp.recip_mut();

            t = (&pow * &h).complete(precision);
            t.shr_assign(1);
            t *= &tmp;
            acc += &t;

            pow *= &z;
            h += &tmp;

            tmp = t / &acc;
            tmp.abs_mut();
            if tmp
                .cmp_abs(&(Float::with_val(precision, 1) / Float::with_val(precision, 10)))
                .unwrap()
                .is_le()
            {
                // Convergence check
                break;
            }
        }

        acc.shr_assign(2 * m);
        acc.ceil_mut();
        acc.to_integer().unwrap()
    }

    fn incr(&self) -> Self {
        (self + Integer::ONE).complete()
    }

    fn shl(&self, n: u32) -> Self {
        (self << n).complete()
    }
}
