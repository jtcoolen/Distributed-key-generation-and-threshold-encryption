use core::cmp::Ordering;
use std::cmp::Ordering::Less;
use std::ops::AddAssign;

use rug::integer::IntegerExt64;
use rug::ops::DivRounding;
use rug::ops::NegAssign;
use rug::Complete;
use rug::Integer;

use crate::z::EuclideanDivResult;

// special compilation suite for wasm with emscripten
// see https://github.com/cryspen/libcrux/blob/0f74f5a6fa7477e5dd553d6d589ab04c6fcca2bd/sys/hacl/wasm.sh

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

    fn from_bytes(b: Vec<u8>) -> Self {
        Integer::from_digits(&b, rug::integer::Order::Lsf)
    }

    fn to_bytes(&self) -> Vec<u8> {
        Integer::to_digits::<u8>(self, rug::integer::Order::Lsf)
    }

    fn random<R: rand_core::CryptoRng + rand_core::RngCore>(_rng: &mut R) -> Self {
        let mut rand = rug::rand::RandState::new();
        Integer::random_bits(4096, &mut rand).into()
    }

    fn eq_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Equal
    }

    fn eq(&self, rhs: &Self) -> bool {
        self == rhs
    }

    fn less_than_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Less
    }

    fn less_than(&self, rhs: &Self) -> bool {
        self <= rhs
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

    fn divide_by_2(&mut self) {
        self.div_exact_u_mut(2)
    }

    fn divide_by_4(&mut self) {
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

    fn div_floor(&self, other: Self) -> Self {
        DivRounding::div_floor(self, other)
    }

    fn bit_size(&self) -> u64 {
        Integer::significant_bits_64(&self)
    }

    fn get_bit(&self, index: u64) -> bool {
        Integer::get_bit_64(&self, index)
    }

    fn sqrt(&self) -> Self {
        Integer::sqrt_ref(self).complete()
    }

    fn next_prime(&self) -> Self {
        Integer::next_prime_ref(&self).complete()
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

        //buffer.reverse();
        // Ensure the integer has exactly nbits by masking the most significant bits
        if nbits % 8 != 0 {
            let mask = (1u8 << (nbits % 8)) - 1;
            let last_byte_index = buffer.len() - 1;
            buffer[last_byte_index] &= mask;
            integer = Integer::from_digits(&buffer, rug::integer::Order::Lsf);
        }

        println!("nbits: {}, buffer={:?}", nbits, buffer);

        integer
    }

    fn kronecker(&self, other: &Self) -> i32 {
        Integer::kronecker(self, other)
    }

    fn invert_mod(&self, modulo: &Self) -> Option<Self>
    where
        Self: Sized,
    {
        self.invert_ref(modulo).and_then(|b| Some(b.into()))
    }

    fn compare(&self, other: &Self) -> Ordering {
        Integer::cmp(self, other)
    }

    fn remove(&self, factor: &Self) -> (Self, u32)
    where
        Self: Sized,
    {
        Integer::remove_factor_ref(self, factor).complete()
    }

    fn sqrt_mod_prime(&self, prime: &Self) -> Option<Self> {
        // Shanks-Tonelli
        // TODO clean this code
        if self == &Integer::ZERO {
            return Some(Integer::ZERO);
        }

        // Check if a is a quadratic residue modulo p
        if self.legendre(prime) != 1 {
            return None;
        }

        if prime.take_mod(&<Integer as From<u32>>::from(4u32))
            == (<Integer as From<u32>>::from(3u32))
        {
            let mut exp: Integer = <Integer as From<_>>::from(prime + 1).div_exact_u(4);
            return Some(<Integer as From<_>>::from(
                self.pow_mod_ref(&exp, prime).unwrap(),
            ));
        }

        // Find n such that n is a quadratic non-residue modulo p
        let mut n = <Integer as From<u32>>::from(2u32);
        while n.compare(prime) == Less {
            if n.legendre(prime) == -1 {
                break;
            }
            n += 1;
        }

        // p - 1 = 2^s * q
        let mut q = prime.sub(Integer::ONE);
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
        let mut b = <Integer as From<_>>::from(n.pow_mod_ref(&q, prime).unwrap());
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

    fn pow_mod(&self, exponent: &Self, modulo: &Self) -> Self {
        Integer::pow_mod_ref(self, exponent, modulo).unwrap().into()
    }

    fn from_i64(i: i64) -> Self {
        <Integer as From<_>>::from(i)
    }

    fn from_string(s: &str, base: u64) -> Self {
        Integer::from_str_radix(s, base as i32).unwrap()
    }

    fn from_bytes_be(b: Vec<u8>) -> Self {
        Integer::from_bytes(b)
    }
}
