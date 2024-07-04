use crate::z::EuclideanDivResult;
use core::cmp::Ordering;
use rug::integer::IntegerExt64;
use rug::ops::DivRounding;
use rug::ops::NegAssign;
use rug::rand::MutRandState;
use rug::Complete;
use rug::Integer;

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

    fn random<R: rand_core::CryptoRng>(_rng: &mut R) -> Self {
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
        let (q, r) = self.div_rem_floor_ref(other).complete();
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

    fn sample_bits<R: rand_core::CryptoRng + MutRandState>(nbits: u32, rng: &mut R) -> Self {
        <Self as From<_>>::from(Integer::random_bits(nbits, rng))
    }

    fn kronecker(&self, other: &Self) -> i32 {
        Integer::kronecker(self, other)
    }
}
