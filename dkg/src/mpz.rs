use crate::z::EuclideanDivResult;
use core::cmp::Ordering;
use rug::ops::DivRounding;
use rug::ops::NegAssign;
use rug::Complete;
use rug::Integer;

impl crate::z::Z for rug::Integer {
    fn default() -> Self {
        Integer::new()
    }

    fn add(&self, rhs: &Self) -> Self {
        (self + rhs).complete()
    }

    fn add_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self + other).complete().modulo(modulo)
    }

    fn clone(&self) -> Self {
        Clone::clone(&self)
    }

    fn divide_by_2(&mut self) {
        self.div_exact_u_mut(2)
    }

    fn divide_by_4(&mut self) {
        self.div_exact_u_mut(4)
    }

    fn divide_exact(&self, other: &Self) -> Self {
        self.div_exact_ref(other).complete()
    }

    fn divides(&self, other: &Self) -> bool {
        other.is_divisible(self)
    }

    fn eq(&self, rhs: &Self) -> bool {
        self == rhs
    }

    fn eq_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Equal
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

    fn gcd(&self, other: &Self) -> Self {
        self.gcd_ref(other).complete()
    }

    fn is_odd(&self) -> bool {
        self.is_odd()
    }

    fn is_positive(&self) -> bool {
        self.is_positive()
    }

    fn mul(&self, rhs: &Self) -> Self {
        (self * rhs).complete()
    }

    fn less_than(&self, rhs: &Self) -> bool {
        self <= rhs
    }

    fn oppose(&mut self) {
        self.neg_assign()
    }

    fn sqr(&self) -> Self {
        self.square_ref().complete()
    }

    fn sub(&self, rhs: &Self) -> Self {
        (self - rhs).complete()
    }

    fn random<R: rand_core::CryptoRng>(_rng: &mut R) -> Self {
        let mut rand = rug::rand::RandState::new();
        Integer::random_bits(4096, &mut rand).into()
    }

    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self * other).complete().modulo_ref(modulo).complete()
    }

    fn sub_mod(&self, other: &Self, modulo: &Self) -> Self {
        (self - other).complete().modulo_ref(modulo).complete()
    }

    fn neg(&self) -> Self {
        (-self).complete()
    }

    fn from(n: u64) -> Self {
        <Self as From<u64>>::from(n)
    }

    fn zero() -> Self {
        <Self as From<u64>>::from(0u64)
    }

    fn root(&self, n: u32) -> Self {
        self.root_ref(n).complete()
    }

    fn take_mod(&self, modulo: &Self) -> Self {
        self.modulo_ref(modulo).complete()
    }

    fn set_sign(&mut self, positive: bool) {
        match (self.is_positive(), positive) {
            (true, true) => (),
            (true, false) => self.oppose(),
            (false, true) => self.oppose(),
            (false, false) => (),
        }
    }

    fn less_than_abs(&self, rhs: &Self) -> bool {
        self.cmp_abs(rhs) == Ordering::Less
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized,
    {
        self.extended_gcd_ref(other).complete()
    }

    fn div_floor(&self, other: Self) -> Self {
        DivRounding::div_floor(self, other)
    }
}
