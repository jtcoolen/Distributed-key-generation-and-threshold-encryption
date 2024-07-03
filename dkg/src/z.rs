use rand_core::CryptoRng;

pub struct EuclideanDivResult<Z> {
    pub quotient: Z,
    pub remainder: Z,
}

pub struct ExtendedGCDResult<Z> {
    pub bezout_coeff_1: Z,
    pub bezout_coeff_2: Z,
}

pub trait Z {
    fn zero() -> Self;

    fn default() -> Self;

    fn from(n: u64) -> Self;

    fn from_bytes(b: Vec<u8>) -> Self;

    fn to_bytes(&self) -> Vec<u8>;

    fn random<R: CryptoRng>(rng: &mut R) -> Self;

    fn sample_range<R: CryptoRng>(rng: &mut R, lower: &Self, upper: &Self) -> Self
    where
        Self: Sized,
    {
        let rd = Self::random(rng);
        rd.take_mod(&(upper.sub(lower).add(&Self::from(1))))
            .add(lower)
    }

    fn eq_abs(&self, rhs: &Self) -> bool;

    fn eq(&self, rhs: &Self) -> bool;

    fn less_than_abs(&self, rhs: &Self) -> bool;

    fn less_than(&self, rhs: &Self) -> bool;

    fn add(&self, rhs: &Self) -> Self;

    fn sub(&self, rhs: &Self) -> Self;

    fn mul(&self, rhs: &Self) -> Self;

    fn sqr(&self) -> Self;

    fn neg(&self) -> Self;

    fn divide_by_2(&mut self);

    fn divide_by_4(&mut self);

    fn is_odd(&self) -> bool;

    // rounds quotient towards +infinity
    // the remainder gets the opposite sign of the denominator
    fn euclidean_div_ceil(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized;

    fn oppose(&mut self);

    fn is_positive(&self) -> bool;

    fn root(&self, n: u32) -> Self;

    fn gcd(&self, other: &Self) -> Self;

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized;

    /*fn partial_extended_gcd(
        &self,
        other: &Self,
        bezout_coefficients_upper_bound: &Self,
    ) -> ExtendedGCDResult<Self>
    where
        Self: Sized;*/

    fn divide_exact(&self, other: &Self) -> Self;

    fn divides(&self, other: &Self) -> bool;

    fn add_mod(&self, other: &Self, modulo: &Self) -> Self;

    fn sub_mod(&self, other: &Self, modulo: &Self) -> Self;

    fn take_mod(&self, modulo: &Self) -> Self;

    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self;

    fn set_sign(&mut self, positive: bool);

    // Solves a congruence self * x = other [modulo] for x
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

    fn div_floor(&self, other: Self) -> Self;

    fn bit_size(&self) -> u64;

    fn get_bit(&self, index: u64) -> bool;
}
