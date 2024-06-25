use rand_core::CryptoRng;
// TODO: Rug is using floats! Replace by "raw" GMP bindings
use rug::ops::NegAssign;
use rug::Complete;
use rug::Integer;

use crate::z::EuclideanDivResult;
use crate::z::ExtendedGCDResult;
use crate::z::Z;

/// Signed integer with range (-2^4096, 2^4096)
/// Does not check for overflows!
#[derive(Debug, Clone)]
pub struct Bignum4096 {
    pub positive: bool,
    // we take 128 = 64 * 2 to hold the results of sqr and mul methods
    pub limbs: [u64; 128],
}

impl PartialEq for Bignum4096 {
    fn eq(&self, rhs: &Self) -> bool {
        self.eq_abs(rhs) && self.positive == rhs.positive
    }
}

fn from_bignum4096(n: &Bignum4096) -> Integer {
    let mut x = Integer::from_digits(&n.limbs, rug::integer::Order::Lsf);
    if !n.positive {
        x.neg_assign()
    }
    x
}

fn to_bignum4096(n: &Integer) -> Bignum4096 {
    let digits = n.to_digits::<u64>(rug::integer::Order::Lsf);
    let mut limbs = [0u64; 128];
    limbs[0..digits.len()].copy_from_slice(&digits);
    Bignum4096 {
        positive: n.is_positive(),
        limbs,
    }
}

// Constant-time functions
impl Z for Bignum4096 {
    fn default() -> Self {
        Bignum4096 {
            positive: true,
            limbs: [0; 128],
        }
    }
    fn zero() -> Self {
        Bignum4096 {
            positive: true,
            limbs: [0; 128],
        }
    }

    fn clone(&self) -> Self {
        Clone::clone(&self)
    }

    // returns a positive integer sampled uniformly at random
    fn random<R: CryptoRng>(rng: &mut R) -> Self {
        let mut limbs = [0u64; 128];
        let mut dest = [0u8; 8 * 64];
        rng.fill_bytes(&mut dest);
        for (i, chunk) in dest.chunks_exact(8).enumerate() {
            limbs[i] = u64::from_le_bytes(chunk.try_into().expect("Chunk size is not 8"));
        }
        Bignum4096 {
            positive: true,
            limbs,
        }
    }

    fn eq_abs(&self, rhs: &Self) -> bool {
        let eq: u64 = unsafe {
            hacl_sys::Hacl_Bignum4096_eq_mask(
                self.limbs.as_ptr() as *mut _,
                rhs.limbs.as_ptr() as *mut _,
            )
        };
        eq != 0
    }

    fn eq(&self, rhs: &Self) -> bool {
        self.eq_abs(rhs) && self.positive == rhs.positive
    }

    fn from(n: u64) -> Self {
        let mut res = Self::zero();
        res.limbs[0] = n;
        res
    }

    // TODO: is the inequality strict or not?
    fn less_than_abs(&self, rhs: &Self) -> bool {
        let lt: u64 = unsafe {
            hacl_sys::Hacl_Bignum4096_lt_mask(
                self.limbs.as_ptr() as *mut _,
                rhs.limbs.as_ptr() as *mut _,
            )
        };
        lt != 0
    }

    fn less_than(&self, rhs: &Self) -> bool {
        let both_positive = self.positive && rhs.positive;
        let both_negative = !self.positive && !rhs.positive;
        let self_negative_rhs_positive = !self.positive && rhs.positive;
        let less_than_abs = self.less_than_abs(&rhs);

        let result_if_both_positive = both_positive && less_than_abs;
        let result_if_both_negative = both_negative && !less_than_abs;
        let result_if_self_negative_rhs_positive = self_negative_rhs_positive;

        result_if_both_positive | result_if_both_negative | result_if_self_negative_rhs_positive
    }

    fn add(&self, rhs: &Self) -> Self {
        let mut res = Self::zero();
        // TODO use subtle crate for booleans, or use u32
        match (self.positive, rhs.positive, self.less_than_abs(&rhs)) {
            (true, true, _) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_add(
                        self.limbs.as_ptr() as *mut _,
                        rhs.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
                res.positive = true;
            }
            (false, false, _) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_add(
                        self.limbs.as_ptr() as *mut _,
                        rhs.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
                res.positive = false;
            }
            (true, false, true) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_sub(
                        rhs.limbs.as_ptr() as *mut _,
                        self.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
                res.positive = false;
            }
            (true, false, false) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_sub(
                        self.limbs.as_ptr() as *mut _,
                        rhs.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
                res.positive = true;
            }
            (false, true, true) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_sub(
                        rhs.limbs.as_ptr() as *mut _,
                        self.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
                res.positive = true;
            }
            (false, true, false) => {
                unsafe {
                    hacl_sys::Hacl_Bignum4096_sub(
                        self.limbs.as_ptr() as *mut _,
                        rhs.limbs.as_ptr() as *mut _,
                        res.limbs.as_mut_ptr(),
                    );
                }
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
        unsafe {
            hacl_sys::Hacl_Bignum4096_mul(
                self.limbs.as_ptr() as *mut _,
                rhs.limbs.as_ptr() as *mut _,
                res.limbs.as_mut_ptr(),
            );
        }
        res.positive = self.positive == rhs.positive;
        res
    }

    fn sqr(&self) -> Self {
        let mut res = Self::zero();
        unsafe {
            hacl_sys::Hacl_Bignum4096_sqr(self.limbs.as_ptr() as *mut _, res.limbs.as_mut_ptr());
        }
        res
    }

    fn neg(&self) -> Self {
        Bignum4096 {
            positive: !self.positive,
            limbs: self.limbs.clone(),
        }
    }

    fn divide_by_2(&mut self) {
        let mut carry: u64 = 0;
        let shift = 63 as u8;
        for limb in self.limbs.iter_mut().rev() {
            let new_carry = *limb & 0b1;
            *limb = (*limb >> 1) | (carry << shift);
            carry = new_carry;
        }
    }

    fn divide_by_4(&mut self) {
        let mut carry: u64 = 0;
        let shift = 62;
        for limb in self.limbs.iter_mut().rev() {
            let new_carry = *limb & 0b11;
            *limb = (*limb >> 2) | (carry << shift);
            carry = new_carry;
        }
    }

    fn is_odd(&self) -> bool {
        self.limbs[0] & 0b1 == 1
    }

    fn oppose(&mut self) {
        self.positive = !self.positive
    }

    fn is_positive(&self) -> bool {
        self.positive
    }

    // make call to mpn_sec_div_qr from gmp
    // use https://docs.rs/gmp-mpfr-sys/latest/i686-unknown-linux-gnu/gmp_mpfr_sys/gmp/fn.mpn_sec_div_qr.html#
    fn euclidean_div_ceil(&self, other: &Self) -> EuclideanDivResult<Self>
    where
        Self: Sized,
    {
        let n = from_bignum4096(self);
        let d = from_bignum4096(other);

        let (quotient, rem) = n.div_rem_ceil(d);
        // TODO: use gmp_mpfr_sys::gmp::mpn_sec

        let quotient = to_bignum4096(&quotient);
        let remainder = to_bignum4096(&rem);
        EuclideanDivResult {
            quotient,
            remainder,
        }
    }

    fn root(&self, n: u32) -> Self {
        let mut x = from_bignum4096(self);
        Integer::root_mut(&mut x, n);
        to_bignum4096(&x)
    }

    fn gcd(&self, other: &Self) -> Self {
        let a = from_bignum4096(self);
        let b = from_bignum4096(other);
        to_bignum4096(&a.gcd(&b))
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized,
    {
        let a = from_bignum4096(self);
        let b = from_bignum4096(other);
        let (g, u, v) = a.extended_gcd_ref(&b).complete();
        (to_bignum4096(&g), to_bignum4096(&u), to_bignum4096(&v))
    }

    /// Implements the partial extended GCD from https://eprint.iacr.org/2022/1466.pdf
    /// Basically Lehmer's variant that's truncated until the remainders get lower than the upper bound
    /// https://gite.lirmm.fr/crypto/bicycl/-/blob/master/src/bicycl/gmp_extras.inl?ref_type=heads#L1232
    /// see https://eprint.iacr.org/2021/1292.pdf (reported to be constant-time)
    /// https://perso.ens-lyon.fr/damien.stehle/downloads/recbinary.pdf (less efficient, just for curiosity)
    /// sounds nice A Double-Digit Lehmer-Euclid Algorithm for Finding the GCD of Long Integers TUDOR JEBELEAN
    /// PROBABLY GO WITH GMP LEHMER'S IMPLEM
    /// For now let's go with a naive extended GCD with euclidean divisions! let's start with something working first
    /// TODO then abstract method using primitives from Z only
    fn partial_extended_gcd(
        &self,
        other: &Self,
        bezout_coefficients_upper_bound: &Self,
    ) -> crate::z::ExtendedGCDResult<Self>
    where
        Self: Sized,
    {
        let upper_bound = from_bignum4096(bezout_coefficients_upper_bound);
        let mut r = from_bignum4096(self); // a
        let mut r_next = from_bignum4096(other); // b
        // intermediate Bézout coefficients for a
        let mut u = Integer::ONE.clone();
        let mut u_next = Integer::ZERO.clone();
        // intermediate Bézout coefficients for b
        let mut v = Integer::ZERO.clone();
        let mut v_next = Integer::ONE.clone();

        loop {
            let (q, rem) = r.div_rem_ref(&r_next).complete();
            r = rem;
            u -= &q * &u_next;
            v -= &q * &v_next;
            std::mem::swap(&mut r, &mut r_next);
            std::mem::swap(&mut u, &mut u_next);
            std::mem::swap(&mut v, &mut v_next);
            if v <= upper_bound && u <= upper_bound {
                break;
            }
        }
        ExtendedGCDResult {
            bezout_coeff_1: to_bignum4096(&u),
            bezout_coeff_2: to_bignum4096(&v),
        }
    }

    fn divide_exact(&self, other: &Self) -> Self {
        let a = from_bignum4096(self);
        let b = from_bignum4096(other);
        to_bignum4096(&a.div_exact(&b))
    }

    fn divides(&self, other: &Self) -> bool {
        let a = from_bignum4096(self);
        let b = from_bignum4096(other);
        b.is_divisible(&a)
    }

    // TODO: recheck if we can use hacl star implem (yes we can)
    fn add_mod(&self, other: &Self, modulo: &Self) -> Self {
        let mut res = Self::zero();
        unsafe {
            hacl_sys::Hacl_Bignum4096_add_mod(
                modulo.limbs.as_ptr() as *mut _,
                self.limbs.as_ptr() as *mut _,
                other.limbs.as_ptr() as *mut _,
                res.limbs.as_mut_ptr(),
            );
        }
        res
    }

    fn take_mod(&self, modulo: &Self) -> Self {
        let mut res = Self::zero();
        unsafe {
            hacl_sys::Hacl_Bignum4096_mod(
                modulo.limbs.as_ptr() as *mut _,
                self.limbs.as_ptr() as *mut _,
                res.limbs.as_mut_ptr(),
            );
        }
        res
    }

    fn mul_mod(&self, other: &Self, modulo: &Self) -> Self {
        self.mul(other).take_mod(modulo)
    }

    fn set_sign(&mut self, positive: bool) {
        self.positive = positive;
    }
}
