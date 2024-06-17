//use rand_core::CryptoRng;

use crate::z;

trait BinaryQuadraticForm<Z>
where
    Z: crate::z::Z,
{
    fn new(a: Z, b: Z, c: Z) -> Self;

    //fn random<R: CryptoRng>(rng: &mut R) -> Self;

    fn discriminant(self) -> Z;

    fn identity(self) -> Self;

    fn normalize(self) -> Self;

    fn reduce(self) -> Self;

    fn compose(self, other: Self) -> Self;

    fn double(self) -> Self;

    fn pow(self, exponent: Z) -> Self;

    fn inverse(self) -> Self;
}

#[warn(dead_code)]
struct BQF<Z>
where
    Z: crate::z::Z,
{
    a: Z,
    b: Z,
    c: Z,
}

impl<Z: z::Z> BQF<Z> {
    fn rho(self) -> Self {
        BQF {
            a: self.c,
            b: self.b.neg(),
            c: self.a,
        }
        .normalize()
    }
}

impl<Z: z::Z> BinaryQuadraticForm<Z> for BQF<Z> {
    fn new(a: Z, b: Z, c: Z) -> Self {
        BQF { a, b, c }
    }

    fn discriminant(self) -> Z {
        self.b.sqr().sub(Z::from(4).mul(self.a.mul(self.c)))
    }

    fn identity(self) -> Self {
        let disc = self.discriminant();
        let b = if self.b.is_odd() {
            Z::from(1)
        } else {
            Z::from(0)
        };
        let mut c = self.b.sub(disc);
        c.divide_by_4();
        BQF {
            a: Z::from(1),
            b,
            c,
        }
    }

    fn inverse(self) -> Self {
        BQF {
            a: self.a,
            b: self.b.neg(),
            c: self.c,
        }
    }

    // TODO check/test correctness and benchmark
    // TODO inplace operations?
    // TODO check for overflow
    fn normalize(self) -> Self {
        let z::EuclideanDivResult {
            mut quotient,
            remainder,
        } = self.b.euclidean_div_ceil(self.a);
        let remainder = if quotient.is_odd() {
            remainder.add(self.a)
        } else {
            remainder
        };
        quotient.divide_by_2();
        let b = remainder;
        let mut remainder = self.b.add(remainder);
        remainder.divide_by_2();
        let c = self.c.sub(quotient.mul(remainder));
        BQF { a: self.a, b, c }
    }

    fn reduce(self) -> Self {
        let mut n = self.normalize();
        while n.a.less_than_abs(n.c) {
            n = n.rho();
        }
        if n.a.eq_abs(n.c) && !n.b.is_positive() {
            n.b.oppose()
        }
        n
    }

    // TODO (harder): implement compose, double, pow
    //https://gite.lirmm.fr/crypto/bicycl/-/blob/master/src/bicycl/qfi.inl?ref_type=heads#L621

    // TODO (less of a priority): compute class group number
}

// TODO instantiate with BQF<Bignum4096> for a discriminant of 1827 bits
// security level of 128 bits
// since |b|<=|a|<= sqrt(|D|/3) and |c|<=|D| the coefficients of the form fit in a 4096-bit big num
// TODO check that the operations do not overflow!
