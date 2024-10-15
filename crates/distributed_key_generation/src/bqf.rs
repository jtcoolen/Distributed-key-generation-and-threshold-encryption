// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use std::cmp::Ordering::{Equal, Greater};
use std::fmt::Debug;

use serde::Serialize;

use crate::z::{self};

pub trait BinaryQuadraticForm<Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone + PartialEq,
{
    fn new(a: &Z, b: &Z, c: &Z) -> Self;

    fn new_with_discriminant(a: &Z, b: &Z, discriminant: &Z) -> Self;

    fn a(&self) -> Z;

    fn b(&self) -> Z;

    fn c(&self) -> Z;

    fn equals(&self, other: &Self) -> bool;

    fn rho(&self) -> Self;

    fn discriminant(&self) -> Z;

    fn identity(&self) -> Self;

    fn normalize(&self) -> Self;

    fn reduce(&self) -> Self;

    fn is_normal(&self) -> bool;

    fn is_reduced(&self) -> bool;

    fn compose(&self, other: &Self) -> Self
    where
        Self: Sized,
    {
        let mut g = self.b().add(&other.b());
        g.divide_by_2_exact();
        let w = self.a().gcd(&other.a()).gcd(&g);
        let mut h = other.b().sub(&self.b());
        h.divide_by_2_exact();
        let j = Clone::clone(&w);
        let s = self.a().divide_exact(&w);
        let t = other.a().divide_exact(&w);
        let u = g.divide_exact(&w);
        let st = s.mul(&t);
        let (mu, nu) = t.mul_mod(&u, &st).solve_congruence(
            &h.mul_mod(&u, &st).add_mod(&s.mul_mod(&self.c(), &st), &st),
            &st,
        );
        let (lambda, _) = t
            .mul_mod(&nu, &s)
            .solve_congruence(&h.sub_mod(&t.mul_mod(&mu, &s), &s), &s);
        let k = mu.add(&nu.mul(&lambda));
        let l = k.mul(&t).sub(&h).divide_exact(&s);
        let m = t
            .mul(&u)
            .mul(&k)
            .sub(&h.mul(&u))
            .sub(&self.c().mul(&s))
            .divide_exact(&st);
        let b = j.mul(&u).sub(&k.mul(&t)).sub(&l.mul(&s));
        let c = k.mul(&l).sub(&j.mul(&m));
        Self::new(&st, &b, &c).reduce()
    }

    // For squaring https://www.michaelstraka.com/posts/classgroups/
    // optimization if discriminant is negative of a prime
    fn double(&self) -> Self
    where
        Self: Sized,
    {
        let (mu, _) = self.b().solve_congruence(&self.c(), &self.a());
        let a = self.a().sqr();
        let b = self.b().sub(&Z::from(2).mul(&self.a()).mul(&mu));
        let rhs = self.b().mul(&mu).sub(&self.c()).divide_exact(&self.a());
        let c = mu.sqr().sub(&rhs);
        Self::new(&a, &b, &c).reduce()
    }

    fn pow(&self, exponent: &Z) -> Self
    where
        Self: Sized + Clone,
    {
        // naive double and add
        let mut res = self.identity();
        let n = exponent.bit_size();
        for i in (0..n).rev() {
            res = res.compose(&res);
            if exponent.get_bit(i) {
                res = res.compose(self);
            }
        }
        res
    }

    fn inverse(self) -> Self;

    fn to_bytes(&self) -> Vec<u8>;

    fn prime_form(discriminant: &Z, prime: &Z) -> Self
    where
        Self: Sized,
    {
        let a = prime.clone();
        let mut b = discriminant.sqrt_mod_prime(prime).unwrap();
        if b.is_odd() != discriminant.is_odd() {
            /* not same parity */
            b = prime.sub(&b);
        }
        Self::new_with_discriminant(&a, &b, discriminant).reduce()
    }

    #[cfg(feature = "gmp")]
    fn class_number_bound(discriminant: &Z) -> Z;

    fn prime_to(&self, l: &Z) -> Self;

    fn to_maximal_order(&self, l: &Z, delta_k: &Z) -> BQF<Z>;

    fn kernel_representative(&self, l: &Z, delta_k: &Z) -> Z;

    fn multiexp(x: &[BQF<Z>], e: &[Z]) -> BQF<Z>;
}

#[derive(Debug, Clone, Serialize)]
#[warn(dead_code)]
#[derive(PartialEq)]
pub struct BQF<Z>
where
    Z: crate::z::Z,
{
    a: Z,
    b: Z,
    c: Z,
}

impl<Z: z::Z + std::fmt::Debug + std::clone::Clone + std::cmp::PartialEq> BinaryQuadraticForm<Z>
    for BQF<Z>
{
    fn new(a: &Z, b: &Z, c: &Z) -> Self {
        BQF {
            a: Clone::clone(a),
            b: Clone::clone(b),
            c: Clone::clone(c),
        }
    }

    fn new_with_discriminant(a: &Z, b: &Z, discriminant: &Z) -> Self {
        let mut c = b.sqr().sub(discriminant).divide_exact(a);
        c.divide_by_4_exact();
        BQF {
            a: Clone::clone(a),
            b: Clone::clone(b),
            c,
        }
    }

    fn a(&self) -> Z {
        self.a.clone()
    }

    fn b(&self) -> Z {
        self.b.clone()
    }

    fn c(&self) -> Z {
        self.c.clone()
    }

    fn equals(&self, other: &Self) -> bool {
        self.a.eq(&other.a) && self.b.eq(&other.b) && self.c.eq(&other.c)
    }

    fn discriminant(&self) -> Z {
        self.b.sqr().sub(&Z::from(4).mul(&self.a.mul(&self.c)))
    }

    fn identity(&self) -> Self {
        let disc = self.discriminant();
        let b = if self.b.is_odd() {
            Z::from(1)
        } else {
            Z::from(0)
        };
        let mut c = Clone::clone(&b.sub(&disc));
        c.divide_by_4_exact();
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
        .reduce()
    }

    fn normalize(&self) -> Self {
        let z::EuclideanDivResult {
            mut quotient,
            mut remainder,
        } = self.b.euclidean_div_ceil(&self.a);
        if quotient.is_odd() {
            remainder = remainder.add(&self.a)
        };
        quotient.divide_by_2();

        let b = remainder.clone();

        remainder = self.b.add(&remainder);
        remainder.divide_by_2();
        let c = self.c.sub(&quotient.mul(&remainder));
        BQF {
            a: self.a.clone(),
            b,
            c,
        }
    }

    fn rho(&self) -> Self {
        BQF {
            a: self.c.clone(),
            b: self.b.clone().neg(),
            c: self.a.clone(),
        }
        .normalize()
    }

    fn reduce(&self) -> BQF<Z> {
        let mut n;
        if !self.is_normal() {
            n = self.normalize()
        } else {
            n = self.clone()
        };
        let mut cmp;
        while {
            cmp = n.a.compare(&n.c);
            cmp == Greater
        } {
            n = n.rho();
        }
        if cmp == Equal && !n.b.is_positive() {
            n.b.oppose();
        }
        n
    }

    fn is_normal(&self) -> bool {
        self.b.compare(&self.a).is_le() && self.b.compare(&self.a.neg()).is_gt()
    }

    fn is_reduced(&self) -> bool {
        self.a.compare(&self.c).is_le()
            && !(self.a.eq(&self.c) && self.b.compare(&Z::zero()).is_lt()) // !self.b.is_positive())
    }

    fn to_bytes(&self) -> Vec<u8> {
        // Convert each field to bytes
        let mut bytes = Vec::new();
        let a = self.c.to_bytes_be();
        bytes.extend_from_slice(&a.0);
        bytes.extend_from_slice(&[u8::from(a.1)]);
        let b = self.c.to_bytes_be();
        bytes.extend_from_slice(&b.0);
        bytes.extend_from_slice(&[u8::from(b.1)]);
        let c = self.c.to_bytes_be();
        bytes.extend_from_slice(&c.0);
        bytes.extend_from_slice(&[u8::from(c.1)]);
        bytes
    }

    #[cfg(feature = "gmp")]
    fn class_number_bound(discriminant: &Z) -> Z {
        use rug::{ops::DivFrom, Float, Integer};

        use crate::z::Z;
        let abs_d = discriminant.abs();
        let sqrt_d = abs_d.sqrt();
        let primebound = Integer::from_bytes_be(abs_d.ceil_abslog_square().to_bytes_be().0, true);
        let precision = discriminant.bit_size() as u32;

        let mut product = Float::with_val(precision, 1);
        let mut l: Integer = Integer::from_i64(2);

        let disc = discriminant.to_bytes_be();
        let disc = Integer::from_bytes_be(disc.0, disc.1);
        let mut tmp: Float;
        while l.compare(&primebound).is_le() {
            let kronecker = disc.kronecker(&l);
            tmp = Float::with_val(precision, l.clone());
            tmp.div_from(&Float::with_val(
                precision,
                l.clone() - Integer::from_i64(kronecker as i64),
            ));
            product *= tmp;
            l = l.next_prime();
        }

        let sqrt_d = sqrt_d.to_bytes_be();
        let sqrt_d = Float::with_val(precision, Integer::from_bytes_be(sqrt_d.0, sqrt_d.1));
        let factor: Float = sqrt_d * Float::with_val(precision, 21);
        let bound: Float = (factor * product / Float::with_val(precision, 88)).floor();
        let (b, pos) = bound.to_integer().unwrap().to_bytes_be();
        Z::from_bytes_be(b, pos)
    }

    fn prime_to(&self, l: &Z) -> Self {
        let mut g = self.a().gcd(l);
        let mut a = self.a();
        let mut b = self.b();
        let mut c = self.c();
        if g.compare(&Z::from(1u64)) == Greater {
            g = self.c().gcd(l);
            if g.compare(&Z::from(1u64)) == Greater {
                // Transform f into (a + b + c, -b - 2a, a)
                c = c.add(&self.a());
                c = c.add(&self.b());
                b = b.add(&self.a());
                b = b.add(&self.a());
                b = b.neg();
                std::mem::swap(&mut a, &mut c);
            } else {
                // c is coprime to l: transform f into (c, -b, a)
                std::mem::swap(&mut a, &mut c);
                b = b.neg();
            }
        }
        // else do nothing if a is already coprime to l
        BQF::new(&a, &b, &c)
    }

    fn to_maximal_order(&self, l: &Z, delta_k: &Z) -> BQF<Z> {
        let new_self = self.prime_to(l);

        let (_, g0, g1) = l.extended_gcd(&new_self.a());

        let mut b = new_self.b().mul(&g0);
        b = b.add(&new_self.a().mul(&g1));

        BQF::new_with_discriminant(&new_self.a(), &b, delta_k)
    }

    fn kernel_representative(&self, l: &Z, delta_k: &Z) -> Z {
        let mut tmp0 = Z::default();
        let mut g0 = Z::from(1u64);
        let mut g1 = Z::from(0u64);

        let ft = self.clone();
        ft.to_maximal_order(l, delta_k);

        let mut ft = ft.normalize();
        while ft.a().abs().compare(&ft.c().abs()) == Greater {
            tmp0 = g1.mul(delta_k);
            g1 = g1.mul(&ft.b());
            g1 = g1.add(&g0);
            g0 = g0.mul(&ft.b());
            g0 = g0.add(&tmp0);

            ft = ft.rho();
        }

        if !ft.a().eq(&Z::from(1u64)) || !ft.b().eq(&Z::from(1u64)) {
            panic!("The form is not in the kernel")
        }

        let tmp1 = g0.gcd(&g1);
        g0 = g0.divide_exact(&tmp1);
        g1 = g1.divide_exact(&tmp1);

        tmp0 = g0.invert_mod(l).unwrap();
        g1 = g1.neg();
        tmp0 = tmp0.mul_mod(&g1, l);

        tmp0
    }

    fn multiexp(x: &[BQF<Z>], e: &[Z]) -> BQF<Z> {
        // Check for size inconsistencies or empty vectors
        if x.len() != e.len() || x.is_empty() {
            panic!("invalid size");
        }

        let blank = BQF::new(&Z::from(1), &Z::from(1), &Z::from(1));
        let n = x.len();
        let mut s = blank.clone();
        let mut r = blank.clone();
        let mut b = vec![blank.clone(); 16];
        let mut p = blank.clone(); // Initialize P

        // Find the largest element in 'e'
        let mut mt = e[0].clone();
        for i in 1..n {
            if e[i].compare(&mt) == Greater {
                mt = e[i].clone();
            }
        }

        let nb = (mt.bit_size() + 3) / 4;

        // Implementing Pippenger's algorithm
        for i in (0..nb).rev() {
            // Reset the array B for each iteration
            for j in 0..16 {
                b[j] = blank.clone();
            }

            for j in 0..n {
                mt = e[j].clone();
                Z::divide_pow2(&mut mt, (4 * i) as usize);
                let mut res = mt.clone();
                res.rem_trunc_inplace(&Z::from(16));

                let vv = res; // Assuming appropriate conversion
                let mut array = [0u8; std::mem::size_of::<usize>()];
                let bb = vv.to_bytes_be();
                if !bb.0.is_empty() {
                    array[0] = bb.0[bb.0.len() - 1];
                }

                let k = usize::from_ne_bytes(array);

                if b[k] == blank {
                    b[k] = x[j].clone();
                } else {
                    b[k] = BQF::compose(&b[k], &x[j]);
                }
            }

            r = blank.clone();
            s = blank.clone();

            for j in (1..=15).rev() {
                if b[j] != blank {
                    if r == blank {
                        r = b[j].clone();
                    } else {
                        r = r.compose(&b[j]);
                    }
                }

                if r != blank {
                    if s == blank {
                        s = r.clone();
                    } else {
                        s = s.compose(&r);
                    }
                }
            }

            if p != blank {
                for _ in 0..4 {
                    p = p.compose(&p);
                }
            }

            if s != blank {
                if p == blank {
                    p = s.clone();
                } else {
                    p = p.compose(&s);
                }
            }
        }

        p.reduce()
    }
}
