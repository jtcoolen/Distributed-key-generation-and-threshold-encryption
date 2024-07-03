//use rand_core::CryptoRng;

use crate::z::{self};

pub trait BinaryQuadraticForm<Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
{
    fn new(a: &Z, b: &Z, c: &Z) -> Self;

    fn a(&self) -> Z;
    fn b(&self) -> Z;
    fn c(&self) -> Z;

    //fn random<R: CryptoRng>(rng: &mut R) -> Self;

    fn equals(&self, other: &Self) -> bool;

    fn discriminant(&self) -> Z;

    fn identity(&self) -> Self;

    fn normalize(&self) -> Self;

    fn reduce(&self) -> Self;

    /// https://www.ams.org/journals/mcom/2003-72-244/S0025-5718-03-01518-7/S0025-5718-03-01518-7.pdf
    /// Quadratic form (u,v,w)=uX^2+vXY+wY^2
    /// How to compose two forms (u1,v1,w1) and (u2,v2,w2) of same discriminant D < 0?
    /// (so in the same class group)
    ///
    /// By computing the 2x4 matrix
    /// M = [Ax Bx Cx Dx; Ay By Cy Dy]
    /// with the following constraints:
    /// - Ax = [gcd(u1,u2,(v1+v2)/2) =: G]
    /// - Ay = 0
    /// - By = u1/Ax
    /// - Cy = u2/Ax
    /// - Dy = [(v1+v2)/2 =: s] / Ax
    /// - [m := -(v1-v2)/2] = Bx*Cy-By*Cx
    /// - w1 = Cx*Dy-Cy*Dx
    /// - w2 = Bx*Dy-By*Dx
    ///
    /// The composition of the two forms is then
    /// (u3,v3,w3) = (By*Cy-Ay*Dy, (Ax*Dy+Ay*Dx)-(Bx*Cy+By*Cx), Bx*Cx-Ax*Dx)
    ///
    /// How to compute M concretely?
    /// Using Atkins' refinement
    ///
    /// 1. Compute G, Bx and By
    /// 2. Extended GCD on Bx and By: Bx*bx+By*by=gcd(Bx,By) until bx and by
    /// both are at most ceil(|D|^1/4)
    /// 3. this gives us the first two columns of M: [G; 0] and [bx; by]
    /// 4. then the remaining coefficients are given by
    /// If bx != 0:
    /// - cx = (bx*u2-m*ax)/u1
    /// - cy = (by*cx+m)/bx
    /// - dx = (bx*s-w2*ax)/u1
    /// - dy = (dx*ay+s)/ax
    /// If bx = 0 then [w1 0 -s u2] is orthogonal to both rows of M:
    /// that is the following equalities on the inner products
    /// ax*w1-cx*s+dx*u2 = 0
    /// ay*w1-cy*s+dy*u2 = 0
    /// give us
    /// dx = (cx*s-ax*w1)/u2
    /// cy = (dy*u2+ay*w1)/s
    /// "Gives a form very close to reduced"
    /// TODO: it would be interesting to see how to make this constant-time
    fn compose(&self, other: &Self) -> Self
    where
        Self: Sized,
    {
        let mut g = self.b().add(&other.b());
        g.divide_by_2();
        let w = self.a().gcd(&other.a()).gcd(&g);
        let mut h = other.b().sub(&self.b());
        h.divide_by_2();
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
        Self::new(&st, &b, &c).reduce().reduce().reduce()
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
        Self::new(&a, &b, &c)
    }

    // We can use double and add with Shamir's window trick (q-ary decomposition)
    // See https://gite.lirmm.fr/crypto/bicycl/-/blob/master/src/bicycl/qfi.inl?ref_type=heads#L1162
    // https://github.com/jtcoolen/GLV_arkworks/blob/main/src/lib.rs#L84
    // https://github.com/jtcoolen/asymmetric_crypto/blob/master/elliptic_curves/ec_elgamal_codage_decodage.gp#L100
    fn pow(&self, exponent: &Z) -> Self
    where
        Self: Sized + Clone,
    {
        // naive double and add
        let n = exponent.bit_size();
        let mut res = self.clone();
        for i in (0..=n).rev() {
            res = res.double();
            if exponent.get_bit(i) {
                res = self.compose(&res).reduce();
            }
        }
        res
    }

    fn inverse(self) -> Self;

    fn to_bytes(&self) -> Vec<u8>;
}

#[derive(Debug, Clone)]
#[warn(dead_code)]
pub struct BQF<Z>
where
    Z: crate::z::Z,
{
    a: Z,
    b: Z,
    c: Z,
}

impl<Z: z::Z + std::fmt::Debug + std::clone::Clone> BQF<Z> {
    fn rho(self) -> Self {
        BQF {
            a: self.c,
            b: self.b.neg(),
            c: self.a,
        }
        .normalize()
    }
}

impl<Z: z::Z + std::fmt::Debug + std::clone::Clone> BinaryQuadraticForm<Z> for BQF<Z> {
    fn new(a: &Z, b: &Z, c: &Z) -> Self {
        BQF {
            a: Clone::clone(&a),
            b: Clone::clone(&b),
            c: Clone::clone(&c),
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
    fn normalize(&self) -> Self {
        let z::EuclideanDivResult {
            mut quotient,
            remainder,
        } = self.b.euclidean_div_ceil(&self.a);
        let remainder = if quotient.is_odd() {
            remainder.add(&self.a)
        } else {
            remainder
        };
        quotient.divide_by_2();
        let b = Clone::clone(&remainder);
        let mut remainder = self.b.add(&remainder);
        remainder.divide_by_2();
        let c = self.c.sub(&quotient.mul(&remainder));
        BQF {
            a: self.a.clone(),
            b,
            c,
        }
    }

    /*fn reduce(&self) -> Self {
        let mut n = self.normalize();
        while n.a.less_than_abs(&n.c) {
            n = n.rho();
        }
        if n.a.eq_abs(&n.c) && !n.b.is_positive() {
            n.b.oppose()
        }
        n
    }*/
    fn reduce(&self) -> BQF<Z> {
        let mut h: BQF<Z>;
        let mut h_new = self.clone();
        if !is_normal2(&self) {
            h_new = self.normalize();
        }
        h = h_new;
        while !is_reduced2(&h) {
            let h_new = rho2(&h);
            h = h_new;
        }
        h
    }

    fn to_bytes(&self) -> Vec<u8> {
        // Convert each field to bytes
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.a.to_bytes());
        bytes.extend_from_slice(&self.b.to_bytes());
        bytes.extend_from_slice(&self.c.to_bytes());
        bytes
    }

    // TODO (harder): implement compose, double, pow
    // reverse-engineer their implementation
    //https://gite.lirmm.fr/crypto/bicycl/-/blob/master/src/bicycl/qfi.inl?ref_type=heads#L621

    // TODO (less of a priority): compute class group number
}

fn normalize2<Z: z::Z + Clone>(x: BQF<Z>) -> BQF<Z> {
    // assume delta<0 and a>0
    let a_sub_b = x.a.sub(&x.b);
    let s_f = a_sub_b.div_floor(Z::from(2).mul(&x.a));

    BQF {
        a: x.a.clone(),
        b: x.b.add(&Z::from(2).mul(&s_f).mul(&x.a)),
        c: x.a.mul(&s_f.sqr()).add(&x.b.mul(&s_f)).add(&x.c),
    }
}

fn rho2<Z: z::Z + Clone>(x: &BQF<Z>) -> BQF<Z> {
    let qf_new = BQF {
        a: x.c.clone(),
        b: x.b.clone().neg(),
        c: x.a.clone(),
    };

    normalize2(qf_new)
}

fn is_normal2<Z: z::Z>(x: &BQF<Z>) -> bool {
    x.b.less_than(&x.a) && !(x.b.less_than(&x.a.neg()))
}

fn is_reduced2<Z: z::Z>(x: &BQF<Z>) -> bool {
    is_normal2(x) && x.a.less_than(&x.c) && !(x.a.eq(&x.c) && x.b.less_than(&Z::zero()))
}

fn reduce2<Z: z::Z + Clone>(x: BQF<Z>) -> BQF<Z> {
    let mut h: BQF<Z>;
    let mut h_new = x.clone();
    if !is_normal2(&x) {
        h_new = normalize2(x);
    }
    h = h_new;
    while !is_reduced2(&h) {
        let h_new = rho2(&h);
        h = h_new;
    }
    h
}

// TODO instantiate with BQF<Bignum4096> for a discriminant of 1827 bits
// security level of 128 bits
// since |b|<=|a|<= sqrt(|D|/3) and |c|<=|D| the coefficients of the form fit in a 4096-bit big num
// TODO check that the operations do not overflow!

#[cfg(test)]
mod tests {

    use bicycl::cpp_std::VectorOfUchar;
    use proptest::prelude::*;
    use proptest_derive::Arbitrary;

    use bicycl::b_i_c_y_c_l::{ClassGroup, Mpz, QFI};
    use bicycl::cpp_core::{self, CppBox, MutRef, Ref};
    use bicycl::{b_i_c_y_c_l, cpp_std::String};
    use rug::ops::NegAssign;
    use rug::Integer;
    use z::Z;

    use std::os::raw::c_char;

    use crate::bignum4096::Bignum4096;

    use super::*;

    pub fn convert(v: Vec<u8>) -> [u64; 128] {
        assert!(v.len() <= 512);
        let mut array: [u8; 512] = [0; 512];
        array[..v.len()].copy_from_slice(&v);
        let cast: [u64; 64] = bytemuck::cast(array);
        let mut res = [0u64; 128];
        res[0..64].copy_from_slice(&cast);
        res
    }

    fn mpz_to_bignum(n: &mut CppBox<Mpz>) -> Bignum4096 {
        let limbs = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *n)) };

        let mut limbs = limbs.clone();

        limbs.reverse();
        Bignum4096 {
            positive: unsafe { Mpz::sgn(&n) == 1 }.clone(),
            limbs: convert(limbs).clone(),
        }
    }

    fn mpz_to_bignum1(n: &Mpz) -> Integer {
        let x = mpz_to_bignum2(n);
        let limbs = x.limbs;

        //limbs.reverse();
        let mut x = Integer::from_digits(&limbs, rug::integer::Order::Lsf);
        let sg = unsafe { Mpz::sgn(&n) == 1 };
        if !sg {
            x.neg_assign()
        }
        x
    }

    fn mpz_to_bignum2(n: &Mpz) -> Bignum4096 {
        let mut a_bytes = unsafe { VectorOfUchar::new() };

        let mutref_a_bytes: cpp_core::MutRef<VectorOfUchar> =
            unsafe { cpp_core::MutRef::from_raw_ref(&mut a_bytes) };

        unsafe { n.mpz_to_vector(mutref_a_bytes) };
        let limbs = unsafe { bicycl::cpp_vec_to_rust(&mutref_a_bytes) };
        let mut limbs = limbs[1..].to_vec();
        println!("LIMBS={:?}", limbs);
        limbs.reverse();
        Bignum4096 {
            positive: unsafe { Mpz::sgn(&n) == 1 },
            limbs: convert(limbs),
        }
    }

    // TODO write proptest version of this test
    // TODO also write benchmarks
    #[test]
    fn test_discriminant() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "219211015245339659606923489058910718059300326777750522511052678189518994793540424066835449513550321156725825999213223313349543556000887051910142135883849284272371752572695727069191024343809964931675197912960671210283614957513396516801587047437150725164417736257899198555560707024341197521616471833363582196266114484743";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "-81545512674410670486936176483359392564511873763104542618989110464161210362234569734599379968174823311559519130742888365338004517383649260007974211906536061358807858053200069189129489104433807401679673870588655431254253093929351305686001886690914729246160866675068623497789888497184276822188325873470578904605036981681";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "397911913619280235397607492118765996458330730701094047301004809429554024195731134066384833319157074257034700248114739847099306289375744183458513982780336763663761492062176352821924273594044341543655974705344319285820641420038042279276523776451336649974742547552579516047919797400613300929758921320569816271189292178311";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let qfi = unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum1(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum1(&mut a),
            &mpz_to_bignum1(&mut b),
            &mpz_to_bignum1(&mut c),
        );

        println!("qfi2={:?}", qfi2);
        println!("qfi2={:?}, disc={:?}", qfi2, qfi2.discriminant());
        assert!(qfi2.discriminant() == mpz_to_bignum1(&mut disc))
    }

    #[test]
    fn test_identity() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "219211015245339659606923489058910718059300326777750522511052678189518994793540424066835449513550321156725825999213223313349543556000887051910142135883849284272371752572695727069191024343809964931675197912960671210283614957513396516801587047437150725164417736257899198555560707024341197521616471833363582196266114484743";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "-81545512674410670486936176483359392564511873763104542618989110464161210362234569734599379968174823311559519130742888365338004517383649260007974211906536061358807858053200069189129489104433807401679673870588655431254253093929351305686001886690914729246160866675068623497789888497184276822188325873470578904605036981681";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "397911913619280235397607492118765996458330730701094047301004809429554024195731134066384833319157074257034700248114739847099306289375744183458513982780336763663761492062176352821924273594044341543655974705344319285820641420038042279276523776451336649974742547552579516047919797400613300929758921320569816271189292178311";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let qfi: CppBox<b_i_c_y_c_l::QFI> =
            unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum1(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum1(&mut a),
            &mpz_to_bignum1(&mut b),
            &mpz_to_bignum1(&mut c),
        );

        let cl = unsafe { b_i_c_y_c_l::ClassGroup::new(&disc) };

        let one = unsafe { cl.one() };
        let a = unsafe { one.a() };
        let b = unsafe { one.b() };
        let c = unsafe { one.c() };

        println!(
            "one={:?},{:?},{:?}",
            mpz_to_bignum2(&a),
            mpz_to_bignum2(&b),
            mpz_to_bignum2(&c),
        );
        assert!(mpz_to_bignum1(&a) == qfi2.identity().a);
        assert!(mpz_to_bignum1(&b) == qfi2.identity().b);
        assert!(mpz_to_bignum1(&c) == qfi2.identity().c);

        println!("qfi2={:?}", qfi2);
        println!("qfi2={:?}, id={:?}", qfi2, qfi2.identity());

        println!(
            "disc qfi2={:?}, disc={:?}",
            qfi2.discriminant(),
            mpz_to_bignum1(&mut disc)
        );
        assert!(qfi2.discriminant() == mpz_to_bignum1(&mut disc))
    }

    #[test]
    fn test_normalize() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "219211015245339659606923489058910718059300326777750522511052678189518994793540424066835449513550321156725825999213223313349543556000887051910142135883849284272371752572695727069191024343809964931675197912960671210283614957513396516801587047437150725164417736257899198555560707024341197521616471833363582196266114484743";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "-81545512674410670486936176483359392564511873763104542618989110464161210362234569734599379968174823311559519130742888365338004517383649260007974211906536061358807858053200069189129489104433807401679673870588655431254253093929351305686001886690914729246160866675068623497789888497184276822188325873470578904605036981681";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "397911913619280235397607492118765996458330730701094047301004809429554024195731134066384833319157074257034700248114739847099306289375744183458513982780336763663761492062176352821924273594044341543655974705344319285820641420038042279276523776451336649974742547552579516047919797400613300929758921320569816271189292178311";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let mut qfi = unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum(&mut a),
            &mpz_to_bignum(&mut b),
            &mpz_to_bignum(&mut c),
        );
        unsafe { qfi.normalize_0a() };
        unsafe { qfi.normalize_0a() };

        let mut aaa = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut aaa, qfi.a()) };

        let mut bbb = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut bbb, qfi.b()) };

        let mut ccc = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut ccc, qfi.c()) };

        let qfi3 = super::BQF::new(
            &mpz_to_bignum(&mut aaa),
            &mpz_to_bignum(&mut bbb),
            &mpz_to_bignum(&mut ccc),
        );

        println!("qfi2={:?}", qfi2);

        println!(
            "norm={:?}, norm={:?}",
            qfi3,
            qfi2.normalize()
                .normalize()
                .normalize()
                .normalize()
                .normalize()
                .normalize()
        );
        assert!(qfi3.equals(&qfi2));
    }

    // TODO test reduce

    #[test]
    fn test_compose() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "219211015245339659606923489058910718059300326777750522511052678189518994793540424066835449513550321156725825999213223313349543556000887051910142135883849284272371752572695727069191024343809964931675197912960671210283614957513396516801587047437150725164417736257899198555560707024341197521616471833363582196266114484743";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "-81545512674410670486936176483359392564511873763104542618989110464161210362234569734599379968174823311559519130742888365338004517383649260007974211906536061358807858053200069189129489104433807401679673870588655431254253093929351305686001886690914729246160866675068623497789888497184276822188325873470578904605036981681";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "397911913619280235397607492118765996458330730701094047301004809429554024195731134066384833319157074257034700248114739847099306289375744183458513982780336763663761492062176352821924273594044341543655974705344319285820641420038042279276523776451336649974742547552579516047919797400613300929758921320569816271189292178311";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let mut qfi = unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum1(&mut a),
            &mpz_to_bignum1(&mut b),
            &mpz_to_bignum1(&mut c),
        );
        unsafe { qfi.normalize_0a() };
        unsafe { qfi.normalize_0a() };

        let mut res = unsafe { QFI::new_0a() };
        //let mutref_res: cpp_core::MutRef<QFI> = unsafe { cpp_core::MutRef::from_raw_ref(&mut res) };

        let cl = unsafe { ClassGroup::new(&disc) };

        unsafe { cl.nucomp(&mut res, &qfi, &qfi) };

        unsafe { res.normalize_0a() };
        let mut d = unsafe { QFI::discriminant(&qfi) };
        assert!(mpz_to_bignum1(&d) == mpz_to_bignum1(&mut disc));

        //let d = mpz_to_bignum(&mut d);

        let mut aaa = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut aaa, res.a()) };

        let mut bbb = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut bbb, res.b()) };

        let mut ccc = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut ccc, res.c()) };

        let qfi3 = super::BQF::new(
            &mpz_to_bignum1(&mut aaa),
            &mpz_to_bignum1(&mut bbb),
            &mpz_to_bignum1(&mut ccc),
        );
        println!("disc={:?}, disc={:?}", d, qfi3.discriminant());

        let comp = qfi2.compose(&qfi2);
        let dd = comp.discriminant();

        assert!(mpz_to_bignum1(&d) == dd);
        println!("D={}", dd);
        println!("qfi2={:?}", qfi2);

        let comp_bicycl = qfi3;

        let comp = reduce2(reduce2(qfi2.compose(&qfi2)));

        println!("compose BICYCL={:?}\ncompose={:?}", comp_bicycl, comp);
        assert!(comp_bicycl.equals(&comp));
    }

    #[test]
    fn test_squaring() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "219211015245339659606923489058910718059300326777750522511052678189518994793540424066835449513550321156725825999213223313349543556000887051910142135883849284272371752572695727069191024343809964931675197912960671210283614957513396516801587047437150725164417736257899198555560707024341197521616471833363582196266114484743";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "-81545512674410670486936176483359392564511873763104542618989110464161210362234569734599379968174823311559519130742888365338004517383649260007974211906536061358807858053200069189129489104433807401679673870588655431254253093929351305686001886690914729246160866675068623497789888497184276822188325873470578904605036981681";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "397911913619280235397607492118765996458330730701094047301004809429554024195731134066384833319157074257034700248114739847099306289375744183458513982780336763663761492062176352821924273594044341543655974705344319285820641420038042279276523776451336649974742547552579516047919797400613300929758921320569816271189292178311";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let mut qfi = unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum1(&mut a),
            &mpz_to_bignum1(&mut b),
            &mpz_to_bignum1(&mut c),
        );
        unsafe { qfi.normalize_0a() };
        unsafe { qfi.normalize_0a() };

        let mut res = unsafe { QFI::new_0a() };
        //let mutref_res: cpp_core::MutRef<QFI> = unsafe { cpp_core::MutRef::from_raw_ref(&mut res) };

        let cl = unsafe { ClassGroup::new(&disc) };

        unsafe { cl.nudupl(&mut res, &qfi) };

        unsafe { res.normalize_0a() };
        let d = unsafe { QFI::discriminant(&qfi) };
        assert!(mpz_to_bignum1(&d) == mpz_to_bignum1(&mut disc));

        //let d = mpz_to_bignum(&mut d);

        let mut aaa = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut aaa, res.a()) };

        let mut bbb = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut bbb, res.b()) };

        let mut ccc = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut ccc, res.c()) };

        let qfi3 = super::BQF::new(
            &mpz_to_bignum1(&mut aaa),
            &mpz_to_bignum1(&mut bbb),
            &mpz_to_bignum1(&mut ccc),
        );
        println!("disc={:?}, disc={:?}", d, qfi3.discriminant());

        let comp = qfi2.compose(&qfi2);
        let dd = comp.discriminant();

        assert!(mpz_to_bignum1(&d) == dd);
        println!("D={}", dd);
        println!("qfi2={:?}", qfi2);

        let dupl_bicycl = qfi3;

        let dupl = reduce2(reduce2(qfi2.double()));

        println!("double BICYCL={:?}\ndouble={:?}", dupl_bicycl, dupl);
        assert!(dupl_bicycl.equals(&dupl));
    }

    #[test]
    fn test_pow() {
        // TODO randomize test, use bicycl keygen function to generate valid binary quadratic form
        let a = "1";

        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(a.as_ptr() as *const c_char, a.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut a = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let b = "1";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(b.as_ptr() as *const c_char, b.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut b = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *b)) };
        println!("sign={:?}", cc);

        let c = "633825300114114700748351602708";
        let s: bicycl::cpp_std::cpp_core::CppBox<String> =
            unsafe { String::from_char_usize(c.as_ptr() as *const c_char, c.len()) };
        let s: Ref<String> = unsafe { Ref::from_raw_ref(&s) };
        let mut c = unsafe { b_i_c_y_c_l::Mpz::from_string(s) };

        let a_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&a) };
        let b_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&b) };
        let c_: cpp_core::Ref<Mpz> = unsafe { cpp_core::Ref::from_raw_ref(&c) };
        let mut qfi = unsafe { bicycl::b_i_c_y_c_l::QFI::new_4a(a_, b_, c_, false) };

        let s = unsafe { Ref::from_raw_ref(&qfi) };

        let mut disc = unsafe { b_i_c_y_c_l::QFI::discriminant(&s) };

        let cc = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *disc)) };
        println!("disc2={:?}", cc);

        let aa = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *a)) };
        println!("aa={:?}", aa);

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            &mpz_to_bignum1(&mut a),
            &mpz_to_bignum1(&mut b),
            &mpz_to_bignum1(&mut c),
        );
        unsafe { qfi.normalize_0a() };
        unsafe { qfi.normalize_0a() };

        let mut res = unsafe { QFI::new_0a() };
        //let mutref_res: cpp_core::MutRef<QFI> = unsafe { cpp_core::MutRef::from_raw_ref(&mut res) };

        let cl = unsafe { ClassGroup::new(&disc) };

        unsafe { cl.nupow_3a(&mut res, &qfi, &disc) };

        unsafe { res.normalize_0a() };
        let d = unsafe { QFI::discriminant(&qfi) };
        assert!(mpz_to_bignum1(&d) == mpz_to_bignum1(&mut disc));

        //let d = mpz_to_bignum(&mut d);

        let mut aaa = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut aaa, res.a()) };

        let mut bbb = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut bbb, res.b()) };

        let mut ccc = unsafe { Mpz::new() };
        let _ = unsafe { Mpz::copy_from_mpz(&mut ccc, res.c()) };

        let qfi3 = super::BQF::new(
            &mpz_to_bignum1(&mut aaa),
            &mpz_to_bignum1(&mut bbb),
            &mpz_to_bignum1(&mut ccc),
        );
        println!("disc={:?}, disc={:?}", d, qfi3.discriminant());

        //let comp = qfi2.compose(&qfi2);
        //let dd = comp.discriminant();

        //assert!(mpz_to_bignum1(&d) == dd);
        //println!("D={}", dd);
        println!("qfi2={:?}", qfi2);

        let dupl_bicycl = qfi3;

        let dupl = reduce2(reduce2(qfi2.pow(&mpz_to_bignum1(&d))));

        println!("double BICYCL={:?}\ndouble={:?}", dupl_bicycl, dupl);
        assert!(dupl_bicycl.equals(&dupl));
    }
}
