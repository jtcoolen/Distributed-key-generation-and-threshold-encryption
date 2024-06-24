//use rand_core::CryptoRng;

use crate::z;

trait BinaryQuadraticForm<Z>
where
    Z: crate::z::Z,
{
    fn new(a: Z, b: Z, c: Z) -> Self;

    //fn random<R: CryptoRng>(rng: &mut R) -> Self;

    fn equals(&self, other: &Self) -> bool;

    fn discriminant(&self) -> Z;

    fn identity(&self) -> Self;

    fn normalize(&self) -> Self;

    fn reduce(&self) -> Self;

    fn compose(self, other: Self) -> Self;

    fn double(self) -> Self;

    fn pow(self, exponent: Z) -> Self;

    fn inverse(self) -> Self;
}

#[derive(Debug)]
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
        let mut c = b.sub(&disc).clone();
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
        let b = remainder.clone();
        let mut remainder = self.b.add(&remainder);
        remainder.divide_by_2();
        let c = self.c.sub(&quotient.mul(&remainder));
        BQF {
            a: self.a.clone(),
            b,
            c,
        }
    }

    fn reduce(&self) -> Self {
        let mut n = self.normalize();
        while n.a.less_than_abs(&n.c) {
            n = n.rho();
        }
        if n.a.eq_abs(&n.c) && !n.b.is_positive() {
            n.b.oppose()
        }
        n
    }

    fn compose(self, other: Self) -> Self {
        // https://www.ams.org/journals/mcom/2003-72-244/S0025-5718-03-01518-7/S0025-5718-03-01518-7.pdf
        todo!()
    }

    fn double(self) -> Self {
        todo!()
    }

    fn pow(self, exponent: Z) -> Self {
        todo!()
    }

    // TODO (harder): implement compose, double, pow
    // reverse-engineer their implementation
    //https://gite.lirmm.fr/crypto/bicycl/-/blob/master/src/bicycl/qfi.inl?ref_type=heads#L621

    // TODO (less of a priority): compute class group number
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
        let mut limbs = unsafe { bicycl::cpp_vec_to_rust(&Mpz::mpz_to_b_i_g_bytes(&mut *n)) };

        limbs.reverse();
        Bignum4096 {
            positive: unsafe { Mpz::sgn(&n) == 1 },
            limbs: convert(limbs),
        }
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

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            mpz_to_bignum(&mut a),
            mpz_to_bignum(&mut b),
            mpz_to_bignum(&mut c),
        );

        println!("qfi2={:?}", qfi2);
        println!("qfi2={:?}, disc={:?}", qfi2, qfi2.discriminant());
        assert!(qfi2.discriminant() == mpz_to_bignum(&mut disc))
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

        println!("disc {:?}", mpz_to_bignum(&mut disc));

        let qfi2 = super::BQF::new(
            mpz_to_bignum(&mut a),
            mpz_to_bignum(&mut b),
            mpz_to_bignum(&mut c),
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
        assert!(mpz_to_bignum2(&a) == qfi2.identity().a);
        assert!(mpz_to_bignum2(&b) == qfi2.identity().b);
        assert!(mpz_to_bignum2(&c) == qfi2.identity().c);

        println!("qfi2={:?}", qfi2);
        println!("qfi2={:?}, id={:?}", qfi2, qfi2.identity());
        assert!(qfi2.discriminant() == mpz_to_bignum(&mut disc))
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
            mpz_to_bignum(&mut a),
            mpz_to_bignum(&mut b),
            mpz_to_bignum(&mut c),
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
            mpz_to_bignum(&mut aaa),
            mpz_to_bignum(&mut bbb),
            mpz_to_bignum(&mut ccc),
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
}
