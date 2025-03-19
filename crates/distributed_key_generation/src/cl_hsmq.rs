// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use std::fmt::Debug;

#[cfg(feature = "random")]
use rand_core::CryptoRng;

use crate::bqf::BinaryQuadraticForm;
use crate::z;

#[cfg(feature = "random")]
use crate::z::Randomizable;

#[derive(Clone, Copy)]
pub enum SecurityLevel {
    SecLvl112,
    SecLvl128,
}

impl From<SecurityLevel> for u32 {
    fn from(val: SecurityLevel) -> Self {
        match val {
            SecurityLevel::SecLvl112 => 112,
            SecurityLevel::SecLvl128 => 128,
        }
    }
}

fn discriminant_bit_size(security_level: SecurityLevel) -> u32 {
    match security_level {
        SecurityLevel::SecLvl112 => 1348,
        SecurityLevel::SecLvl128 => 1827,
    }
}

/// Linear homomorphic encryption based on the hidden subgroup membership problem
/// with the CL framework
pub trait ClHSMq<
    Z: crate::z::Z + std::fmt::Debug + Clone + std::cmp::PartialEq,
    BQF: BinaryQuadraticForm<Z> + Clone,
>
{
    // setup
    #[cfg(feature = "random")]
    fn new<Rng: CryptoRng + rand_core::RngCore, S>(
        security_level: SecurityLevel,
        rng: &mut Rng,
    ) -> Self
    where
        Z: Randomizable,
        S: crate::scalar::Scalar + Clone + Debug;

    // generates an asymmetric key pair (pk, sk)
    #[cfg(feature = "random")]
    fn keygen<Rng: CryptoRng + rand_core::RngCore>(&self, rng: &mut Rng) -> (BQF, Z)
    where
        Z: Randomizable;

    // ciphertext is composed of two binary quadratic forms
    #[cfg(feature = "random")]
    fn encrypt<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_key: &BQF,
        cleartext: &Z,
        rng: &mut Rng,
    ) -> (BQF, BQF)
    where
        Z: Randomizable;
    #[cfg(feature = "random")]
    fn encrypt_batch<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_keys: &[BQF],
        cleartexts: &[Z],
        rng: &mut Rng,
    ) -> (BQF, Vec<BQF>, Z)
    where
        Z: Randomizable;

    fn decrypt(&self, private_key: &Z, ciphertext: &(BQF, BQF)) -> Z;

    fn add_cleartexts(&self, public_key: BQF, cleartext1: Z, cleartext2: Z) -> Z;

    fn scale_cleartext(&self, public_key: BQF, cleartext: Z, scaling_factor: Z) -> Z;
    #[cfg(feature = "random")]
    fn add_ciphertexts<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF)
    where
        Z: Randomizable;

    #[cfg(feature = "random")]
    fn scale_ciphertext<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext: (BQF, BQF),
        scaling_factor: Z,
    ) -> (BQF, BQF)
    where
        Z: Randomizable;

    fn generator_f(&self) -> BQF;

    fn generator_h(&self) -> BQF;

    fn class_number_bound_h(&self) -> Z;
}

// We'll take q=order of BLS12-381 scalar field F_r (multiplicative order)
// so that the messages matches scalar elements (F_r)^
#[derive(Clone)]
pub struct ClHSMqInstance<Z, BQF>
where
    BQF: BinaryQuadraticForm<Z> + Clone,
    Z: z::Z + std::fmt::Debug + std::clone::Clone + std::cmp::PartialEq,
{
    pub(crate) class_number_h_bound: Z,
    pub(crate) generator_f: BQF,
    pub(crate) generator_h: BQF,
    pub(crate) discriminant: Z,
    pub(crate) q: Z,
}

// Note that with a security level of 112 we can support up to 3 multiplications
// of 1348-bit integers fitting inside 4096 bit ints
#[cfg(feature = "random")]
fn sample_random_p<Z: crate::z::Z + PartialEq, R: CryptoRng + rand_core::RngCore>(
    q: &Z,
    q_bit_size: u32,
    disc_bit_size: u32,
    rng: &mut R,
) -> Z
where
    Z: Randomizable,
{
    assert!(disc_bit_size > q_bit_size);
    let mut p = Z::sample_bits(disc_bit_size - q_bit_size, rng);
    let three = Z::from(3);
    let four = Z::from(4);
    let target_residue = if q.take_mod(&four).eq(&three) {
        Z::from(1)
    } else {
        three
    };
    p = p.next_prime();
    while !p.take_mod(&four).eq(&target_residue) || q.kronecker(&p) != -1 {
        p = p.next_prime();
    }
    p
}

impl<Z, BQF> ClHSMq<Z, BQF> for ClHSMqInstance<Z, BQF>
where
    Z: crate::z::Z + std::fmt::Debug + std::clone::Clone + std::cmp::PartialEq,
    BQF: BinaryQuadraticForm<Z> + Clone + std::fmt::Debug,
{
    #[cfg(feature = "random")]
    fn new<Rng: CryptoRng + rand_core::RngCore, S>(
        security_level: SecurityLevel,
        rng: &mut Rng,
    ) -> Self
    where
        Z: Randomizable,
        S: crate::scalar::Scalar + Clone + Debug,
    {
        // "We denote η(λ) the bitsize of a fundamental discriminant ∆K such that"
        // "the computation of the class number h(∆K) and computation of discrete logarithms"
        // "in Cl(∆K ) takes at least 2^λ operations (see Table 2 for concrete sizes)."
        // lambda: security level
        // for lambda = 128, eta(lambda)=1827 bits
        // mu: q is a mu-bits prime
        // with BLS12-381 mu=255=log2up(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        let q: Z = S::modulus_as_z();
        let q_bit_size = q.bit_size() as u32;
        let disc_bit_size = discriminant_bit_size(security_level);
        assert!(q_bit_size >= <SecurityLevel as Into<u32>>::into(security_level));
        assert!(q_bit_size < disc_bit_size);
        let p = sample_random_p::<Z, Rng>(&q, q_bit_size, disc_bit_size, rng);
        let fundamental_discriminant = p.mul(&q).neg();
        let discriminant_conductor_q = q.sqr().mul(&fundamental_discriminant);
        // G_hat commutative group of unknown order M * s_hat
        // G_hat isomorphic to F x H where F has order M
        let generator_h = compute_generator_h(&q, &discriminant_conductor_q);
        let mut c = Z::from(1).sub(&fundamental_discriminant);
        c.divide_by_4_exact();
        let generator_f = BQF::new(&q.sqr(), &q, &c).reduce();
        let class_number_h_bound = BQF::class_number_bound(&fundamental_discriminant);

        let lhs = q.sqr().mul(&Z::from(4));
        let rhs = Z::from(1).sub(&fundamental_discriminant);
        // we're not in the large message variant case
        assert!(lhs.compare(&rhs).is_le());
        ClHSMqInstance {
            class_number_h_bound,
            generator_f, // generator of subgroup of G_hat of order M
            generator_h, // generator of subgroup H of G_hat
            discriminant: fundamental_discriminant,
            q,
        }
    }

    #[cfg(feature = "random")]
    fn keygen<Rng: CryptoRng + rand_core::RngCore>(&self, rng: &mut Rng) -> (BQF, Z)
    where
        Z: Randomizable,
    {
        let secret_key = Z::sample_range(rng, &Z::from(0), &self.class_number_h_bound);
        let public_key = self.generator_h.pow(&secret_key);
        (public_key, secret_key)
    }

    #[cfg(feature = "random")]
    fn encrypt<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_key: &BQF,
        cleartext: &Z,
        rng: &mut Rng,
    ) -> (BQF, BQF)
    where
        Z: Randomizable,
    {
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_h_bound);
        let c1 = self.generator_h.pow(&r);
        let f_m = power_f(&self.generator_f, &self.discriminant, &self.q, cleartext);
        let pk_r = public_key.pow(&r);
        let c2 = f_m.compose(&pk_r);
        (c1, c2)
    }

    #[cfg(feature = "random")]
    fn encrypt_batch<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_keys: &[BQF],
        cleartexts: &[Z],
        rng: &mut Rng,
    ) -> (BQF, Vec<BQF>, Z)
    where
        Z: Randomizable,
    {
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_h_bound);
        let c1 = self.generator_h.pow(&r);

        let c2s = cleartexts
            .iter()
            .zip(public_keys)
            .map(|(c, k)| {
                power_f(&self.generator_f, &self.discriminant, &self.q, c).compose(&k.pow(&r))
            })
            .collect();

        (c1, c2s, r)
    }

    fn decrypt(&self, private_key: &Z, ciphertext: &(BQF, BQF)) -> Z {
        let dlog = ciphertext
            .1
            .compose(&ciphertext.0.reduce().pow(private_key).inverse());
        self.dlog_solve_f(&dlog).take_mod(&self.q)
    }

    fn add_cleartexts(&self, _public_key: BQF, _cleartext1: Z, _cleartext2: Z) -> Z {
        todo!()
    }

    fn scale_cleartext(&self, _public_key: BQF, _cleartext: Z, _scaling_factor: Z) -> Z {
        todo!()
    }

    #[cfg(feature = "random")]
    fn add_ciphertexts<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF)
    where
        Z: Randomizable,
    {
        let c1_prod = ciphertext1.0.compose(&ciphertext2.0);
        let c2_prod = ciphertext1.1.compose(&ciphertext2.1);
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_h_bound);
        let c1 = c1_prod.compose(&self.generator_h.pow(&r));
        let c2 = c2_prod.compose(&public_key.pow(&r));
        (c1, c2)
    }

    #[cfg(feature = "random")]
    fn scale_ciphertext<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext: (BQF, BQF),
        scaling_factor: Z,
    ) -> (BQF, BQF)
    where
        Z: Randomizable,
    {
        let c1_scaled = ciphertext.0.pow(&scaling_factor);
        let c2_scaled = ciphertext.1.pow(&scaling_factor);
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_h_bound);
        let c1 = c1_scaled.compose(&self.generator_h.pow(&r));
        let c2 = c2_scaled.compose(&public_key.pow(&r));
        (c1, c2)
    }

    fn generator_f(&self) -> BQF {
        self.generator_f.clone()
    }

    fn generator_h(&self) -> BQF {
        self.generator_h.clone()
    }

    fn class_number_bound_h(&self) -> Z {
        self.class_number_h_bound.clone()
    }
}

pub fn power_f<Z, BQF>(_f: &BQF, disc: &Z, q: &Z, cleartext: &Z) -> BQF
where
    BQF: BinaryQuadraticForm<Z> + Clone + std::fmt::Debug,
    Z: std::clone::Clone + std::fmt::Debug + z::Z + std::cmp::PartialEq,
{
    let a;
    let b;
    let mut c = cleartext.invert_mod(q).unwrap();
    if !c.is_odd() {
        c = c.sub(q);
    }

    /*    [ q^2, Lm*q, ((Lm*q)^2-Delta_q)/(4*q^2) ]
     * =  [ q^2, Lm*q, ((Lm*q)^2-q^2*Delta_K)/(4*q^2) ]
     * =  [ q^2, Lm*q, (Lm^2-Delta_K)/4 ]
     */
    a = q.sqr(); /* a = q^2 */
    b = c.mul(q); /* b = Lm*q */
    c = c.sqr();
    c = c.sub(disc);
    c.divide_by_4_exact(); /* c = (Lm^2-Delta_K)/4 */

    BQF::new(&a, &b, &c).reduce()
}

impl<Z, BQF> ClHSMqInstance<Z, BQF>
where
    Z: crate::z::Z + Debug + Clone + std::cmp::PartialEq,
    BQF: BinaryQuadraticForm<Z> + Clone,
{
    pub fn dlog_solve_f(&self, fm: &BQF) -> Z {
        let mut m = Z::from(0u64);
        if !fm.equals(&fm.identity()) {
            m = fm.b().divide_exact(&self.q).invert_mod(&self.q).unwrap();
        }
        m
    }
}

fn compute_generator_h<
    Z: std::clone::Clone + std::fmt::Debug + z::Z + std::cmp::PartialEq,
    BQF: BinaryQuadraticForm<Z> + Clone,
>(
    q: &Z,
    discriminant_conductor_q: &Z,
) -> BQF {
    let mut l = Z::from(2);
    while discriminant_conductor_q.kronecker(&l) != 1 {
        l = l.next_prime();
    }
    BQF::prime_form(discriminant_conductor_q, &l)
        .double()
        .pow(q)
}

#[test]
#[cfg(feature = "random")]
#[cfg(feature = "gmp")]
fn test_encrypt_decrypt_identity() {
    use crate::bqf::BQF;
    use crate::scalar::Scalar;
    use rand_core::OsRng;
    use rug::Integer;

    let security_level = SecurityLevel::SecLvl112;
    let mut rng = &mut OsRng;

    // println!("before new");
    let instance = ClHSMqInstance::new::<&mut OsRng, blstrs::Scalar>(security_level, &mut rng);

    // println!("before keygen");

    let (public_key, private_key) = instance.keygen::<&mut OsRng>(&mut rng);

    // println!("generated keys = pk={:?} sk={}", public_key, private_key);
    let rd = blstrs::Scalar::random(&mut OsRng);
    // println!("rd={:?}", rd);
    let cleartext = blstrs::Scalar::to_z(&rd); // Sample cleartext
    // println!("cleartext={:?}", cleartext);
    let ciphertext: (BQF<Integer>, BQF<Integer>) =
        instance.encrypt::<&mut OsRng>(&public_key, &cleartext, &mut rng);

    // println!("generated ciphertext {:?}", ciphertext);

    // println!("before decrypt");

    let decrypted = instance.decrypt(&private_key, &ciphertext);

    assert_eq!(
        decrypted, cleartext,
        "Decrypted text does not match the original cleartext"
    );

    let pk_sk: Vec<(BQF<Integer>, Integer)> = (1..10)
        .map(|_| instance.keygen::<&mut OsRng>(&mut rng))
        .collect();
    let public_keys: Vec<BQF<Integer>> = pk_sk.iter().map(|(k, _)| k.clone()).collect();
    let private_keys: Vec<Integer> = pk_sk.iter().map(|(_, s)| s.clone()).collect();
    let cleartexts: Vec<Integer> = (1..10).map(|_| blstrs::Scalar::to_z(&rd)).collect(); // Sample cleartext
    // println!("cleartext={:?}", cleartext);
    let (common, ciphertexts, _) =
        instance.encrypt_batch::<&mut OsRng>(&public_keys, &cleartexts, &mut rng);

    // println!("generated ciphertexts {:?}", ciphertexts);

    // println!("before decrypt");

    for i in 0..9 {
        let decrypted =
            instance.decrypt(&private_keys[i], &(common.clone(), ciphertexts[i].clone()));

        assert_eq!(
            decrypted, cleartext,
            "Decrypted text does not match the original cleartext"
        );
    }
}
