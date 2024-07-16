// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use std::fmt::Debug;

use blst::{blst_fr, blst_fr_from_scalar};
use ff::{Field, PrimeField};
use group::Group;
use rand_core::{CryptoRng, RngCore};
use serde::Serialize;

use crate::bqf::BinaryQuadraticForm;
use crate::cl_hsmq::ClHSMq;
use crate::z;
use crate::z::Z;

pub trait Scalar {
    fn random<R: RngCore>(rng: &mut R) -> Self;

    fn add_assign(&mut self, rhs: &Self);

    fn mul_assign(&mut self, rhs: &Self);

    fn pow(&self, exp: u64) -> Self;

    fn modulus_as_z<Z: z::Z>() -> Z;

    fn to_z<Z: z::Z>(&self) -> Z;

    fn from_bytes_be(_: Vec<u8>) -> Self;

    fn sqr(&self) -> Self;

    fn zero() -> Self;
    fn from(n: u64) -> Self;

    fn hash(b: Vec<u8>) -> Self;
}

impl Scalar for blstrs::Scalar {
    fn random<R: RngCore + rand_core::RngCore>(rng: &mut R) -> Self {
        <blstrs::Scalar as Field>::random(rng)
    }

    fn add_assign(&mut self, rhs: &Self) {
        *self += rhs
    }

    fn mul_assign(&mut self, rhs: &Self) {
        *self *= rhs
    }

    fn pow(&self, exp: u64) -> Self {
        <blstrs::Scalar as Field>::pow(self, &[exp])
    }

    fn modulus_as_z<Z: z::Z>() -> Z {
        Z::from_string(&blstrs::Scalar::MODULUS[2..], 16)
    }

    fn to_z<Z: z::Z>(&self) -> Z {
        Z::from_bytes_be(self.to_bytes_be().to_vec())
    }

    fn from_bytes_be(b: Vec<u8>) -> Self {
        assert!(b.len() <= 32, "Input byte array is too long");
        let mut padded = [0u8; 32];
        padded[32 - b.len()..].copy_from_slice(&b);
        Self::from_bytes_be(&padded).unwrap()
    }

    fn sqr(&self) -> Self {
        self.square()
    }

    fn zero() -> Self {
        blstrs::Scalar::ZERO.clone()
    }

    fn from(n: u64) -> Self {
        <blstrs::Scalar as From<u64>>::from(n)
    }

    fn hash(b: Vec<u8>) -> Self {
        let dst = vec![];
        let hash = unsafe { blst::blst_scalar::hash_to(&b, &dst).unwrap() };
        let mut ret = blst::blst_fr::default();
        unsafe {
            blst_fr_from_scalar(&mut ret, &hash);
        }
        <blstrs::Scalar as From<blst_fr>>::from(ret)
    }
}

pub trait EllipticCurve<S>
where
    S: Scalar,
{
    type Scalar;

    fn generator() -> Self;

    fn add_assign(&mut self, other: &Self);

    fn mul_assign(&mut self, scalar: &Self::Scalar);

    fn to_bytes(&self) -> Vec<u8>;

    fn equals(&self, other: &Self) -> bool;

    // neutral element
    fn zero() -> Self;
}

impl<S> EllipticCurve<S> for blstrs::G1Projective
where
    S: Scalar,
{
    type Scalar = blstrs::Scalar;

    fn generator() -> Self {
        group::Group::generator()
    }

    fn add_assign(&mut self, other: &Self) {
        *self += other
    }

    fn mul_assign(&mut self, scalar: &Self::Scalar) {
        *self *= scalar;
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.to_compressed().to_vec()
    }

    fn equals(&self, other: &Self) -> bool {
        self.eq(other)
    }

    fn zero() -> Self {
        blstrs::G1Projective::identity()
    }
}

mod ProofOfCorrectSharing {
    use std::fmt::Debug;
    use std::marker::PhantomData;

    use rand_core::{CryptoRng, RngCore};
    use serde::Serialize;

    use crate::bqf::{BinaryQuadraticForm, BQF};
    use crate::cl_hsmq::power_f;
    use crate::vss::{EllipticCurve, Scalar};
    use crate::z;

    #[derive(Serialize)]
    pub(crate) struct Config<Z>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone,
    {
        n: u32,
        threshold: u32,
        q: Z,
        discriminant: Z,
        generator_H: BQF<Z>,
        generator_F: BQF<Z>,
        _integer_type: PhantomData<Z>,
    }

    pub(crate) struct Witness<Z, S>
    where
        Z: crate::z::Z,
        S: Scalar + Clone,
    {
        pub(crate) s: Vec<S>,
        pub(crate) r: Z,
    }

    fn vec_z_tobytes<S, Z>(v: &Vec<Z>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
        Z: z::Z,
    {
        let res = v.iter().flat_map(|e| e.to_bytes_be()).collect::<Vec<u8>>();

        serializer.serialize_bytes(&res)
    }

    fn vec_ec_tobytes<S, E, Z>(v: &Vec<E>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
        Z: Scalar,
        E: EllipticCurve<Z>,
    {
        let res = v.iter().flat_map(|e| e.to_bytes()).collect::<Vec<u8>>();

        serializer.serialize_bytes(&res)
    }

    #[derive(Serialize)]
    pub struct Instance<Z, E, S>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S> + Clone,
        S: Scalar,
    {
        pub(crate) public_keys: Vec<BQF<Z>>,
        pub(crate) ciphertexts_common: BQF<Z>,
        pub(crate) ciphertexts: Vec<BQF<Z>>,
        #[serde(serialize_with = "vec_ec_tobytes")]
        pub(crate) polynomial_coefficients_commitments: Vec<E>,
        pub(crate) _scalar_type: PhantomData<S>,
        pub(crate) _int_type: PhantomData<Z>,
    }

    #[derive(Clone)]
    pub(crate) struct Proof<Z, E, S>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S>,
        S: Scalar,
    {
        w: BQF<Z>,
        x: E,
        y: BQF<Z>,
        z_r: Z, // integer
        z_s: S, // scalar field
    }

    pub(crate) fn new<
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S, Scalar = S> + std::clone::Clone,
        S: Scalar + Clone,
    >(
        n: u32,
        threshold: u32,
        q: Z,
        discriminant: Z,
        generator_H: BQF<Z>,
        generator_F: BQF<Z>,
    ) -> Config<Z> {
        Config {
            n,
            threshold,
            q,
            discriminant,
            generator_H,
            generator_F,
            _integer_type: Default::default(),
        }
    }

    pub(crate) fn prove<
        R: CryptoRng + RngCore,
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + std::clone::Clone,
        S: Scalar + Clone + Debug,
    >(
        config: &Config<Z>,
        instance: &Instance<Z, E, S>,
        witness: &Witness<Z, S>,
        rng: &mut R,
    ) -> Proof<Z, E, S> {
        let alpha = S::random(rng);
        let rho = Z::sample_range(rng, &Z::zero(), &S::modulus_as_z()); // TODO correct range
        let w = config.generator_H.pow(&rho).reduce();

        let mut x = E::generator();
        x.mul_assign(&alpha);

        let gamma_s = S::hash(bincode::serialize(&instance).unwrap());
        println!("gamma_s 1 = {:?}", gamma_s);
        let gamma_z: Z = S::to_z(&gamma_s);
        let mut gamma_pow = gamma_z.clone();

        let y: BQF<Z> = instance
            .public_keys
            .clone()
            .into_iter()
            .fold(config.generator_H.identity(), |acc, pk| {
                let res = acc.compose(&pk.pow(&gamma_pow)).reduce();
                gamma_pow = gamma_pow.mul(&gamma_z);
                res
            })
            .pow(&rho)
            .reduce()
            .compose(&power_f(
                &config.generator_F,
                &config.discriminant,
                &config.q,
                &S::to_z(&alpha),
            ))
            .reduce();

        let mut b: Vec<u8> = vec![];
        b.append(&mut gamma_z.to_bytes_be());
        b.append(&mut bincode::serialize(&w).unwrap());
        b.append(&mut x.to_bytes());
        b.append(&mut bincode::serialize(&y).unwrap());
        let gamma_prime_s = S::hash(b);
        let gamma_prime_z = S::to_z(&gamma_prime_s);
        println!("gamma_prime_s 1 = {:?}", gamma_prime_s);

        let z_r = witness.r.mul(&gamma_prime_z).add(&rho);
        let mut gamma_pow = gamma_s.clone();
        let mut z_s = witness
            .s
            .clone()
            .into_iter()
            .fold(S::zero(), |mut acc, mut s| {
                s.mul_assign(&gamma_pow);
                gamma_pow.mul_assign(&gamma_s);
                acc.add_assign(&s);
                acc
            });
        z_s.mul_assign(&gamma_prime_s);
        z_s.add_assign(&alpha);

        Proof { w, x, y, z_r, z_s }
    }

    pub(crate) fn verify<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + std::clone::Clone,
        S: Scalar + Clone + Debug,
    >(
        config: &Config<Z>,
        instance: &Instance<Z, E, S>,
        proof: &Proof<Z, E, S>,
    ) -> bool {
        let mut b = vec![];
        b.append(&mut bincode::serialize(&instance).unwrap());
        let gamma_s = S::hash(b);
        println!("gamma_s 2 = {:?}", gamma_s);
        let gamma_z: Z = S::to_z(&gamma_s);

        // TODO refactor
        let mut b: Vec<u8> = vec![];
        b.append(&mut gamma_z.to_bytes_be());
        b.append(&mut bincode::serialize(&proof.w).unwrap());
        b.append(&mut proof.x.to_bytes());
        b.append(&mut bincode::serialize(&proof.y).unwrap());
        let gamma_prime_s = S::hash(b);
        let gamma_prime_z = S::to_z(&gamma_prime_s);
        println!("gamma_prime_s 2 = {:?}", gamma_prime_s);

        let lhs = proof
            .w
            .compose(&instance.ciphertexts_common.pow(&gamma_prime_z))
            .reduce();
        let rhs = config.generator_H.pow(&proof.z_r).reduce();

        if !lhs.equals(&rhs) {
            return false;
        }

        /*let (_, mut lhs) = instance
            .polynomial_coefficients_commitments
            .clone()
            .into_iter()
            .enumerate()
            .reduce(|(_, acc), (j, mut a_j)| {
                let pow_gamma = gamma_s.pow(j as u64);
                let exp = (1..=config.n).fold(S::zero(), |mut acc, i| {
                    acc.add_assign(&pow_gamma);
                    acc.mul_assign(&S::from(i.pow(j as u32) as u64));
                    acc
                });
                a_j.mul_assign(&exp);
                a_j.add_assign(&acc);
                (j, a_j)
            })
            .unwrap();
        lhs.mul_assign(&gamma_prime_s);
        lhs.add_assign(&proof.x);

        let mut rhs = E::generator();
        rhs.mul_assign(&proof.z_s);
        if !lhs.equals(&rhs) {
            return false;
        }

        let mut gamma_pow = gamma_z.clone();
        let (_, mut lhs) = instance
            .ciphertexts
            .clone()
            .into_iter()
            .enumerate()
            .reduce(|(_, acc), (i, c)| {
                let res = c.pow(&gamma_pow);
                gamma_pow = gamma_pow.mul(&gamma_z);
                (i, acc.compose(&res))
            })
            .unwrap();
        lhs = lhs.pow(&gamma_prime_z).compose(&proof.y);

        let (_, mut rhs) = instance
            .public_keys
            .clone()
            .into_iter()
            .enumerate()
            .reduce(|(_, acc), (i, k)| {
                let res = k.pow(&gamma_pow);
                gamma_pow = gamma_pow.mul(&gamma_z);
                (i, acc.compose(&res))
            })
            .unwrap();
        rhs = rhs
            .pow(&proof.z_r)
            .compose(&config.generator_F.pow(&S::to_z(&proof.z_s))); // TODO: One can use faster power_of_f here
        if !lhs.equals(&rhs) {
            return false;
        }*/
        true
    }
}

pub trait Polynomial<S>
where
    S: Scalar,
{
    type Scalar;
    fn random<R: RngCore + CryptoRng>(rng: &mut R, degree: usize) -> Self;

    fn set_coefficient(&mut self, index: usize, scalar: &Self::Scalar);

    fn evaluate(&self, point: &Self::Scalar) -> Self::Scalar;

    fn coefficients(&self) -> Vec<Self::Scalar>;
}

impl<S> Polynomial<S> for Vec<blstrs::Scalar>
where
    S: Scalar,
{
    type Scalar = blstrs::Scalar;

    fn random<R: RngCore + CryptoRng>(rng: &mut R, degree: usize) -> Self {
        //vec![<blstrs::Scalar as Field>::random(rng); degree + 1]
        (0..=degree).map(|_| Scalar::random(rng)).collect()
    }

    fn set_coefficient(&mut self, index: usize, scalar: &Self::Scalar) {
        self[index] = *scalar
    }

    fn evaluate(&self, point: &Self::Scalar) -> Self::Scalar {
        // with Horner's method
        let mut result = Scalar::zero();
        for coef in self.iter().rev() {
            result = result * point + coef;
        }
        result
    }

    fn coefficients(&self) -> Vec<Self::Scalar> {
        self.clone()
    }
}

pub mod VerifiableSecretSharingPrimitives {
    use std::fmt::Debug;

    use rand_core::{CryptoRng, RngCore};
    use serde::Serialize;

    use crate::bqf::{BinaryQuadraticForm, BQF};
    use crate::cl_hsmq::{ClHSMq, ClHSMqInstance};
    use crate::vss::ProofOfCorrectSharing::Proof;
    use crate::vss::{EllipticCurve, Polynomial, ProofOfCorrectSharing, Scalar};

    pub fn share<
        R: RngCore + CryptoRng,
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S> + Debug,
    >(
        rng: &mut R,
        n: usize,
        threshold: usize,
        s: &S,
    ) -> (Vec<S>, Vec<E>) {
        let mut p = P::random(rng, threshold);
        p.set_coefficient(0, s);
        let s_i = (1..=n).map(|x| p.evaluate(&S::from(x as u64))).collect();
        let cmt = p
            .coefficients()
            .iter()
            .take(threshold + 1)
            .map(|c| {
                let mut g = E::generator();
                g.mul_assign(c);
                g
            })
            .collect();
        (s_i, cmt)
    }

    pub fn encrypt<
        R: RngCore + CryptoRng,
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        rng: &mut R,
        enc_scheme: &ClHSMqInstance<Z, BQF<Z>>,
        pocs_config: &ProofOfCorrectSharing::Config<Z>,
        cmt: &Vec<E>,
        shares: &Vec<S>,
        pub_keys: &Vec<BQF<Z>>,
    ) -> (BQF<Z>, Vec<BQF<Z>>, ProofOfCorrectSharing::Proof<Z, E, S>) {
        let (common, encryptions, r) = enc_scheme.encrypt_batch(
            pub_keys,
            &shares.iter().map(|s| s.to_z()).collect::<Vec<Z>>(),
            rng,
        );
        let instance = ProofOfCorrectSharing::Instance {
            public_keys: pub_keys.clone(),
            ciphertexts_common: common.clone(),
            ciphertexts: encryptions.clone(),
            polynomial_coefficients_commitments: cmt.clone(),
            _scalar_type: Default::default(),
            _int_type: Default::default(),
        };
        let witness = ProofOfCorrectSharing::Witness {
            s: shares.clone(),
            r,
        };
        let proof: Proof<Z, E, S> =
            ProofOfCorrectSharing::prove(pocs_config, &instance, &witness, rng);

        (common, encryptions, proof)
    }

    pub fn verify<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        pocs_config: &ProofOfCorrectSharing::Config<Z>,
        cmt: &Vec<E>,
        common_enc: &BQF<Z>,
        enc_shares: &Vec<BQF<Z>>,
        pub_keys: &Vec<BQF<Z>>,
        proof: &Proof<Z, E, S>,
    ) -> bool {
        let instance = ProofOfCorrectSharing::Instance {
            public_keys: pub_keys.clone(),
            ciphertexts_common: common_enc.clone(),
            ciphertexts: enc_shares.clone(),
            polynomial_coefficients_commitments: cmt.clone(),
            _scalar_type: Default::default(),
            _int_type: Default::default(),
        };
        ProofOfCorrectSharing::verify(pocs_config, &instance, &proof)
    }

    pub fn decrypt<
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone,
        P: Polynomial<S, Scalar = S>,
    >(
        enc_scheme: &ClHSMqInstance<Z, BQF<Z>>,
        sk: &Z,
        common_enc: &BQF<Z>,
        enc_share: &BQF<Z>,
    ) -> S {
        S::from_bytes_be(Z::to_bytes_be(
            &enc_scheme.decrypt(sk, &(common_enc.clone(), enc_share.clone())),
        ))
    }

    pub fn verify_commitment<
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone,
        P: Polynomial<S, Scalar = S>,
    >(
        cmt: &Vec<E>,
        share_index: usize,
        share: &S,
    ) -> bool {
        let mut lhs = E::generator();
        lhs.mul_assign(share);

        let rhs = cmt
            .clone()
            .into_iter()
            .enumerate()
            .fold(E::zero(), |acc, (j, mut a_j)| {
                a_j.mul_assign(&S::from(share_index.pow(j as u32) as u64));
                a_j.add_assign(&acc);
                a_j
            });

        lhs.equals(&rhs)
    }
}

pub mod NIVSS {
    use std::fmt::Debug;

    use rand_core::{CryptoRng, RngCore};
    use serde::Serialize;

    use crate::bqf::{BinaryQuadraticForm, BQF};
    use crate::cl_hsmq::{ClHSMq, ClHSMqInstance, SecurityLevel};
    use crate::vss::{
        EllipticCurve, Polynomial, ProofOfCorrectSharing, Scalar, VerifiableSecretSharingPrimitives,
    };

    pub struct InitialConfig<Z>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone,
    {
        pub(crate) n: usize,
        pub(crate) threshold: usize,
        pub(crate) security_level: SecurityLevel,
        pub(crate) q: Z,
        pub(crate) encryption_scheme: ClHSMqInstance<Z, BQF<Z>>,
        pub(crate) index: usize,
    }

    pub struct Config<Z>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
    {
        pub(crate) public_parameters: PublicParameters<Z>,
        pub(crate) index: usize,
        pub(crate) public_key: BQF<Z>,
        pub(crate) secret_key: Z,
    }

    #[derive(Clone)]
    pub struct PublicParameters<Z>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
    {
        pub(crate) n: usize,
        pub(crate) threshold: usize,
        pub(crate) public_keys: Vec<BQF<Z>>,
        pub(crate) security_level: SecurityLevel,
        pub(crate) q: Z,
        pub(crate) discriminant: Z,
        pub(crate) encryption_scheme: ClHSMqInstance<Z, BQF<Z>>,
    }

    #[derive(Clone)]
    pub struct Dealing<E, S, Z>
    where
        Z: crate::z::Z + std::fmt::Debug + Clone,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone,
    {
        common_encryption: BQF<Z>,
        pub(crate) encryptions: Vec<BQF<Z>>,
        correct_sharing_proof: ProofOfCorrectSharing::Proof<Z, E, S>,
        pub(crate) cmt: Vec<E>,
    }

    pub fn generate_dealing<
        R: RngCore + CryptoRng,
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S> + Debug,
    >(
        rng: &mut R,
        pp: &PublicParameters<Z>,
    ) -> Dealing<E, S, Z> {
        let s = S::random(rng);
        let (shares, cmt) =
            VerifiableSecretSharingPrimitives::share::<R, Z, E, S, P>(rng, pp.n, pp.threshold, &s);
        let pocs = ProofOfCorrectSharing::new::<Z, E, S>(
            pp.n as u32,
            pp.threshold as u32,
            pp.q.clone(),
            pp.discriminant.clone(),
            pp.encryption_scheme.generator_h(),
            pp.encryption_scheme.generator_f(),
        );
        println!(
            "pk len = {}, shares len={}, cmt len={}",
            pp.public_keys.len(),
            shares.len(),
            cmt.len()
        );
        let (common, encryptions, correct_sharing_proof) =
            VerifiableSecretSharingPrimitives::encrypt::<R, Z, E, S, P>(
                rng,
                &pp.encryption_scheme,
                &pocs,
                &cmt,
                &shares,
                &pp.public_keys,
            );
        Dealing {
            common_encryption: common,
            encryptions,
            correct_sharing_proof,
            cmt,
        }
    }

    pub fn verify_dealing<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        pp: &PublicParameters<Z>,
        dealing: &Dealing<E, S, Z>,
    ) -> bool {
        let pocs_cfg = ProofOfCorrectSharing::new::<Z, E, S>(
            pp.n as u32,
            pp.threshold as u32,
            pp.q.clone(),
            pp.discriminant.clone(),
            pp.encryption_scheme.generator_h(),
            pp.encryption_scheme.generator_f(),
        );
        VerifiableSecretSharingPrimitives::verify::<Z, E, S, P>(
            &pocs_cfg,
            &dealing.cmt,
            &dealing.common_encryption,
            &dealing.encryptions,
            &pp.public_keys,
            &dealing.correct_sharing_proof,
        )
    }

    pub fn verify_dealing_output_secret_key_share<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        cfg: &Config<Z>,
        dealing: &Dealing<E, S, Z>,
    ) -> Option<S> {
        if !verify_dealing::<Z, E, S, P>(&cfg.public_parameters, &dealing) {
            return None;
        }

        let s = VerifiableSecretSharingPrimitives::decrypt::<Z, E, S, P>(
            &cfg.public_parameters.encryption_scheme,
            &cfg.secret_key,
            &dealing.common_encryption,
            &dealing.encryptions[cfg.index],
        );

        Some(s)
    }
}

// TODO assert threshold < n/2
// TODO modify interface to compile valid dealings into a single transcript
// then provide the ability to derive the master public key and share of the corresponding
// master secret key?
// TODO key resharing with Shamir secret sharing of the existing shares,
// this allows to keep the existing master public key
// TODO implement section mitigating the biasing of master public key
mod NIDKG {
    use std::fmt::Debug;

    use rand_core::{CryptoRng, RngCore};
    use serde::Serialize;

    use crate::bqf::BQF;
    use crate::cl_hsmq::ClHSMq;
    use crate::nizk_dlog::nizk_dlog_verify;
    use crate::vss::NIVSS::{Config, Dealing, InitialConfig, PublicParameters};
    use crate::vss::{EllipticCurve, Polynomial, Scalar, NIVSS};

    // key pair+proof of possession
    pub fn generate_key_pair<R: RngCore + CryptoRng, Z: crate::z::Z + std::fmt::Debug + Clone>(
        cfg: &InitialConfig<Z>,
        rng: &mut R,
    ) -> (BQF<Z>, Z, crate::nizk_dlog::NizkDlogProof<Z>) {
        let (public_key, secret_key) = cfg.encryption_scheme.keygen(rng);
        let proof = crate::nizk_dlog::nizk_dlog_prove(
            &cfg.encryption_scheme.generator_h(),
            &public_key,
            &secret_key,
            &cfg.encryption_scheme.class_number_bound_h(), // TODO fix bound
        );
        println!("k1 = {:?}, p1 = {:?}", public_key, proof);
        (public_key, secret_key, proof)
    }

    pub fn verify_public_keys<Z: crate::z::Z + std::fmt::Debug + Clone>(
        cfg: &InitialConfig<Z>,
        public_keys: &Vec<BQF<Z>>,
        proofs_of_possession: &Vec<crate::nizk_dlog::NizkDlogProof<Z>>,
    ) -> bool {
        public_keys.iter().zip(proofs_of_possession).all(|(k, p)| {
            println!("k = {:?}, p = {:?}", k, p);
            nizk_dlog_verify(
                p,
                &cfg.encryption_scheme.generator_h(),
                k,
                &cfg.encryption_scheme.class_number_bound_h(),
            ) // TODO fix bound
        })
    }

    pub fn generate_dealing<
        R: RngCore + CryptoRng,
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S> + Debug,
    >(
        pp: &PublicParameters<Z>,
        rng: &mut R,
    ) -> Dealing<E, S, Z> {
        crate::vss::NIVSS::generate_dealing::<R, Z, E, S, P>(rng, pp)
    }

    pub fn verify_dealing<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        pp: &PublicParameters<Z>,
        dealing: &Dealing<E, S, Z>,
    ) -> bool {
        crate::vss::NIVSS::verify_dealing::<Z, E, S, P>(pp, dealing)
    }

    // TODO add function that computes the master public key and public key shares only
    pub fn aggregate_dealings<
        Z: crate::z::Z + std::fmt::Debug + Clone + Serialize,
        E: EllipticCurve<S, Scalar = S> + Clone,
        S: Scalar + Clone + Debug,
        P: Polynomial<S, Scalar = S>,
    >(
        cfg: &Config<Z>,
        dealings: &Vec<(usize, Dealing<E, S, Z>)>,
    ) -> (E, E, S, Vec<E>, Vec<E>) {
        let res: Vec<S> = dealings
            .iter()
            .filter_map(|(i, d)| {
                NIVSS::verify_dealing_output_secret_key_share::<Z, E, S, P>(cfg, d)
            })
            .collect();
        let secret_key_share = res
            .into_iter()
            .reduce(|acc, x| {
                let mut acc = acc.clone();
                acc.add_assign(&x);
                acc
            })
            .unwrap()
            .clone(); // TODO better error handling
        let mut public_key_share = E::generator();
        public_key_share.mul_assign(&secret_key_share);

        let master_public_key = dealings.iter().fold(E::zero(), |mut acc, (_, d)| {
            acc.add_assign(&d.cmt[0]);
            acc
        });

        /*let b_j = dealings.iter().fold(E::zero(), |mut acc, (_, d)| {
            acc.add_assign(&d.cmt[cfg.index]);
            acc
        });*/

        let mut public_poly = vec![E::zero(); cfg.public_parameters.threshold + 1];
        for (_, d) in dealings {
            for (i, x) in d.cmt.iter().enumerate() {
                public_poly[i].add_assign(x)
            }
        }

        let indices: Vec<usize> = dealings.iter().map(|(i, _)| *i).collect();
        // TODO we can do Pippenger here, though in our case the amount of nodes is small (< 10)
        // TODO horner method to evaluate public polynomial at j
        let public_key_shares: Vec<E> = (1..=cfg.public_parameters.n)
            .map(|j| {
                let mut j_pow = 1usize;
                public_poly.iter().fold(E::zero(), |acc, x| {
                    let mut acc = acc.clone();
                    let mut x = x.clone();
                    x.mul_assign(&S::from(j_pow as u64));
                    acc.add_assign(&x);
                    j_pow *= j;
                    acc
                })
            })
            .collect::<Vec<E>>();

        // TODO store in struct
        (
            master_public_key,
            public_key_share,
            secret_key_share,
            public_key_shares,
            public_poly,
        )
    }
}

#[cfg(test)]
mod tests {
    use blstrs::{G1Projective, Scalar};
    use rand_core::OsRng;
    use rug::Integer;

    use crate::cl_hsmq::{ClHSMq, ClHSMqInstance};
    use crate::vss::NIVSS::{Dealing, PublicParameters};

    #[test]
    fn f() {
        let number_of_participants: usize = 10;
        let threshold: usize = 4; // < number_of_participants / 2
        let security_level = crate::cl_hsmq::SecurityLevel::SecLvl112;
        let q = rug::Integer::from_str_radix(
            "52435875175126190479447740508185965837690552500527637822603658699938581184513",
            10,
        )
        .unwrap();

        let encryption_scheme = ClHSMqInstance::new(q.clone(), security_level, &mut OsRng);
        let vss_config_0 = crate::vss::NIVSS::InitialConfig {
            n: number_of_participants,
            threshold,
            security_level: security_level.clone(),
            q: q.clone(),
            encryption_scheme: encryption_scheme.clone(),
            index: 1,
        };

        let mut pks = vec![];
        let mut sks = vec![];
        let mut pops = vec![];

        let mut configs = vec![];

        for i in 0..number_of_participants {
            let vss_config_i = crate::vss::NIVSS::InitialConfig {
                n: number_of_participants,
                threshold,
                security_level: security_level.clone(),
                q: q.clone(),
                encryption_scheme: encryption_scheme.clone(),
                index: i + 1,
            };

            let (pk, sk, pop) = crate::vss::NIDKG::generate_key_pair(&vss_config_i, &mut OsRng);

            pks.push(pk);
            sks.push(sk);
            pops.push(pop);
        }

        let pp: PublicParameters<Integer> = PublicParameters {
            n: number_of_participants,
            threshold,
            public_keys: pks.clone(),
            security_level: security_level.clone(),
            q: q.clone(),
            discriminant: encryption_scheme.discriminant.clone(),
            encryption_scheme: encryption_scheme.clone(),
        };

        for i in 0..number_of_participants {
            let config_i = crate::vss::NIVSS::Config {
                public_parameters: pp.clone(),
                public_key: pks[i].clone(),
                secret_key: sks[i].clone(),
                index: i + 1,
            };
            configs.push(config_i);
        }

        assert!(crate::vss::NIDKG::verify_public_keys(
            &vss_config_0,
            &vec![pks[0].clone()],
            &vec![pops[0].clone()],
        ));

        let mut dealings: Vec<Dealing<G1Projective, Scalar, Integer>> = Vec::new();
        //generating dealings for nodes
        for i in 0..=threshold {
            let dealing: Dealing<G1Projective, Scalar, Integer> =
                crate::vss::NIDKG::generate_dealing::<
                    _,
                    Integer,
                    blstrs::G1Projective,
                    blstrs::Scalar,
                    Vec<blstrs::Scalar>,
                >(&pp, &mut OsRng);
            dealings.push(dealing.clone());

            assert!(crate::vss::NIDKG::verify_dealing::<
                Integer,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(&pp, &dealing));
        }

        for i in 0..=threshold {
            let (
                master_public_key,
                public_key_share,
                secret_key_share,
                public_key_shares,
                public_poly,
            ) = crate::vss::NIDKG::aggregate_dealings::<
                Integer,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(
                &configs[i],
                &dealings
                    .iter()
                    .enumerate()
                    .map(|(i, d)| (i + 1, d.clone()))
                    .collect(),
            );

            println!(
                "mpk={:?},pks={:?},sks={:?},pkss={:?}",
                master_public_key, public_key_share, secret_key_share, public_key_shares
            );

            // g^a0 g^a1 g^a2 g^a3 g^a4
        }
    }
}
