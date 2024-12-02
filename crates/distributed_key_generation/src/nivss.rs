// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use std::fmt::Debug;
use std::time::Instant;

#[cfg(feature = "random")]
use rand_core::{CryptoRng, RngCore};

#[cfg(feature = "random")]
use crate::z::Randomizable;

use serde::Serialize;

use crate::bqf::BQF;
use crate::cl_hsmq::{ClHSMq, ClHSMqInstance};
use crate::elliptic_curve::EllipticCurve;
use crate::polynomial::Polynomial;

/// Returns an
#[cfg(feature = "random")]
pub fn share<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Debug,
>(
    rng: &mut R,
    n: usize,
    threshold: usize,
    s: &S,
) -> (Vec<S>, Vec<E>, P) {
    let mut p: P = P::random(rng, threshold);
    let p_len = p.coefficients().len();
    assert_eq!(p_len, threshold + 1);
    p.set_coefficient(0, s);
    let s_i = (1..=n).map(|i| p.evaluate(&S::from(i as u64))).collect();
    let cmt: Vec<E> = p
        .coefficients()
        .iter()
        .map(|a| {
            let mut g = E::generator();
            g.mul_assign(a);
            g
        })
        .collect();
    assert_eq!(cmt.len(), p_len);
    for (p, c) in p.coefficients().iter().zip(&cmt) {
        let mut e = E::generator().clone();
        e.mul_assign(p);
    }
    (s_i, cmt, p)
}

#[cfg(feature = "random")]
pub fn encrypt<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S>,
>(
    rng: &mut R,
    enc_scheme: &ClHSMqInstance<Z, BQF<Z>>,
    pocs_config: &crate::nizkp_secret_sharing::Config<Z>,
    cmt: &[E],
    shares: &[S],
    pub_keys: &[BQF<Z>],
) -> (
    BQF<Z>,
    Vec<BQF<Z>>,
    crate::nizkp_secret_sharing::Proof<Z, E, S>,
)
where
    Z: Randomizable,
{
    let (common, encryptions, r) = enc_scheme.encrypt_batch::<R>(
        pub_keys,
        &shares.iter().map(|s| s.to_z()).collect::<Vec<Z>>(),
        rng,
    );
    let instance = crate::nizkp_secret_sharing::Instance {
        public_keys: pub_keys.to_owned(),
        ciphertexts_common: common.clone(),
        ciphertexts: encryptions.clone(),
        polynomial_coefficients_commitments: cmt.to_owned(),
    };
    let witness = crate::nizkp_secret_sharing::Witness {
        s: shares.to_owned(),
        r,
    };
    let proof: crate::nizkp_secret_sharing::Proof<Z, E, S> =
        crate::nizkp_secret_sharing::prove::<R, Z, E, S>(pocs_config, &instance, &witness, rng);

    (common, encryptions, proof)
}

pub fn verify<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S>,
>(
    pocs_config: &crate::nizkp_secret_sharing::Config<Z>,
    cmt: &[E],
    common_enc: &BQF<Z>,
    enc_shares: &[BQF<Z>],
    pub_keys: &[BQF<Z>],
    proof: &crate::nizkp_secret_sharing::Proof<Z, E, S>,
) -> bool {
    let instance = crate::nizkp_secret_sharing::Instance {
        public_keys: pub_keys.to_owned(),
        ciphertexts_common: common_enc.clone(),
        ciphertexts: enc_shares.to_owned(),
        polynomial_coefficients_commitments: cmt.to_owned(),
    };
    crate::nizkp_secret_sharing::verify(pocs_config, &instance, proof)
}

#[cfg(feature = "random")]
pub fn decrypt<
    Z: crate::z::Z + std::fmt::Debug + Clone + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone,
    P: Polynomial<Scalar = S>,
>(
    enc_scheme: &ClHSMqInstance<Z, BQF<Z>>,
    sk: &Z,
    common_enc: &BQF<Z>,
    enc_share: &BQF<Z>,
) -> S {
    // number is positive
    S::from_bytes_be(
        Z::to_bytes_be(&enc_scheme.decrypt(sk, &(common_enc.clone(), enc_share.clone()))).0,
    )
}

pub fn verify_commitment<
    Z: crate::z::Z + std::fmt::Debug + Clone,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone,
    P: Polynomial<Scalar = S>,
>(
    cmt: &[E],
    share_index: usize,
    share: &S,
) -> bool {
    let mut lhs = E::generator();
    lhs.mul_assign(share);

    let rhs = cmt
        .iter()
        .cloned()
        .enumerate()
        .fold(E::zero(), |acc, (j, mut a_j)| {
            a_j.mul_assign(&S::from((share_index + 1).pow(j as u32) as u64));
            a_j.add_assign(&acc);
            a_j
        });

    lhs.equals(&rhs)
}

use crate::cl_hsmq::SecurityLevel;

#[derive(Clone)]
pub struct PublicParameters<Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
{
    pub(crate) n: usize,
    pub(crate) threshold: usize,
    pub(crate) security_level: SecurityLevel,
    pub(crate) q: Z,
    pub(crate) discriminant: Z,
    pub(crate) encryption_scheme: ClHSMqInstance<Z, BQF<Z>>,
}

#[derive(Clone)]
pub struct Dealing<E, S, Z, P>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone,
    P: Polynomial<Scalar = S> + Clone + Debug,
{
    pub(crate) common_encryption: BQF<Z>,
    pub(crate) encryptions: Vec<BQF<Z>>,
    pub(crate) correct_sharing_proof: crate::nizkp_secret_sharing::Proof<Z, E, S>,
    pub(crate) cmt: Vec<E>,
    pub(crate) p: P,
}

#[cfg(feature = "random")]
pub fn generate_dealing_resharing<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    rng: &mut R,
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
    s: &S,
) -> Dealing<E, S, Z, P>
where
    Z: Randomizable,
{
    let (shares, cmt, p) = share::<R, Z, E, S, P>(rng, pp.n, pp.threshold, s);
    let pocs = crate::nizkp_secret_sharing::Config::new::<E, S>(
        pp.n as u32,
        pp.threshold as u32,
        pp.q.clone(),
        pp.discriminant.clone(),
        pp.encryption_scheme.generator_h(),
        pp.encryption_scheme.generator_f(),
    );
    let (common, encryptions, correct_sharing_proof) = encrypt::<R, Z, E, S, P>(
        rng,
        &pp.encryption_scheme,
        &pocs,
        &cmt,
        &shares,
        public_keys,
    );
    Dealing {
        common_encryption: common,
        encryptions,
        correct_sharing_proof,
        cmt,
        p,
    }
}

#[cfg(feature = "random")]
pub fn generate_dealing<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    rng: &mut R,
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
) -> Dealing<E, S, Z, P>
where
    Z: Randomizable,
{
    let s = S::random(rng);
    generate_dealing_resharing::<R, Z, E, S, P>(rng, pp, public_keys, &s)
}

pub fn verify_dealing<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
    dealing: &Dealing<E, S, Z, P>,
) -> bool {
    let pocs_cfg = crate::nizkp_secret_sharing::Config::new::<E, S>(
        pp.n as u32,
        pp.threshold as u32,
        pp.q.clone(),
        pp.discriminant.clone(),
        pp.encryption_scheme.generator_h(),
        pp.encryption_scheme.generator_f(),
    );
    verify::<Z, E, S, P>(
        &pocs_cfg,
        &dealing.cmt,
        &dealing.common_encryption,
        &dealing.encryptions,
        public_keys,
        &dealing.correct_sharing_proof,
    )
}

#[cfg(feature = "random")]
pub fn verify_dealing_output_secret_key_share<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    index: usize,
    secret_key: &Z,
    public_keys: &[BQF<Z>],
    dealing: &Dealing<E, S, Z, P>,
) -> Option<S> {
    let now = Instant::now();
    if !verify_dealing::<Z, E, S, P>(pp, public_keys, dealing) {
        return None;
    }
    println!("verify dealing {:?}", now.elapsed());

    let now = Instant::now();
    let s = decrypt::<Z, E, S, P>(
        &pp.encryption_scheme,
        secret_key,
        &dealing.common_encryption,
        &dealing.encryptions[index],
    );
    println!("decrypt for dealing {:?}", now.elapsed());

    let now = Instant::now();
    if !verify_commitment::<Z, E, S, P>(&dealing.cmt, index, &s) {
        return None;
    }
    println!("verify comm for dealing {:?}", now.elapsed());

    Some(s)
}
