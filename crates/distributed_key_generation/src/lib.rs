// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

//! Client for a non-interactive distributed key generation.
//!
//! Multiple parties collectively compute an asymmetric key pair
//! whose private key can be recovered with a threshold of private
//! key shares.
//!
//! The process is non-interactive, meaning that there is no back-and-forth
//! exchanges between participants: they share their dealing to the other
//! participants, and once a threshold of dealings has been exchanged, and
//! compiled in a transcript they agree on, each party can derive a common public
//! key and a share of the associated private key.
//!
//! In addition, a resharing protocol allows to refresh the private key
//! shares set when new participants join, or to update the threshold,
//! while keeping the existing public key. Another reason for resharing
//! is proactive security: in case shares are leaked one wants to reset
//! those keys to not compromise the security of the system.
//!
//! ## Usage
//!
//! ```rust
//! use rand::rngs::OsRng;
//! use rug::Integer;
//! use blstrs::G1Projective;
//! use distributed_key_generation::setup;
//! use distributed_key_generation::generate_key_pair;
//! use distributed_key_generation::verify_public_key;
//! use distributed_key_generation::generate_dealing;
//! use distributed_key_generation::verify_dealing;
//! use distributed_key_generation::SecurityLevel;
//! 
//! let number_of_participants: usize = 10;
//! let threshold: usize = 4; // < number_of_participants / 2
//! let security_level = SecurityLevel::SecLvl112;
//! let pp = setup::<_, blstrs::Scalar, _>(
//!     number_of_participants,
//!     threshold,
//!     security_level,
//!     &mut OsRng,
//! );

//! let mut pks = vec![];
//! let mut sks = vec![];
//! let mut pops = vec![];

//! for _ in 0..number_of_participants {
//!     let (pk, sk, pop) = generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);

//!     pks.push(pk);
//!     sks.push(sk);
//!     pops.push(pop);
//! }

//! // verify keys
//! assert!(verify_public_key(&pp, &pks[0], &pops[0],));

//! let mut dealings = Vec::new();

//! //generating dealings for nodes
//! for _ in 0..=threshold {
//!     let dealing = generate_dealing::<
//!         _,
//!         _,
//!         blstrs::G1Projective,
//!         blstrs::Scalar,
//!         Vec<blstrs::Scalar>,
//!     >(&pp, &pks, &mut OsRng);
//!     dealings.push(dealing);
//! }

//! for i in 0..=0 {
//!     assert!(verify_dealing::<
//!         Integer,
//!         G1Projective,
//!         blstrs::Scalar,
//!         Vec<blstrs::Scalar>,
//!     >(&pp, &pks, &dealings[i]));
//! }
//! ```

mod bqf;
mod cl_hsmq;

#[cfg(feature = "gmp")]
pub mod mpz;

mod elliptic_curve;
mod nivss;
mod nizk_dlog;
mod nizkp_secret_sharing;
mod polynomial;
mod scalar;
pub mod signed4096;
pub mod z;
use crate::cl_hsmq::ClHSMq;

use std::fmt::Debug;

use crate::bqf::BQF;
use crate::elliptic_curve::EllipticCurve;

#[doc(inline)]
pub use crate::nivss::{Dealing, PublicParameters};

use crate::polynomial::Polynomial;
#[cfg(feature = "random")]
use crate::z::Randomizable;
pub use cl_hsmq::SecurityLevel;
#[cfg(feature = "random")]
use rand_core::{CryptoRng, RngCore};
use rand_core::{CryptoRngCore, OsRng};
use serde::Serialize;

/// Public Key Infrastructure (PKI) setup for the DKG:

/// Randomized key generation that returns an asymmetric key pair
/// for a given participant identified by `index`.
#[cfg(feature = "random")]
pub fn generate_key_pair<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + std::cmp::PartialEq + Serialize,
>(
    pp: &PublicParameters<Z>,
    rng: &mut R,
) -> (BQF<Z>, Z, crate::nizk_dlog::NizkDlogProof<Z>)
where
    Z: Randomizable,
{
    let (public_key, secret_key) = pp.encryption_scheme.keygen::<R>(rng);
    let proof = crate::nizk_dlog::nizk_dlog_prove::<Z, R>(
        &pp.encryption_scheme.generator_h(),
        &public_key,
        &secret_key,
        &pp.encryption_scheme.class_number_bound_h(),
        rng,
    );
    (public_key, secret_key, proof)
}

/// Returns true if and only if (with high probability) the public keys are accompanied
/// with valid proofs of possession, guaranteeing that each issuer of a public key knows
/// the corresponding secret key.
pub fn verify_public_key<
    Z: crate::z::Z + std::fmt::Debug + Clone + std::cmp::PartialEq + Serialize,
>(
    pp: &PublicParameters<Z>,
    public_key: &BQF<Z>,
    proof_of_possession: &crate::nizk_dlog::NizkDlogProof<Z>,
) -> bool {
    crate::nizk_dlog::nizk_dlog_verify(
        proof_of_possession,
        &pp.encryption_scheme.generator_h(),
        public_key,
        &pp.encryption_scheme.class_number_bound_h(),
    )
}

/// Randomized dealing algorithm that given a `pp.threshold` and a set of public keys
/// `pp.public_keys` generates a dealing.
#[cfg(feature = "random")]
pub fn generate_dealing<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
    rng: &mut R,
) -> Dealing<E, S, Z, P>
where
    Z: Randomizable,
{
    crate::nivss::generate_dealing::<R, Z, E, S, P>(rng, pp, public_keys)
}

/// Randomized dealing algorithm that given a `pp.threshold` and a set of public keys
/// `pp.public_keys` generates a dealing for a participant by doing a secret sharing of its
/// secret key share `dkg_result.secret_key_share`.
/// This allows participants to a previous distributed key generation to generate
/// a dealing for a new set of participants, so that they can in turn derive their
/// own secret key share of the existing master public key.
#[cfg(feature = "random")]
pub fn generate_dealing_resharing<
    R: RngCore + CryptoRng,
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
    rng: &mut R,
    dkg_result: &DKGPrivateResult<E, S>,
) -> Dealing<E, S, Z, P>
where
    Z: Randomizable,
{
    crate::nivss::generate_dealing_resharing::<R, Z, E, S, P>(
        rng,
        pp,
        public_keys,
        &dkg_result.secret_key_share,
    )
}

/// Deterministic dealing verification.
pub fn verify_dealing<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    public_keys: &[BQF<Z>],
    dealing: &Dealing<E, S, Z, P>,
) -> bool {
    crate::nivss::verify_dealing::<Z, E, S, P>(pp, public_keys, dealing)
}

pub struct DKGPublicResult<E> {
    master_public_key: E,
    public_key_shares: Vec<E>,
    public_poly: Vec<E>,
}

/// Derives from all valid dealings the result of the distributed
/// key generation, including the master public key
/// and master public key shares.
pub fn aggregate_verified_dealings_public<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    dealings: &[Dealing<E, S, Z, P>],
) -> DKGPublicResult<E> {
    let master_public_key2 = dealings.iter().fold(E::zero(), |mut acc, d| {
        acc.add_assign(&d.cmt[0]);
        acc
    });

    let mut public_poly = vec![E::zero(); pp.threshold + 1];
    for (_i, d) in dealings.iter().enumerate() {
        for (j, x) in d.cmt.iter().enumerate() {
            public_poly[j].add_assign(x);
        }
    }

    let mut coeffs = P::random(&mut OsRng, pp.threshold);
    for i in 0..=pp.threshold {
        coeffs.set_coefficient(i, &S::zero());
    }

    for i in 0..=pp.threshold {
        for d in dealings {
            let mut c = coeffs.coefficients()[i].clone();
            c.add_assign(&d.p.coefficients()[i].clone());
            coeffs.set_coefficient(i, &c)
        }
    }
    // println!("sum poly = {:?}", coeffs);

    for i in 1..=pp.n {
        let sk_i = coeffs.evaluate(&S::from(i as u64));
        // println!("eval(sum poly,index={})={:?}", i, sk_i);
        let mut e = E::generator().clone();
        e.mul_assign(&sk_i);
        // println!("pk={:?}", e);
    }

    let sk_i = coeffs.evaluate(&S::zero());
    // println!("sum poly(0) = {:?}", sk_i);
    let mut e = E::generator().clone();
    e.mul_assign(&sk_i);
    // println!(
    //     "MPK={:?}, mpk={:?}",
    //     e.to_bytes(),
    //     master_public_key2.to_bytes()
    // );

    // println!("len dealings = {}", dealings.len());

    let master_public_key = public_poly[0].clone();
    // println!(
    //     "mpk1={:?}, mpk2={:?}",
    //     master_public_key, master_public_key2
    // );

    let public_key_shares: Vec<E> = (1..=pp.n)
        .map(|j| {
            let scalars: Vec<S> = (0..=pp.threshold)
                .map(|i| S::from(j.pow(i as u32) as u64))
                .collect();
            E::multiexp(&public_poly, &scalars)
        })
        .collect();

    DKGPublicResult {
        master_public_key,
        public_key_shares,
        public_poly,
    }
}

#[cfg(feature = "random")]
pub struct DKGPrivateResult<E, S> {
    master_public_key: E,
    public_key_share: E,
    secret_key_share: S,
    public_key_shares: Vec<E>,
    public_poly: Vec<E>,
}

/// Derives from all valid dealings the result of the distributed
/// key generation, including the master public key,
/// master public key shares, and share of the master secret key
/// for the given participant identified by `index`.
#[cfg(feature = "random")]
pub fn aggregate_dealings<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + std::cmp::PartialEq,
    E: EllipticCurve<S = S> + Clone + Debug,
    S: crate::scalar::Scalar + Clone + Debug,
    P: Polynomial<Scalar = S> + Clone + Debug,
>(
    pp: &PublicParameters<Z>,
    index: usize,
    secret_key: &Z,
    public_keys: &[BQF<Z>],
    dealings: &[(usize, Dealing<E, S, Z, P>)],
) -> DKGPrivateResult<E, S> {
    let extracted_secret_key_shares_with_dealings = dealings
        .to_owned()
        .clone()
        .into_iter()
        .filter_map(|(_i, d)| {
            // TODO get rid of index
            crate::nivss::verify_dealing_output_secret_key_share::<Z, E, S, P>(
                pp,
                index,
                secret_key,
                public_keys,
                &d,
            )
            .map(|s| (s, d))
        })
        .collect::<Vec<_>>();

    // println!(
    //     "extracted_secret_key_shares_with_dealings.len={}",
    //     extracted_secret_key_shares_with_dealings.len()
    // );
    assert!(extracted_secret_key_shares_with_dealings.len() == pp.n);

    let (decrypted_evaluations, verified_dealings): (Vec<S>, Vec<Dealing<E, S, Z, P>>) =
        extracted_secret_key_shares_with_dealings
            .into_iter()
            .unzip();

    // sum decrypted evaluations of pp.n polynomials at the same point `index`
    let secret_key_share = decrypted_evaluations
        .into_iter()
        .enumerate()
        .fold(S::zero(), |mut acc, (_i, x)| {
            acc.add_assign(&x);
            acc
        })
        .clone();
    let mut public_key_share = E::generator();
    public_key_share.mul_assign(&secret_key_share);

    let DKGPublicResult {
        master_public_key,
        public_key_shares,
        public_poly,
    } = aggregate_verified_dealings_public::<Z, E, S, P>(pp, &verified_dealings);

    //let mut public_key_shares = public_key_shares.clone();
    //public_key_shares[index] = public_key_share.clone();

    DKGPrivateResult {
        master_public_key,
        public_key_share,
        secret_key_share,
        public_key_shares,
        public_poly,
    }
}

use crate::cl_hsmq::ClHSMqInstance;

pub fn setup<
    R: CryptoRng + CryptoRngCore,
    S: crate::scalar::Scalar + Clone + Debug,
    Z: z::Z + Randomizable + Serialize + Debug + Clone + PartialEq,
>(
    number_of_participants: usize,
    threshold: usize,
    security_level: SecurityLevel,
    rng: &mut R,
) -> PublicParameters<Z> {
    let encryption_scheme: ClHSMqInstance<Z, BQF<Z>> =
        ClHSMqInstance::new::<R, S>(security_level, rng);
    PublicParameters {
        n: number_of_participants,
        threshold,
        security_level,
        q: S::modulus_as_z(),
        encryption_scheme: encryption_scheme.clone(),
        discriminant: encryption_scheme.discriminant.clone(),
    }
}

#[cfg(test)]
#[cfg(feature = "random")]
#[cfg(feature = "gmp")]
pub(crate) mod tests {
    use blstrs::{G1Affine, G1Projective, Scalar};
    use ff::Field;
    use group::Group;
    use rand::RngCore;
    use rand_core::OsRng;
    use rug::Integer;
    use threshold_encryption::ciphertext::Ciphertext;

    use crate::bqf::{BinaryQuadraticForm, BQF};
    use crate::cl_hsmq::{ClHSMqInstance, SecurityLevel};
    use crate::nivss::{Dealing, PublicParameters};
    use crate::nizkp_secret_sharing::Proof;
    use crate::polynomial::Polynomial;
    use crate::signed4096::Bignum;
    use crate::z::Z;
    use crate::{setup, DKGPrivateResult, DKGPublicResult};
    use std::ops::Mul;
    use std::time::Instant;

    fn int_to_bignum<Z: crate::z::Z>(p0: &Z) -> Bignum {
        let (p0, sgn) = p0.to_bytes_be();

        Bignum::from_bytes_be(p0, sgn)
    }

    #[cfg(feature = "gmp")]
    fn convert_bqf(p0: &BQF<Integer>) -> BQF<Bignum> {
        BQF::new(
            &int_to_bignum(&p0.a()),
            &int_to_bignum(&p0.b()),
            &int_to_bignum(&p0.c()),
        )
    }

    #[cfg(feature = "gmp")]
    fn convert_bqf2(p0: &BQF<Bignum>) -> BQF<Integer> {
        let a = Bignum::to_bytes_be(&p0.a());
        let b = Bignum::to_bytes_be(&p0.b());
        let c = Bignum::to_bytes_be(&p0.c());
        BQF::new(
            &Integer::from_bytes_be(a.0, a.1),
            &Integer::from_bytes_be(b.0, b.1),
            &Integer::from_bytes_be(c.0, c.1),
        )
    }

    fn recover_master_secret_key_from_secret_key_shares(
        secret_key_shares: &[Scalar],
        threshold: usize,
        number_of_participants: usize,
    ) -> Scalar {
        let mut res = Scalar::ZERO.clone();
        let indices: Vec<Scalar> = (1..=(threshold + 1))
            .map(|e| {
                let mut e = e.to_le_bytes().to_vec();
                e.resize(32, 0);
                Scalar::from_bytes_le(&e.try_into().unwrap()).unwrap()
            })
            .collect();
        let lagrange_coeffs: Vec<Scalar> = indices
            .iter()
            .map(|j| threshold_encryption::helpers::lagrange_coeff(&indices, j))
            .collect();
        for (i, k) in secret_key_shares.iter().take(threshold + 1).enumerate() {
            res += k * lagrange_coeffs[i]
        }
        res
    }

    #[test]
    fn test_dkg() {
        let v: Vec<u8> = vec![
            13, 231, 152, 219, 237, 17, 174, 101, 21, 100, 62, 5, 114, 48, 186, 249, 169, 146, 201,
            195, 232, 93, 14, 88, 82, 39, 98, 122, 15, 141, 56, 255,
        ];

        let s = Scalar::from_bytes_be(&v.try_into().unwrap()).unwrap();
        // println!("s={:?}", s);
        // println!(
        //     "id={:?}, mpk={:?}",
        //     G1Projective::generator(),
        //     G1Projective::generator().mul(s)
        // );

        let p: Vec<blstrs::Scalar> = Polynomial::random(&mut OsRng, 4);
        // println!("p={:?}", p);
        let mut evals = vec![];
        for i in 1..=11 {
            let mut v = [0u8; 32];
            v[0] = i;
            evals.push(p.evaluate(&Scalar::from_bytes_be(&v).unwrap()))
        }
        // println!("evals={:?}", evals);

        let p0 = recover_master_secret_key_from_secret_key_shares(&evals, 4, 10);
        // println!("p0={:?}", p0);

        let number_of_participants: usize = 10;
        let threshold: usize = 3; // < number_of_participants / 2
        let security_level = crate::cl_hsmq::SecurityLevel::SecLvl112;

        let start = Instant::now();
        let pp = setup::<_, blstrs::Scalar, _>(
            number_of_participants,
            threshold,
            security_level,
            &mut OsRng,
        );
        let duration = start.elapsed();

        println!("Setup: {:?}", duration);

        let mut pks = vec![];
        let mut sks = vec![];
        let mut pops = vec![];

        for _ in 0..number_of_participants {
            let start = Instant::now();
            let (pk, sk, pop) = crate::generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);
            let duration = start.elapsed();

            println!("Generate key: {:?}", duration);
            pks.push(pk);
            sks.push(sk);
            pops.push(pop);
        }

        // verify keys
        for i in 0..number_of_participants {
            let start = Instant::now();
            assert!(crate::verify_public_key(&pp, &pks[i], &pops[i]));
            let duration = start.elapsed();

            println!("Verify key: {:?}", duration);
        }

        let mut dealings = Vec::new();

        //generating dealings for nodes
        for _ in 0..number_of_participants {
            let start = Instant::now();
            let dealing = crate::generate_dealing::<
                _,
                _,
                blstrs::G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(&pp, &pks, &mut OsRng);
            let duration = start.elapsed();

            println!("Generate dealing: {:?}", duration);
            dealings.push(dealing);
        }

        // Done in aggregate anyway (not public aggregate though)
        /*for i in 0..number_of_participants {
            assert!(crate::verify_dealing::<
                Integer,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(&pp, &pks, &dealings[i]));
        }*/

        let mut sks_: Vec<Scalar> = vec![];
        let mut mpks = vec![];
        for i in 0..number_of_participants {
            //threshold {
            let start = Instant::now();
            let DKGPrivateResult {
                master_public_key,
                public_key_share,
                secret_key_share,
                public_key_shares,
                public_poly,
            } = crate::aggregate_dealings::<
                Integer,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(
                &pp,
                i,
                &sks[i],
                &pks,
                &dealings
                    .clone()
                    .into_iter()
                    .enumerate()
                    .collect::<Vec<(usize, Dealing<G1Projective, Scalar, Integer, _>)>>(),
            );

            let duration = start.elapsed();

            println!("Aggregate dealing: {:?}", duration);

            assert_eq!(
                G1Affine::from(master_public_key),
                G1Affine::from_compressed_unchecked(&master_public_key.to_compressed()).unwrap()
            );
            // println!(
            //     "mpk={:?}\n\npks={:?}\n\nsks={:?}\n\npkss={:?}\n\npub poly={:?}",
            //     G1Affine::from(master_public_key),
            //     G1Affine::from(public_key_share),
            //     secret_key_share,
            //     public_key_shares
            //         .iter()
            //         .map(|e| G1Affine::from(e))
            //         .collect::<Vec<G1Affine>>(),
            //     public_poly
            //         .iter()
            //         .map(|e| G1Affine::from(e))
            //         .collect::<Vec<G1Affine>>()
            // );

            // println!(
            //     "mpk={:?}\n\npks={:?}\n\nsks={:?}\n\npkss={:?}\n\npub poly={:?}",
            //     G1Affine::from(master_public_key).to_compressed(),
            //     G1Affine::from(public_key_share).to_compressed(),
            //     secret_key_share,
            //     public_key_shares
            //         .iter()
            //         .map(|e| G1Affine::from(e).to_compressed())
            //         .collect::<Vec<[u8; 48]>>(),
            //     public_poly
            //         .iter()
            //         .map(|e| G1Affine::from(e).to_compressed())
            //         .collect::<Vec<[u8; 48]>>(),
            // );

            mpks.push(master_public_key.clone());

            sks_.push(secret_key_share.clone());
            // TODO check mpk by reconstructing it from public key shares
        }
        assert_eq!(sks_.len(), number_of_participants);
        // Check mpk/associated master secret key by reconstructing it from secret key shares
        let sk = recover_master_secret_key_from_secret_key_shares(
            &sks_,
            threshold,
            number_of_participants,
        );

        let mpk2 = G1Projective::generator().mul(sk);

        // println!(
        //     "sk={:?} , g^sk={:?}, mpk={:?}",
        //     sk,
        //     G1Affine::from(mpk2),
        //     G1Affine::from(mpks[0])
        // );
        //assert_eq!(G1Projective::generator().mul(sk), mpks[0]);

        let mut key = [0u8; 32];
        OsRng::fill_bytes(&mut OsRng, &mut key);
        let key_hash = threshold_encryption::helpers::keccak_256(&key);

        let enc2 = Ciphertext::encrypt(&key, key_hash, &mpks[0].into()).unwrap();
        let uhw = enc2.verify_decode().unwrap();

        // println!("sks len = {}", sks_.len());

        let dsks: Vec<threshold_encryption::decryption_share::DSH> = sks_
            .clone()
            .iter()
            .enumerate()
            .map(|(i, sk)| {
                let sk = threshold_encryption::key_shares::SecretKeyShare::new(
                    i as u8, *sk, //int_to_scalar(sk),
                );
                threshold_encryption::decryption_share::DecryptionShare::create(&uhw, &sk)
                    .decode_unchecked()
                    .unwrap()
            })
            .collect();

        let valid_key_shares: Vec<threshold_encryption::decryption_share::DSH> =
            dsks.into_iter().take(threshold + 1 as usize).collect();
        // println!("len key shares = {}", valid_key_shares.len());

        let combined_key = threshold_encryption::decryption_share::DecryptionShare::combine(
            valid_key_shares.as_slice(),
        )
        .unwrap();

        // println!(
        //     "key={:?}",
        //     threshold_encryption::helpers::elgamal_apply(&enc2.encrypted_key, &combined_key)
        // );

        let dec_key = enc2.decrypt_checked(&combined_key).unwrap();
        // println!("dec={:?}", dec_key);
        assert_eq!(dec_key, key);
    }

    fn int_to_scalar(i: &Integer) -> Scalar {
        let mut e = i.to_bytes_be().0;
        // println!("e before={:?}, e sign={}", e, i.to_bytes_be().1);
        e.reverse();
        e.resize(32, 0);
        e.reverse();
        // println!("e after={:?}", e);
        Scalar::from_bytes_le(&e.try_into().unwrap()).unwrap()
    }

    /*#[test]
    fn wasm() {
        let number_of_participants: usize = 10;
        let threshold: usize = 4; // < number_of_participants / 2
        let security_level = crate::cl_hsmq::SecurityLevel::SecLvl112;

        let pp = setup::<_, blstrs::Scalar, _>(
            number_of_participants,
            threshold,
            security_level,
            &mut OsRng,
        );

        let mut pks = vec![];
        let mut sks = vec![];
        let mut pops = vec![];

        for _i in 0..number_of_participants {
            let (pk, sk, pop) = crate::generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);

            pks.push(pk);
            sks.push(sk);
            pops.push(pop);
        }

        for i in 0..number_of_participants {
            let bn = int_to_bignum(&sks[i]);
            let res = bn.to_bytes_be();
            let res = Integer::from_bytes_be(res.0, res.1);

            assert_eq!(res, sks[i]);

            let pk = convert_bqf(&pks[i]);
            let res = convert_bqf2(&pk);
            assert_eq!(res, pks[i]);
        }

        assert!(crate::verify_public_key(&pp, &pks[0], &pops[0],));

        let mut dealings: Vec<Dealing<G1Projective, Scalar, Integer, _>> = Vec::new();
        //generating dealings for nodes
        for _i in 0..=0 {
            //=threshold {
            let dealing: Dealing<G1Projective, Scalar, Integer, _> =
                crate::generate_dealing::<
                    _,
                    _,
                    blstrs::G1Projective,
                    blstrs::Scalar,
                    Vec<blstrs::Scalar>,
                >(&pp, &pks, &mut OsRng);
            dealings.push(dealing.clone());

            let d = dealing.clone();
            let d = Dealing {
                common_encryption: convert_bqf(&d.common_encryption),
                encryptions: d.encryptions.iter().map(convert_bqf).collect(),
                correct_sharing_proof: crate::nizkp_secret_sharing::Proof {
                    w: convert_bqf(&d.correct_sharing_proof.w),
                    x: d.correct_sharing_proof.x,
                    y: convert_bqf(&d.correct_sharing_proof.y),
                    z_r: int_to_bignum(&d.correct_sharing_proof.z_r),
                    z_s: d.correct_sharing_proof.z_s,
                },
                cmt: d.cmt.clone(),
                p: d.p.clone(),
            };

            let pks_: Vec<BQF<Bignum>> = pks.iter().map(convert_bqf).collect();
            assert!(crate::verify_dealing::<
                Bignum,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(&convert_pp(&pp), &pks_, &d));
        }

        {
            let dealings: Vec<Dealing<G1Projective, Scalar, Bignum, _>> = dealings
                .iter()
                .map(|d| Dealing {
                    common_encryption: convert_bqf(&d.common_encryption),
                    encryptions: d.encryptions.iter().map(|e| convert_bqf(&e)).collect(),
                    correct_sharing_proof: Proof {
                        w: convert_bqf(&d.correct_sharing_proof.w),
                        x: d.correct_sharing_proof.x.clone(),
                        y: convert_bqf(&d.correct_sharing_proof.y),
                        z_r: int_to_bignum(&d.correct_sharing_proof.z_r),
                        z_s: d.correct_sharing_proof.z_s.clone(),
                    },
                    cmt: d.cmt.clone(),
                    p: d.p.clone(),
                })
                .collect();
            // for i in 0..=threshold {
            let DKGPublicResult {
                master_public_key,
                public_key_shares,
                public_poly: _,
            } = crate::aggregate_verified_dealings_public::<
                Bignum,
                G1Projective,
                blstrs::Scalar,
                Vec<blstrs::Scalar>,
            >(&convert_pp(&pp), &dealings);

            println!("mpk={:?},pks={:?}", master_public_key, public_key_shares);

            // g^a0 g^a1 g^a2 g^a3 g^a4
        }
    }

    fn convert_pp(p0: &PublicParameters<Integer>) -> PublicParameters<Bignum> {
        PublicParameters {
            n: p0.n,
            threshold: p0.threshold,
            security_level: SecurityLevel::SecLvl112,
            q: int_to_bignum(&p0.q),
            discriminant: int_to_bignum(&p0.discriminant),
            encryption_scheme: ClHSMqInstance {
                class_number_h_bound: int_to_bignum(&p0.encryption_scheme.class_number_h_bound),
                generator_f: convert_bqf(&p0.encryption_scheme.generator_f),
                generator_h: convert_bqf(&p0.encryption_scheme.generator_h),
                discriminant: int_to_bignum(&p0.discriminant),
                q: int_to_bignum(&p0.q),
            },
        }
    }*/
}
