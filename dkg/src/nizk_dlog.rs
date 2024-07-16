// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use rand_core::OsRng;
use std::cmp::Ordering::Greater;

use crate::bqf::{BinaryQuadraticForm, BQF};

#[derive(Clone, Debug)]
pub(crate) struct NizkDlogProof<Z: crate::z::Z> {
    public_coin: Z,
    blinded_log: Z,
}

// Public coin (computed via Fiat-Shamir)
fn nizk_dlog_challenge<Z: crate::z::Z + std::fmt::Debug + std::clone::Clone>(
    base: &BQF<Z>,
    h: &BQF<Z>,
    a: &BQF<Z>,
) -> Z {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&base.to_bytes());
    hasher.update(&h.to_bytes());
    hasher.update(&a.to_bytes());
    let hash = hasher.finalize();
    Z::from_bytes_be(hash.as_bytes().to_vec())
}

// Schnorr's sigma protocol for proof of discrete logarithm,
// transformed into a non-interactive proof through Fiat-Shamir
pub fn nizk_dlog_prove<Z: crate::z::Z + std::fmt::Debug + std::clone::Clone>(
    base: &BQF<Z>,
    h: &BQF<Z>,
    log: &Z,
    bound: &Z,
) -> NizkDlogProof<Z> {
    let r = Z::from(1); //Z::sample_range(&mut OsRng, &Z::from(1), bound);
    let a = base.pow(&r).reduce();
    let c = nizk_dlog_challenge(base, h, &a);
    let s = r.add(&log.mul(&c));
    println!("a1 disc ={:?}", a.discriminant());
    println!("log={:?}, c={:?}", log, c);
    println!("base1={:?}, h1={:?}, a1={:?}", base, h, a);
    NizkDlogProof {
        public_coin: c,
        blinded_log: s,
    }
}

pub fn nizk_dlog_verify<Z: crate::z::Z + std::fmt::Debug + std::clone::Clone>(
    proof: &NizkDlogProof<Z>,
    base: &BQF<Z>,
    h: &BQF<Z>,
    bound: &Z,
) -> bool {
    /*if proof.blinded_log.compare(bound) == Greater {
        return false;
    }`*/
    let b = h.pow(&proof.public_coin).reduce().inverse().reduce();
    let a = base.pow(&proof.blinded_log).reduce().compose(&b).reduce();
    println!("a disc ={:?}", a.discriminant());
    let c = nizk_dlog_challenge(base, h, &a);
    println!(
        "blinded log={:?}, pc={:?}",
        proof.blinded_log, proof.public_coin
    );
    println!("base={:?}, h={:?}, a={:?}", base, h, a);
    println!("c={:?}, pc={:?}", c, proof.public_coin);
    c.eq(&proof.public_coin)
}

#[cfg(test)]
mod tests {
    use crate::bqf::BinaryQuadraticForm;
    use crate::cl_hsmq::{ClHSMq, ClHSMqInstance};
    use crate::nizk_dlog::nizk_dlog_verify;
    use crate::vss::NIVSS::{Config, Dealing, PublicParameters};
    use crate::z::Z;
    use blstrs::{G1Projective, Scalar};
    use rand_core::OsRng;
    use rug::Integer;

    #[test]
    fn test_poe() {
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
            index: 0,
        };
        let (public_key, secret_key) = encryption_scheme.keygen(&mut OsRng);

        let base = encryption_scheme.generator_h();

        assert!(base.pow(&secret_key).equals(&public_key));

        let c = Integer::random(&mut OsRng);
        let blinded_log = <Integer as Z>::from(1u64).add(&secret_key.mul(&c));
        let b = base
            .pow(&secret_key)
            .reduce()
            .pow(&c)
            .reduce()
            .inverse()
            .reduce();
        let a = base.pow(&blinded_log).reduce().compose(&b).reduce();
        println!("base={:?}\na   ={:?}", base, a);
        assert!(a.equals(&base));

        let bound = encryption_scheme.class_number_bound_h(); // TODO fix bound

        let proof = crate::nizk_dlog::nizk_dlog_prove(&base, &public_key, &secret_key, &bound);

        let v = nizk_dlog_verify(&proof, &base, &public_key, &bound);
        assert!(v)
    }
}
