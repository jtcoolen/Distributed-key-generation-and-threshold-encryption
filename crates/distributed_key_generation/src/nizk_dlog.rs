// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use crate::bqf::{BinaryQuadraticForm, BQF};
#[cfg(feature = "random")]
use crate::z::Randomizable;

#[derive(Clone, Debug)]
pub struct NizkDlogProof<Z: crate::z::Z> {
    public_coin: Z,
    blinded_log: Z,
}

pub const LAMBDA_COMP_BITS: usize = 255; // computational security parameter
pub const LAMBDA_STAT_BITS: usize = 40; // statistical security parameter

// Public coin (computed via Fiat-Shamir)
fn nizk_dlog_challenge<Z: crate::z::Z + std::fmt::Debug + std::clone::Clone + PartialEq>(
    base: &BQF<Z>,
    h: &BQF<Z>,
    a: &BQF<Z>,
) -> Z {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&base.to_bytes());
    hasher.update(&h.to_bytes());
    hasher.update(&a.to_bytes());
    let hash = hasher.finalize();
    Z::from_bytes_be(hash.as_bytes().to_vec(), true)
}

// Schnorr's sigma protocol for proof of discrete logarithm,
// transformed into a non-interactive proof through Fiat-Shamir
pub fn nizk_dlog_prove<
    Z: crate::z::Z + std::fmt::Debug + std::clone::Clone + PartialEq,
    Rng: rand_core::CryptoRng + rand_core::RngCore,
>(
    base: &BQF<Z>,
    h: &BQF<Z>,
    log: &Z,
    bound: &Z,
    rng: &mut Rng,
) -> NizkDlogProof<Z>
where
    Z: Randomizable,
{
    let r = Z::sample_bits(
        (bound.bit_size() as usize + LAMBDA_COMP_BITS + LAMBDA_STAT_BITS) as u32,
        rng,
    );
    let a = base.pow(&r);
    let c = nizk_dlog_challenge(base, h, &a);
    let s = r.add(&log.mul(&c));
    NizkDlogProof {
        public_coin: c,
        blinded_log: s,
    }
}

pub fn nizk_dlog_verify<Z: crate::z::Z + std::fmt::Debug + std::clone::Clone + PartialEq>(
    proof: &NizkDlogProof<Z>,
    base: &BQF<Z>,
    h: &BQF<Z>,
    bound: &Z,
) -> bool {
    let bound = bound
        .mul(&Z::from(1).shl(LAMBDA_COMP_BITS as u32))
        .mul(&Z::from(1).shl((LAMBDA_STAT_BITS + 1usize) as u32));
    if proof.blinded_log.compare(&bound) == std::cmp::Ordering::Greater {
        return false;
    }
    let b = h.pow(&proof.public_coin).inverse();
    let a = base.pow(&proof.blinded_log).compose(&b);
    let c = nizk_dlog_challenge(base, h, &a);
    c == proof.public_coin
}

#[cfg(test)]
#[cfg(feature = "random")]
#[cfg(feature = "gmp")]
mod tests {
    use rand_core::OsRng;
    use rug::Integer;

    use crate::bqf::BinaryQuadraticForm;
    use crate::cl_hsmq::{ClHSMq, ClHSMqInstance};
    use crate::nizk_dlog::nizk_dlog_verify;

    #[test]
    fn test_poe() {
        let _number_of_participants: usize = 10;
        let _threshold: usize = 4; // < number_of_participants / 2
        let security_level = crate::cl_hsmq::SecurityLevel::SecLvl112;

        let encryption_scheme =
            ClHSMqInstance::new::<OsRng, blstrs::Scalar>(security_level, &mut OsRng);

        let (public_key, secret_key) = encryption_scheme.keygen::<OsRng>(&mut OsRng);

        let base: crate::bqf::BQF<Integer> = encryption_scheme.generator_h();

        assert!(base.pow(&secret_key).equals(&public_key));

        let bound = encryption_scheme.class_number_bound_h();

        let proof = crate::nizk_dlog::nizk_dlog_prove::<rug::Integer, OsRng>(
            &base,
            &public_key,
            &secret_key,
            &bound,
            &mut OsRng,
        );

        let v = nizk_dlog_verify(&proof, &base, &public_key, &bound);
        assert!(v)
    }
}
