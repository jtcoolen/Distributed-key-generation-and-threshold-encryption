use crate::bqf::BinaryQuadraticForm;
use rand::rngs::OsRng;

#[derive(Clone, Debug)]
struct NizkDlogProof<Z: crate::z::Z> {
    public_coin: Z,
    blinded_log: Z,
}

// Public coin (computed via Fiat-Shamir)
fn nizk_dlog_challenge<
    Z: crate::z::Z + std::fmt::Debug + std::clone::Clone,
    BQF: BinaryQuadraticForm<Z>,
>(
    base: &BQF,
    h: &BQF,
    a: &BQF,
) -> Z {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&base.to_bytes());
    hasher.update(&h.to_bytes());
    hasher.update(&a.to_bytes());
    let hash = hasher.finalize();
    Z::from_bytes(hash.as_bytes().to_vec())
}

// Schnorr's sigma protocol for proof of discrete logarithm,
// transformed into a non-interactive proof through Fiat-Shamir
fn nizk_dlog_prove<
    Z: crate::z::Z + std::fmt::Debug + std::clone::Clone,
    BQF: BinaryQuadraticForm<Z> + std::clone::Clone,
>(
    base: &BQF,
    h: &BQF,
    log: &Z,
    bound: &Z,
) -> NizkDlogProof<Z> {
    let r = Z::sample_range(&mut OsRng, &Z::from(1), bound);
    let a = base.pow(&r);
    let c = nizk_dlog_challenge(base, h, &a);
    let s = r.add(&log.mul(&c));
    NizkDlogProof {
        public_coin: c,
        blinded_log: s,
    }
}

fn nizk_dlog_verify<
    Z: crate::z::Z + std::fmt::Debug + std::clone::Clone,
    BQF: BinaryQuadraticForm<Z> + std::clone::Clone,
>(
    proof: NizkDlogProof<Z>,
    base: &BQF,
    h: &BQF,
    bound: &Z,
) -> bool {
    /*if &proof.blinded_log >= &bound {
        return false;
    }*/
    let b = h.pow(&proof.public_coin).inverse();
    let a = base.pow(&proof.blinded_log).compose(&b).reduce();
    let c = nizk_dlog_challenge(base, h, &a);
    c.eq(&proof.public_coin)
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_nizk_dlog() {
        /*let hsmcl = HSMCL::keygen(&q, &1600);
        let u = hsmcl.pk.gq;
        let group = u;
        println!("pf={:?}, disc={:?}", group, group.discriminant());

        let bound = BigInt::from_str_radix("134731066417946855817371727982960102805371574927252724399119343247182932538452304549609704350360058405827948976558722087559341859252338031258062288910984654814255199874816496621961922792890687089794760104660404141195904459619180668507135317125790028783030121033883873501532619563806411495141846196437", 10).unwrap();
        let log = BigInt::sample_range(&BigInt::from_bytes(&[1u8]), &bound);

        let base = group.clone();
        let h = base.exp(&log);
        let proof = nizk_dlog_prove(&base, &h, &log, &bound);*/
    }
}
