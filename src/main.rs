//  Ideal class group cryptography
//  Copyright (C) 2024  Julien Coolen
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>

use class_group::{primitives::cl_dl_lcm::HSMCL, ABDeltaTriple, BinaryQF};
use curv::{
    arithmetic::{Converter, Samplable},
    BigInt,
};
use num_traits::One;

#[derive(Clone, Debug)]
struct NizkDlogProof {
    public_coin: BigInt,
    blinded_log: BigInt,
}

// Public coin (computed via Fiat-Shamir)
fn nizk_dlog_challenge(base: &BinaryQF, h: &BinaryQF, a: &BinaryQF) -> BigInt {
    let mut hasher = blake3::Hasher::new();
    hasher.update(base.to_bytes().as_slice());
    hasher.update(h.to_bytes().as_slice());
    hasher.update(a.to_bytes().as_slice());
    let hash = hasher.finalize();
    BigInt::from_bytes(hash.as_bytes())
}

// Schnorr's sigma protocol for proof of discrete logarithm,
// transformed into a non-interactive proof through Fiat-Shamir
fn nizk_dlog_prove(base: &BinaryQF, h: &BinaryQF, log: &BigInt, bound: &BigInt) -> NizkDlogProof {
    let r = BigInt::sample_range(&BigInt::one(), bound);
    let a = base.exp(&r);
    let c = nizk_dlog_challenge(base, h, &a);
    let s = r + log * &c;
    NizkDlogProof {
        public_coin: c,
        blinded_log: s,
    }
}

fn nizk_dlog_verify(proof: NizkDlogProof, base: &BinaryQF, h: &BinaryQF, bound: &BigInt) -> bool {
    /*if &proof.blinded_log >= &bound {
        return false;
    }*/
    let b = h.exp(&proof.public_coin).inverse();
    let a = base.exp(&proof.blinded_log).compose(&b).reduce();
    let c = nizk_dlog_challenge(base, h, &a);
    c.eq(&proof.public_coin)
}

fn main() {
    unsafe {
        class_group::pari_init(50000000, 2);
    }
    let a: BigInt = BigInt::from_str_radix("1347310664179468558147371727982960102805371574927252724399119343247182932538452304549609704350360058405827948976558722087559341859252338031258062288910984654814255199874816496621961922792890687089794760104660404141195904459619180668507135317125790028783030121033883873501532619563806411495141846196437", 10).unwrap();

    let b = BigInt::from(2);
    let delta = -BigInt::from(3) * BigInt::from(201);
    let abd = ABDeltaTriple { a, b, delta };
    let pf = BinaryQF::binary_quadratic_form_disc(&abd);
    let pari_qf = pf.qf_to_pari_qf();
    let pf2 = BinaryQF::pari_qf_to_qf(pari_qf);
    assert_eq!(pf, pf2);

    let q = BigInt::from_str_radix(
        "115792089210356248762697446949407573529996955224135760342422259061068512044369",
        10,
    )
    .unwrap();
    let hsmcl = HSMCL::keygen(&q, &1600);
    let u = hsmcl.pk.gq;
    let group = u;

    let bound = BigInt::from_str_radix("134731066417946855817371727982960102805371574927252724399119343247182932538452304549609704350360058405827948976558722087559341859252338031258062288910984654814255199874816496621961922792890687089794760104660404141195904459619180668507135317125790028783030121033883873501532619563806411495141846196437", 10).unwrap();
    let log = BigInt::sample_range(&BigInt::from_bytes(&[1u8]), &bound);

    let base = group.clone();
    let h = base.exp(&log);
    let proof = nizk_dlog_prove(&base, &h, &log, &bound);
    assert!(nizk_dlog_verify(proof, &base, &h, &bound))
}
