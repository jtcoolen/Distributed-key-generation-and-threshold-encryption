// TODO implement compact variant

use gmp_mpfr_sys::mpfr::sec;
use rand_core::CryptoRng;
use rug::rand::MutRandState;

use crate::bqf::{BinaryQuadraticForm, BQF};
use crate::z;

#[derive(Clone, Copy)]
pub enum SecurityLevel {
    SecLvl112,
    SecLvl128,
}

impl Into<u32> for SecurityLevel {
    fn into(self) -> u32 {
        match self {
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
    Z: crate::z::Z + std::fmt::Debug + Clone,
    Rng: CryptoRng,
    BQF: BinaryQuadraticForm<Z> + Clone,
>
{
    // setup
    fn new(q: Z, security_level: SecurityLevel, rng: &mut Rng) -> Self;

    // generates an asymmetric key pair (pk, sk)
    fn keygen(rng: &mut Rng) -> (BQF, Z);

    // ciphertext is composed of two binary quadratic forms
    fn encrypt(public_key: BQF, cleartext: Z, rng: &mut Rng) -> (BQF, BQF);

    fn decrypt(private_key: Z, ciphertext: (BQF, BQF));

    fn add_cleartexts(public_key: BQF, cleartext1: Z, cleartext2: Z) -> Z;

    fn scale_cleartext(public_key: BQF, cleartext: Z, scaling_factor: Z) -> Z;

    fn add_ciphertexts(
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF);

    fn scale_ciphertext(public_key: BQF, ciphertext: (BQF, BQF), scaling_factor: Z) -> (BQF, BQF);
}

// We'll take q=order of BLS12-381 scalar field F_r (multiplicative order)
// so that the messages matches scalar elements (F_r)^*
struct ClHSMqInstance<Z, BQF>
where
    BQF: BinaryQuadraticForm<Z> + Clone,
    Z: z::Z + std::fmt::Debug + std::clone::Clone,
{
    generator_F: BQF,
    generator_H: BQF,
    discriminant: Z,
    q: Z,
}

// TODO: note that with a security level of 112 we can support up to 3 multiplications
// of 1348-bit integers fitting inside 4096 bit ints

fn sample_random_p<Z: z::Z, R: CryptoRng + MutRandState>(
    q: &Z,
    q_bit_size: u32,
    disc_bit_size: u32,
    rng: &mut R,
) -> Z {
    let mut p = Z::sample_bits(disc_bit_size - q_bit_size, rng);
    let three = Z::from(3);
    let four = Z::from(4);
    let target_residue = if q.take_mod(&four).eq(&three) {
        Z::from(1)
    } else {
        three
    };
    while !p.take_mod(&four).eq(&target_residue) && q.kronecker(&p) != -1 {
        p = p.next_prime();
    }
    p
}

impl<Z, Rng, BQF> ClHSMq<Z, Rng, BQF> for ClHSMqInstance<Z, BQF>
where
    Z: z::Z + std::fmt::Debug + std::clone::Clone,
    Rng: CryptoRng + MutRandState,
    BQF: BinaryQuadraticForm<Z> + Clone,
{
    fn new(q: Z, security_level: SecurityLevel, rng: &mut Rng) -> Self {
        // "We denote η(λ) the bitsize of a fundamental discriminant ∆K such that"
        // "the computation of the class number h(∆K) and computation of discrete logarithms"
        // "in Cl(∆K ) takes at least 2^λ operations (see Table 2 for concrete sizes)."
        // lambda: security level
        // for lambda = 128, eta(lambda)=1827 bits
        // mu: q is a mu-bits prime
        // with BLS12-381 mu=254=log2(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        let q_bit_size = q.bit_size() as u32; // TODO add assert on bit size of q
        let disc_bit_size = discriminant_bit_size(security_level);
        assert!(q_bit_size >= <SecurityLevel as Into<u32>>::into(security_level));
        assert!(q_bit_size < disc_bit_size);
        let p = sample_random_p(&q, q_bit_size, disc_bit_size, rng);
        let fundamental_discriminant = p.mul(&q).neg();
        let discriminant_conductor_q = q.sqr().mul(&fundamental_discriminant);
        // G_hat commutative group of unknown order M * s_hat
        // G_hat isomorphic to F x H where F has order M
        let generator_H = compute_generator_H(&discriminant_conductor_q);
        let mut c = Z::from(1).sub(&fundamental_discriminant);
        c.divide_by_4();
        let generator_F = BQF::new(&q.sqr(), &q, &c);

        ClHSMqInstance {
            generator_F, // generator of subgroup of G_hat of order M
            generator_H, // generator of subgroup H of G_hat
            discriminant: fundamental_discriminant,
            q,
        }
    }

    fn keygen(rng: &mut Rng) -> (BQF, Z) {
        todo!()
    }

    fn encrypt(public_key: BQF, cleartext: Z, rng: &mut Rng) -> (BQF, BQF) {
        todo!()
    }

    fn decrypt(private_key: Z, ciphertext: (BQF, BQF)) {
        todo!()
    }

    fn add_cleartexts(public_key: BQF, cleartext1: Z, cleartext2: Z) -> Z {
        todo!()
    }

    fn scale_cleartext(public_key: BQF, cleartext: Z, scaling_factor: Z) -> Z {
        todo!()
    }

    fn add_ciphertexts(
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF) {
        todo!()
    }

    fn scale_ciphertext(public_key: BQF, ciphertext: (BQF, BQF), scaling_factor: Z) -> (BQF, BQF) {
        todo!()
    }
}

fn compute_generator_H<
    Z: std::clone::Clone + std::fmt::Debug + z::Z,
    BQF: BinaryQuadraticForm<Z> + Clone,
>(
    discriminant_conductor_q: &Z,
) -> BQF {
    let mut l = Z::from(2);
    while discriminant_conductor_q.kronecker(&l) != 1 {
        l = l.next_prime();
    }
    BQF::prime_form(discriminant_conductor_q, &l)
        .double()
        .pow(&discriminant_conductor_q)
        .reduce()
}

#[cfg(test)]
mod tests {}
