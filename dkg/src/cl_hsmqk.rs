// TODO implement compact variant

use crate::bqf::BinaryQuadraticForm;
use rand_core::CryptoRng;

/// Linear homomorphic encryption based on the hidden subgroup membership problem
/// with the CL framework
pub trait ClHSMqk<
    Z: crate::z::Z + std::fmt::Debug + Clone,
    Rng: CryptoRng,
    BQF: BinaryQuadraticForm<Z>,
>
{
    // setup
    fn new(q: Z, k: u64, security_level: u64, rng: &mut Rng) -> Self;

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

#[cfg(test)]
mod tests {}
