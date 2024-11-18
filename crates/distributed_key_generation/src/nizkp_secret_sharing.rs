// SPDX-FileCopyrightText: 2024 Nomadic Labs <contact@nomadic-labs.com>
//
// SPDX-License-Identifier: MIT

use std::fmt::Debug;
use std::marker::PhantomData;
use std::time::Instant;

#[cfg(feature = "random")]
use crate::z::Randomizable;
#[cfg(feature = "random")]
use rand_core::{CryptoRng, RngCore};

use serde::Serialize;

use crate::bqf::{BinaryQuadraticForm, BQF};
use crate::cl_hsmq::power_f;
use crate::elliptic_curve::EllipticCurve;
use crate::signed4096::Bignum;
use crate::z;
use crate::z::Z;

///! Implementation of Non-Interactive Zero-Knowledge (NIZK) proofs of correct Shamir Secret Sharing

/// Converts an instance of type `Z` to a `Bignum` by converting the integer value to bytes
/// and then constructing a `Bignum` from those bytes.
///
/// # Parameters
///
/// - `n`: A value of type `Z` to be converted.
///
/// # Returns
///
/// A `Bignum` representing the integer value of `n`.
fn int_to_bignum<Z: crate::z::Z>(n: &Z) -> Bignum {
    let (n, sgn) = n.to_bytes_be();
    Bignum::from_bytes_be(n, sgn)
}

/// Converts a `BinaryQuadraticForm` with coefficients of type `Z` to a `BinaryQuadraticForm`
/// with coefficients of type `Bignum`.
///
/// # Parameters
///
/// - `qf`: A `BinaryQuadraticForm` with coefficients of type `Z`.
///
/// # Returns
///
/// A `BinaryQuadraticForm` with coefficients converted to `Bignum`.
fn convert_bqf_z<Z: crate::z::Z + Clone + PartialEq + Debug>(qf: &BQF<Z>) -> BQF<Bignum> {
    BQF::new(
        &int_to_bignum(&qf.a()),
        &int_to_bignum(&qf.b()),
        &int_to_bignum(&qf.c()),
    )
}

/// Configuration for the cryptographic protocol, including parameters such as the degree,
/// threshold, and generators for the protocol.
///
/// # Type Parameters
///
/// - `Z`: The type representing integers, implementing `crate::z::Z`.
///
/// # Fields
///
/// - `n`: The number of participants in the protocol.
/// - `threshold`: The threshold of tolerated misbehaving participants
///   (< n / 2, or equivalently, n >= 2*t+1) for the success of the protocol.
/// - `q`: The order of the subgroup of order q.
/// - `discriminant`: The value of the fundamental discriminant used in the protocol.
/// - `generator_h`: A generator for the binary quadratic form `H` used in the protocol.
/// - `generator_f`: A generator for the binary quadratic form `F` used in the protocol.
/// - `_integer_type`: A phantom data marker for the integer type.
#[derive(Serialize)]
pub(crate) struct Config<Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
{
    n: u32,
    threshold: u32,
    q: Z,
    discriminant: Z,
    generator_h: BQF<Z>,
    generator_f: BQF<Z>,
    _integer_type: PhantomData<Z>,
}

impl<Z> Config<Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
{
    pub(crate) fn new<
        E: EllipticCurve<S = S> + std::clone::Clone,
        S: crate::scalar::Scalar + Clone,
    >(
        n: u32,
        threshold: u32,
        q: Z,
        discriminant: Z,
        generator_h: BQF<Z>,
        generator_f: BQF<Z>,
    ) -> Self {
        Config {
            n,
            threshold,
            q,
            discriminant,
            generator_h,
            generator_f,
            _integer_type: Default::default(),
        }
    }
}

/// Represents the witness for the corresponding secret sharing `Instance`,
/// containing the secret scalars and an integer `r`.
///
/// # Type Parameters
///
/// - `Z`: The type representing integers, implementing `crate::z::Z`.
/// - `S`: The type representing scalars, implementing `crate::scalar::Scalar`.
///
/// # Fields
///
/// - `s`: A vector of scalars representing the secret values.
/// - `r`: A "random coin".
pub(crate) struct Witness<Z, S>
where
    Z: crate::z::Z,
    S: crate::scalar::Scalar + Clone,
{
    pub(crate) s: Vec<S>,
    pub(crate) r: Z,
}

/// Serializes a vector of elliptic curve elements into bytes.
///
/// # Parameters
///
/// - `v`: A slice of elliptic curve elements to serialize.
/// - `serializer`: A serializer to use for the serialization process.
///
/// # Returns
///
/// A result containing the serialized bytes or an error if serialization fails.
///
/// # Type Parameters
///
/// - `S`: The type implementing `serde::Serializer`.
/// - `Z`: The type representing scalars used in elliptic curves.
/// - `E`: The elliptic curve type implementing `EllipticCurve<S = Z>`.
fn vec_ec_tobytes<S, E, Z>(v: &[E], serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    Z: crate::scalar::Scalar,
    E: EllipticCurve<S = Z>,
{
    let res = v.iter().flat_map(|e| e.to_bytes()).collect::<Vec<u8>>();

    serializer.serialize_bytes(&res)
}

/// Represents the instance for the protocol, containing public keys, ciphertexts, and other
/// relevant data.
///
/// # Type Parameters
///
/// - `Z`: The type representing integers, implementing `crate::z::Z`.
/// - `E`: The elliptic curve type implementing `EllipticCurve<S = S>`.
/// - `S`: The type representing scalars, implementing `crate::scalar::Scalar`.
///
/// # Fields
///
/// - `public_keys`: A vector of `BinaryQuadraticForm` instances representing public keys.
/// - `ciphertexts_common`: A `BinaryQuadraticForm` representing common ciphertexts.
/// - `ciphertexts`: A vector of `BinaryQuadraticForm` instances representing individual ciphertexts.
/// - `polynomial_coefficients_commitments`: A vector of elliptic curve elements representing
///   polynomial coefficients commitments.
#[derive(Serialize)]
pub struct Instance<Z, E, S>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    E: EllipticCurve<S = S> + Clone,
    S: crate::scalar::Scalar,
{
    pub(crate) public_keys: Vec<BQF<Z>>,
    pub(crate) ciphertexts_common: BQF<Z>,
    pub(crate) ciphertexts: Vec<BQF<Z>>,
    #[serde(serialize_with = "vec_ec_tobytes")]
    pub(crate) polynomial_coefficients_commitments: Vec<E>,
    pub(crate) _scalar_type: PhantomData<S>,
    pub(crate) _int_type: PhantomData<Z>,
}

/// Represents the proof for the protocol, containing values for verification.
///
/// # Type Parameters
///
/// - `Z`: The type representing integers, implementing `crate::z::Z`.
/// - `E`: The elliptic curve type implementing `EllipticCurve<S = S>`.
/// - `S`: The type representing scalars, implementing `crate::scalar::Scalar`.
///
/// # Fields
///
/// - `w`: A `BinaryQuadraticForm` used in the proof.
/// - `x`: An elliptic curve element used in the proof.
/// - `y`: A `BinaryQuadraticForm` used in the proof.
/// - `z_r`: An integer value used in the proof.
/// - `z_s`: A scalar value used in the proof.
#[derive(Clone)]
pub(crate) struct Proof<Z, E, S>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    E: EllipticCurve<S = S>,
    S: crate::scalar::Scalar,
{
    pub(crate) w: BQF<Z>,
    pub(crate) x: E,
    pub(crate) y: BQF<Z>,
    pub(crate) z_r: Z, // integer
    pub(crate) z_s: S, // scalar field
}

/// Converts an `Instance` to use `Bignum` types instead of the original integer type `Z`.
///
/// # Parameters
///
/// - `instance`: The `Instance` to be converted.
///
/// # Returns
///
/// A new `Instance` with all integer values converted to `Bignum`.
fn convert_instance_to_bignum<Z, E, S>(instance: &Instance<Z, E, S>) -> Instance<Bignum, E, S>
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
{
    Instance {
        public_keys: instance.public_keys.iter().map(convert_bqf_z).collect(),
        ciphertexts_common: convert_bqf_z(&instance.ciphertexts_common),
        ciphertexts: instance.ciphertexts.iter().map(convert_bqf_z).collect(),
        polynomial_coefficients_commitments: instance.polynomial_coefficients_commitments.clone(),
        _scalar_type: Default::default(),
        _int_type: Default::default(),
    }
}

/// Generates a proof for the protocol based on the provided configuration, instance, and witness.
///
/// # Parameters
///
/// - `config`: The configuration for the protocol.
/// - `instance`: The instance containing public keys, ciphertexts, and other data.
/// - `witness`: The witness containing secret values and an integer `r`.
/// - `rng`: A random number generator.
///
/// # Returns
///
/// A `Proof` instance containing values used for verification.
#[cfg(feature = "random")]
pub(crate) fn prove<
    R: CryptoRng + RngCore,
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
>(
    config: &Config<Z>,
    instance: &Instance<Z, E, S>,
    witness: &Witness<Z, S>,
    rng: &mut R,
) -> Proof<Z, E, S>
where
    Z: Randomizable,
{
    // Generate random values

    let alpha = S::random(rng);
    let rho = Z::sample_range(rng, &Z::zero(), &S::modulus_as_z());
    let w = config.generator_h.pow(&rho);

    // Calculate x
    let mut x = E::generator();
    x.mul_assign(&alpha);

    // Convert instance to Bignum
    let instance_bn = convert_instance_to_bignum(instance);

    // Compute gamma values
    let gamma_s = S::hash(bincode::serialize(&instance_bn).unwrap());
    let gamma_z = S::to_z(&gamma_s);

    // Compute powers of gamma_z
    let points = compute_gamma_powers(config.n as usize, &gamma_z);

    // Compute y
    let y = compute_y(config, instance, &points, &rho, &alpha);

    // Compute gamma_prime
    let gamma_prime_s = compute_gamma_prime_s(&gamma_z, &w, &x, &y);
    let gamma_prime_z = S::to_z(&gamma_prime_s);

    // Compute z_r and z_s
    let z_r = witness.r.mul(&gamma_prime_z).add(&rho);
    let z_s = compute_z_s::<Z, S, E>(witness, &gamma_s, &gamma_prime_s, &alpha);

    Proof { w, x, y, z_r, z_s }
}

/// Computes powers of a value `gamma_z` up to `n` times.
///
/// # Parameters
///
/// - `n`: The number of powers to compute.
/// - `gamma_z`: The base value to compute powers of.
///
/// # Returns
///
/// A vector containing the computed powers of `gamma_z`.
fn compute_gamma_powers<Z>(n: usize, gamma_z: &Z) -> Vec<Z>
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
{
    let mut gamma_pow = Z::from(1);
    (0..n)
        .map(|_| {
            gamma_pow = gamma_pow.mul(gamma_z);
            gamma_pow.clone()
        })
        .collect()
}

/// Computes the value `y` based on the provided configuration, instance, and gamma powers.
///
/// # Parameters
///
/// - `config`: The configuration for the protocol.
/// - `instance`: The instance containing public keys and ciphertexts.
/// - `points`: Powers of `gamma_z` used for computation.
/// - `rho`: A random integer value used in the computation.
/// - `alpha`: A random scalar value used in the computation.
///
/// # Returns
///
/// The computed `BinaryQuadraticForm` value for `y`.
fn compute_y<Z, E, S>(
    config: &Config<Z>,
    instance: &Instance<Z, E, S>,
    points: &[Z],
    rho: &Z,
    alpha: &S,
) -> BQF<Z>
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
{
    let y_base = BQF::multiexp(&instance.public_keys, points);
    let y_exponent = power_f(
        &config.generator_f,
        &config.discriminant,
        &config.q,
        &S::to_z(alpha),
    );
    y_base.pow(rho).compose(&y_exponent)
}

/// Computes the scalar value `gamma_prime_s` used in the proof.
///
/// # Parameters
///
/// - `gamma_z`: The gamma value used in the computation.
/// - `w`: The `BinaryQuadraticForm` value `w` from the proof.
/// - `x`: The elliptic curve element `x` from the proof.
/// - `y`: The `BinaryQuadraticForm` value `y` from the proof.
///
/// # Returns
///
/// The computed scalar `gamma_prime_s`.
fn compute_gamma_prime_s<Z, E, S>(gamma_z: &Z, w: &BQF<Z>, x: &E, y: &BQF<Z>) -> S
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
{
    let mut b: Vec<u8> = vec![];
    b.append(&mut int_to_bignum(gamma_z).to_bytes_be().0);
    b.append(&mut bincode::serialize(&convert_bqf_z(w)).unwrap());
    b.append(&mut x.to_bytes());
    b.append(&mut bincode::serialize(&convert_bqf_z(y)).unwrap());
    S::hash(b)
}

/// Computes the scalar value `z_s` based on the witness and other parameters.
///
/// # Parameters
///
/// - `witness`: The witness containing secret scalars and integer `r`.
/// - `gamma_s`: The scalar value `gamma_s` computed from the instance.
/// - `gamma_prime_s`: The scalar value `gamma_prime_s` computed from the proof.
/// - `alpha`: The random scalar value used in the proof.
///
/// # Returns
///
/// The computed scalar `z_s`.
fn compute_z_s<Z, S, E>(witness: &Witness<Z, S>, gamma_s: &S, gamma_prime_s: &S, alpha: &S) -> S
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
{
    let mut gamma_pow = gamma_s.clone();
    let mut z_s = witness
        .s
        .clone()
        .into_iter()
        .fold(S::zero(), |mut acc, mut s| {
            s.mul_assign(&gamma_pow);
            gamma_pow.mul_assign(gamma_s);
            acc.add_assign(&s);
            acc
        });
    z_s.mul_assign(gamma_prime_s);
    z_s.add_assign(alpha);
    z_s
}

/// Verifies the given proof against the provided configuration and instance.
///
/// # Parameters
///
/// - `config`: The configuration for the protocol.
/// - `instance`: The instance containing public keys, ciphertexts, and other data.
/// - `proof`: The proof to be verified.
///
/// # Returns
///
/// `true` if the proof is valid, `false` otherwise (with overwhelming probability).
pub(crate) fn verify<
    Z: crate::z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
>(
    config: &Config<Z>,
    instance: &Instance<Z, E, S>,
    proof: &Proof<Z, E, S>,
) -> bool {
    let now = Instant::now();
    // Convert instance to Bignum
    let instance_bn = convert_instance_to_bignum(instance);

    // Compute gamma values
    let gamma_s = S::hash(bincode::serialize(&instance_bn).unwrap());
    let gamma_z = S::to_z(&gamma_s);

    // Compute gamma_prime
    let gamma_prime_s = compute_gamma_prime_s(&gamma_z, &proof.w, &proof.x, &proof.y);
    let gamma_prime_z = S::to_z(&gamma_prime_s);
    println!("challenge {:?}", now.elapsed());

    let now = Instant::now();
    // Verify first equation: w * instance.ciphertexts_common^gamma_prime_z == generator_h^z_r
    let lhs = proof
        .w
        .compose(&instance.ciphertexts_common.pow(&gamma_prime_z));
    let rhs = config.generator_h.pow(&proof.z_r);
    if !lhs.equals(&rhs) {
        return false;
    }
    println!("first check {:?}", now.elapsed());

    let now = Instant::now();
    // Verify second equation: multiexp(instance.polynomial_coefficients_commitments, scalars) * gamma_prime_s + x == generator * z_s
    let scalars = compute_scalars(config, instance, &gamma_s);
    let mut lhs = E::multiexp(&instance.polynomial_coefficients_commitments, &scalars);
    lhs.mul_assign(&gamma_prime_s);
    lhs.add_assign(&proof.x);

    let mut rhs = E::generator();
    rhs.mul_assign(&proof.z_s);
    if !lhs.equals(&rhs) {
        return false;
    }
    println!("second check {:?}", now.elapsed());

    let now = Instant::now();
    // Verify third equation: multiexp(instance.ciphertexts, points)^gamma_prime_z * proof.y == multiexp(instance.public_keys, points)^proof.z_r * f^z_s
    let points = compute_gamma_powers(config.n as usize, &gamma_z);
    
    let lhs = BQF::multiexp(&instance.ciphertexts, &points)
        .pow(&gamma_prime_z)
        .compose(&proof.y.reduce());

        
    let mut rhs = BQF::multiexp(&instance.public_keys, &points);

    let f_pow = power_f(
        &config.generator_f,
        &config.discriminant,
        &config.q,
        &S::to_z(&proof.z_s),
    );
    rhs = rhs.pow(&proof.z_r).compose(&f_pow);

    let res = lhs.equals(&rhs);
  
    println!("third check {:?}", now.elapsed());
    
    res
}

fn compute_scalars<Z, S, E>(config: &Config<Z>, instance: &Instance<Z, E, S>, gamma_s: &S) -> Vec<S>
where
    Z: z::Z + std::fmt::Debug + Clone + Serialize + PartialEq,
    E: EllipticCurve<S = S> + std::clone::Clone,
    S: crate::scalar::Scalar + Clone + Debug,
{
    instance
        .polynomial_coefficients_commitments
        .clone()
        .into_iter()
        .enumerate()
        .map(|(j, _)| {
            (1..=config.n).fold(S::zero(), |mut acc, i| {
                let mut pow_gamma = gamma_s.pow(i as u64);
                pow_gamma.mul_assign(&S::from(i.pow(j as u32) as u64));
                acc.add_assign(&pow_gamma);
                acc
            })
        })
        .collect()
}
