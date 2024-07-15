use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::Add;

use ff::{Field, PrimeField};
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

    fn from_bytes(_: Vec<u8>) -> Self;

    fn sqr(&self) -> Self;

    fn zero() -> Self;
    fn from(n: u64) -> Self;
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
        Z::from_bytes(self.to_bytes_be().to_vec())
    }

    fn from_bytes(b: Vec<u8>) -> Self {
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
}

/*fn hash() {
    let mut hasher = blake3::Hasher::new();
    hasher.update(base.to_bytes().as_slice());
    hasher.update(h.to_bytes().as_slice());
    hasher.update(a.to_bytes().as_slice());
    let hash = hasher.finalize();
    BigInt::from_bytes(hash.as_bytes())
}*/

pub trait NIZKProofOfCorrectSharing<Z, BQF, E, S>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z> + Clone,
    E: EllipticCurve<S>,
    S: Scalar,
{
    type Instance;
    type Proof;
    type Witness;

    /// Creates a new instance with the given threshold.
    /// Since the threshold is small, FFT is not required.
    fn new<R: CryptoRng + RngCore>(
        n: u32,
        threshold: u32,
        generator_H: BQF,
        generator_F: BQF,
        rng: &mut R,
    ) -> Self;

    /// Generates a proof for the given instance and witness.
    fn prove<R: CryptoRng + RngCore>(
        &self,
        instance: &Self::Instance,
        witness: &Self::Witness,
        rng: &mut R,
    ) -> Self::Proof;

    /// Verifies the given proof against the provided instance.
    fn verify(&self, instance: &Self::Instance, proof: &Self::Proof) -> bool;
}

struct Witness<Z, S>
where
    Z: crate::z::Z,
    S: Scalar + Clone,
{
    s: Vec<S>,
    r: Z,
}

fn vec_z_tobytes<S, Z>(v: &Vec<Z>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    Z: z::Z,
{
    let res = v.iter().flat_map(|e| e.to_bytes()).collect::<Vec<u8>>();

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

fn bqf_tobytes<S, BQF, Z>(v: &BQF, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    Z: z::Z + std::fmt::Debug + std::clone::Clone,
    BQF: BinaryQuadraticForm<Z>,
{
    let res = v.to_bytes();

    serializer.serialize_bytes(&res)
}

#[derive(Serialize)]
pub struct Instance<Z, BQF, E, S>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z> + std::clone::Clone + std::fmt::Debug,
    E: EllipticCurve<S> + Clone,
    S: Scalar,
{
    public_keys: Vec<BQF>,
    ciphertexts_common: BQF,
    ciphertexts: Vec<BQF>,
    #[serde(serialize_with = "vec_ec_tobytes")]
    polynomial_coefficients_commitments: Vec<E>,
    _scalar_type: PhantomData<S>,
    _int_type: PhantomData<Z>,
}

struct Proof<Z, BQF, E, S>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z>,
    E: EllipticCurve<S>,
    S: Scalar,
{
    w: BQF,
    x: E,
    y: BQF,
    z_r: Z, // integer
    z_s: S, // scalar field
}

struct PoCS<BQF, Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z>,
{
    n: u32,
    threshold: u32,
    generator_H: BQF,
    generator_F: BQF,
    _integer_type: PhantomData<Z>,
}

impl<Z, BQF, E, S> NIZKProofOfCorrectSharing<Z, BQF, E, S> for PoCS<BQF, Z>
where
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z> + Clone + std::fmt::Debug + Serialize,
    E: EllipticCurve<S, Scalar = S> + std::clone::Clone,
    S: Scalar + Clone,
{
    type Instance = Instance<Z, BQF, E, S>;
    type Proof = Proof<Z, BQF, E, S>;
    type Witness = Witness<Z, S>;

    fn new<R: CryptoRng + RngCore>(
        n: u32,
        threshold: u32,
        generator_H: BQF,
        generator_F: BQF,
        rng: &mut R,
    ) -> Self {
        todo!()
    }

    fn prove<R: CryptoRng + RngCore>(
        &self,
        instance: &Self::Instance,
        witness: &Self::Witness,
        rng: &mut R,
    ) -> Self::Proof {
        let alpha = S::random(rng);
        let rho = Z::sample_range(rng, &Z::zero(), &S::modulus_as_z()); // TODO correct range
        let w = self.generator_H.pow(&rho);

        let mut x = E::generator();
        x.mul_assign(&alpha);

        let mut hasher = blake3::Hasher::new();
        hasher.update(bincode::serialize(&instance).unwrap().as_slice());
        let gamma = hasher.finalize().as_bytes().to_vec();
        let gamma_z = Z::from_bytes(gamma.clone());
        let gamma_s = S::from_bytes(gamma);
        let mut gamma_pow = gamma_z.clone();

        let mut y: BQF = instance
            .public_keys
            .clone()
            .into_iter()
            .reduce(|acc, e| {
                let res = acc.compose(&e.pow(&gamma_pow));
                gamma_pow.sqr();
                res
            })
            .unwrap()
            .pow(&rho);
        y = y.compose(&self.generator_F.pow(&S::to_z(&alpha))); // TODO: One can use faster power_of_f here

        let mut hasher = blake3::Hasher::new();
        hasher.update(gamma_z.to_bytes().as_slice());
        hasher.update(bincode::serialize(&w).unwrap().as_slice());
        hasher.update(&x.to_bytes());
        hasher.update(bincode::serialize(&y).unwrap().as_slice());
        let gamma_prime = Z::from_bytes(hasher.finalize().as_bytes().to_vec());
        let z_r = witness.r.mul(&gamma_prime).add(&rho);
        let mut gamma_pow = gamma_s.clone();
        let z_s = witness
            .s
            .clone()
            .into_iter()
            .reduce(|acc, mut s| {
                s.mul_assign(&gamma_pow);
                gamma_pow.sqr();
                s
            })
            .unwrap();

        Proof { w, x, y, z_r, z_s }
    }

    fn verify(&self, instance: &Self::Instance, proof: &Self::Proof) -> bool {
        let mut hasher = blake3::Hasher::new();
        hasher.update(bincode::serialize(&instance).unwrap().as_slice());
        let mut gamma = hasher.finalize().as_bytes().to_vec();
        let gamma_z = Z::from_bytes(gamma.clone());
        let gamma_s = S::from_bytes(gamma.clone());

        // TODO refactor
        let mut hasher = blake3::Hasher::new();
        hasher.update(gamma_z.to_bytes().as_slice());
        hasher.update(bincode::serialize(&proof.w).unwrap().as_slice());
        hasher.update(&proof.x.to_bytes());
        hasher.update(bincode::serialize(&proof.y).unwrap().as_slice());
        let gamma_prime = hasher.finalize().as_bytes().to_vec();
        let gamma_prime_z = Z::from_bytes(gamma_prime.clone());
        let gamma_prime_s = S::from_bytes(gamma_prime);

        let lhs = proof
            .w
            .compose(&instance.ciphertexts_common.pow(&gamma_prime_z));
        let rhs = self.generator_H.pow(&proof.z_r);
        if !lhs.equals(&rhs) {
            return false;
        }

        let (_, mut lhs) = instance
            .polynomial_coefficients_commitments
            .clone()
            .into_iter()
            .enumerate()
            .reduce(|(_, acc), (j, mut a_j)| {
                let pow_gamma = gamma_s.pow(j as u64);
                let exp = (1..=self.n).fold(S::zero(), |mut acc, i| {
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
                gamma_pow = gamma_pow.sqr();
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
                gamma_pow = gamma_pow.sqr();
                (i, acc.compose(&res))
            })
            .unwrap();
        rhs = rhs
            .pow(&proof.z_r)
            .compose(&self.generator_F.pow(&S::to_z(&proof.z_s))); // TODO: One can use faster power_of_f here
        if !lhs.equals(&rhs) {
            return false;
        }
        true
    }
}

pub trait VerifiableSecretSharingPrimitives<LinearHomomorphicEncryption, Z, BQF>
where
    LinearHomomorphicEncryption: ClHSMq<Z, BQF>,
    Z: crate::z::Z + std::fmt::Debug + Clone,
    BQF: BinaryQuadraticForm<Z> + Clone,
{
    type PublicParams;
    type Instance;
    type ProofOfCorrectSharing;
    type Share;
    type Commitment;
    type Ciphertext;
    type PublicKey;
    type SecretKey;

    fn share(pp: &Self::PublicParams, s: Z) -> (Vec<Self::Share>, Self::Commitment);

    fn encrypt(
        pp: &Self::PublicParams,
        cmt: &Self::Commitment,
        shares: Vec<Self::Share>,
        pub_keys: Vec<Self::PublicKey>,
    ) -> (
        Self::Ciphertext,
        Vec<Self::Ciphertext>,
        Self::ProofOfCorrectSharing,
    );

    fn verify(
        pp: &Self::PublicParams,
        cmt: &Self::Commitment,
        ciphertext: &Self::Ciphertext,
        enc_shares: Vec<Self::Ciphertext>,
        pub_keys: Vec<Self::PublicKey>,
        proof: &Self::ProofOfCorrectSharing,
    ) -> bool;

    fn decrypt(
        pp: &Self::PublicParams,
        sk: &Self::SecretKey,
        enc_share: &Self::Ciphertext,
    ) -> Self::Share;

    fn verify_commitment(cmt: &Self::Commitment, i: usize, share: &Self::Share) -> bool;
}
