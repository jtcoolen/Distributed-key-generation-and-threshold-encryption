// TODO implement compact variant

use std::cmp::Ordering::Less;
use std::fmt::Debug;

use rand_core::{CryptoRng, OsRng};
use rug::Integer;
use thiserror::Error;

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
pub trait ClHSMq<Z: crate::z::Z + std::fmt::Debug + Clone, BQF: BinaryQuadraticForm<Z> + Clone> {
    // setup
    fn new<Rng: CryptoRng + rand_core::RngCore>(
        q: Z,
        security_level: SecurityLevel,
        rng: &mut Rng,
    ) -> Self;

    // generates an asymmetric key pair (pk, sk)
    fn keygen<Rng: CryptoRng + rand_core::RngCore>(&self, rng: &mut Rng) -> (BQF, Z);

    // ciphertext is composed of two binary quadratic forms
    fn encrypt<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_key: &BQF,
        cleartext: &Z,
        rng: &mut Rng,
    ) -> (BQF, BQF);

    fn encrypt_batch<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_keys: &Vec<BQF>,
        cleartexts: &Vec<Z>,
        rng: &mut Rng,
    ) -> (BQF, Vec<BQF>);

    fn decrypt(&self, private_key: &Z, ciphertext: &(BQF, BQF)) -> Z;

    fn add_cleartexts(&self, public_key: BQF, cleartext1: Z, cleartext2: Z) -> Z;

    fn scale_cleartext(&self, public_key: BQF, cleartext: Z, scaling_factor: Z) -> Z;

    fn add_ciphertexts<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF);

    fn scale_ciphertext<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext: (BQF, BQF),
        scaling_factor: Z,
    ) -> (BQF, BQF);

    fn generator_f(&self) -> BQF;

    fn generator_h(&self) -> BQF;
}

// We'll take q=order of BLS12-381 scalar field F_r (multiplicative order)
// so that the messages matches scalar elements (F_r)^*
pub struct ClHSMqInstance<Z, BQF>
where
    BQF: BinaryQuadraticForm<Z> + Clone,
    Z: z::Z + std::fmt::Debug + std::clone::Clone,
{
    class_number_H_bound: Z,
    generator_F: BQF,
    generator_H: BQF,
    discriminant: Z,
    q: Z,
}

// TODO: note that with a security level of 112 we can support up to 3 multiplications
// of 1348-bit integers fitting inside 4096 bit ints

fn sample_random_p<Z: z::Z, R: CryptoRng + rand_core::RngCore>(
    q: &Z,
    q_bit_size: u32,
    disc_bit_size: u32,
    rng: &mut R,
) -> Z {
    println!(
        "disc_bit_size - q_bit_size ={}, qbitsize={}",
        disc_bit_size, q_bit_size
    );
    assert!(disc_bit_size > q_bit_size);
    let mut p = Z::sample_bits(disc_bit_size - q_bit_size, rng);
    let three = Z::from(3);
    let four = Z::from(4);
    let target_residue = if q.take_mod(&four).eq(&three) {
        Z::from(1)
    } else {
        three
    };
    p = p.next_prime();
    while !p.take_mod(&four).eq(&target_residue) || q.kronecker(&p) != -1 {
        p = p.next_prime();
    }
    p
}

impl<Z, BQF> ClHSMq<Z, BQF> for ClHSMqInstance<Z, BQF>
where
    Z: z::Z + std::fmt::Debug + std::clone::Clone,
    BQF: BinaryQuadraticForm<Z> + Clone + std::fmt::Debug,
{
    fn new<Rng: CryptoRng + rand_core::RngCore>(
        q: Z,
        security_level: SecurityLevel,
        rng: &mut Rng,
    ) -> Self {
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
        println!("before generator H");
        println!(
            "disc funda {:?}, disc conduct {:?}, p={:?}",
            fundamental_discriminant, discriminant_conductor_q, p
        );
        let generator_H = compute_generator_H(&q, &discriminant_conductor_q);
        println!("after generator H");
        let mut c = Z::from(1).sub(&fundamental_discriminant);
        c.divide_by_4();
        let generator_F = BQF::new(&q.sqr(), &q, &c).reduce();
        let class_number_H_bound = BQF::class_number_bound(rng, &fundamental_discriminant);

        let lhs = q.sqr().mul(&Z::from(4));
        let rhs = Z::from(1).sub(&fundamental_discriminant);
        println!("lhs = {:?}, rhs = {:?}", lhs, rhs);
        // we're not in the large message variant case
        assert_eq!(lhs.compare(&rhs), Less); // TODO it should be less or equal
        ClHSMqInstance {
            class_number_H_bound,
            generator_F, // generator of subgroup of G_hat of order M
            generator_H, // generator of subgroup H of G_hat
            discriminant: fundamental_discriminant,
            q,
        }
    }

    // TODO typed secret/public keys
    fn keygen<Rng: CryptoRng + rand_core::RngCore>(&self, rng: &mut Rng) -> (BQF, Z) {
        let secret_key = Z::sample_range(rng, &Z::from(0), &self.class_number_H_bound);
        let public_key = self.generator_H.pow(&secret_key).reduce();
        (public_key, secret_key)
    }

    fn encrypt<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_key: &BQF,
        cleartext: &Z,
        rng: &mut Rng,
    ) -> (BQF, BQF) {
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_H_bound);
        let c1 = self.generator_H.pow(&r).reduce();
        let f_m = power_f(&self.generator_F, &self.discriminant, &self.q, &cleartext);
        let pk_r = public_key.pow(&r).reduce();
        let c2 = f_m.compose(&pk_r).reduce();
        (c1, c2)
    }

    fn encrypt_batch<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        public_keys: &Vec<BQF>,
        cleartexts: &Vec<Z>,
        rng: &mut Rng,
    ) -> (BQF, Vec<BQF>) {
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_H_bound);
        let c1 = self.generator_H.pow(&r).reduce();

        let c2s = cleartexts
            .iter()
            .zip(public_keys)
            .map(|(c, k)| {
                power_f(&self.generator_F, &self.discriminant, &self.q, &c)
                    .compose(&k.pow(&r).reduce())
                    .reduce()
            })
            .collect();

        (c1, c2s)
    }

    fn decrypt(&self, private_key: &Z, ciphertext: &(BQF, BQF)) -> Z {
        let dlog = ciphertext
            .1
            .reduce()
            .compose(&ciphertext.0.reduce().pow(&private_key).inverse().reduce())
            .reduce();
        self.dlog_solve_F(&dlog)
    }

    fn add_cleartexts(&self, public_key: BQF, cleartext1: Z, cleartext2: Z) -> Z {
        todo!()
    }

    fn scale_cleartext(&self, public_key: BQF, cleartext: Z, scaling_factor: Z) -> Z {
        todo!()
    }

    fn add_ciphertexts<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext1: (BQF, BQF),
        ciphertext2: (BQF, BQF),
    ) -> (BQF, BQF) {
        let c1_prod = ciphertext1.0.compose(&ciphertext2.0);
        let c2_prod = ciphertext1.1.compose(&ciphertext2.1);
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_H_bound);
        let c1 = c1_prod.compose(&self.generator_H.pow(&r)).reduce();
        let c2 = c2_prod.compose(&public_key.pow(&r)).reduce();
        (c1, c2)
    }

    fn scale_ciphertext<Rng: CryptoRng + rand_core::RngCore>(
        &self,
        rng: &mut Rng,
        public_key: BQF,
        ciphertext: (BQF, BQF),
        scaling_factor: Z,
    ) -> (BQF, BQF) {
        let c1_scaled = ciphertext.0.pow(&scaling_factor);
        let c2_scaled = ciphertext.1.pow(&scaling_factor);
        let r = Z::sample_range(rng, &Z::from(0), &self.class_number_H_bound);
        let c1 = c1_scaled.compose(&self.generator_H.pow(&r)).reduce();
        let c2 = c2_scaled.compose(&public_key.pow(&r)).reduce();
        (c1, c2)
    }

    fn generator_f(&self) -> BQF {
        self.generator_F.clone()
    }

    fn generator_h(&self) -> BQF {
        self.generator_H.clone()
    }
}

fn power_f<Z, BQF>(f: &BQF, disc: &Z, q: &Z, cleartext: &Z) -> BQF
where
    BQF: BinaryQuadraticForm<Z> + Clone + std::fmt::Debug,
    Z: std::clone::Clone + std::fmt::Debug + z::Z,
{
    let mut a;
    let mut b;
    let mut c = cleartext.invert_mod(q).unwrap();
    if !c.is_odd() {
        c = c.sub(q);
    }

    /*    [ q^2, Lm*q, ((Lm*q)^2-Delta_q)/(4*q^2) ]
     * =  [ q^2, Lm*q, ((Lm*q)^2-q^2*Delta_K)/(4*q^2) ]
     * =  [ q^2, Lm*q, (Lm^2-Delta_K)/4 ]
     */
    a = q.sqr(); /* a = q^2 */
    b = c.mul(q); /* b = Lm*q */
    c = c.sqr();
    c = c.sub(disc);
    c.divide_by_4(); /* c = (Lm^2-Delta_K)/4 */

    return BQF::new(&a, &b, &c);
}

#[derive(Error, Debug, Clone)]
pub enum KernelError {
    #[error("The form is not in the kernel")]
    NotInKernel,
    #[error("an error occurred: {0}")]
    Other(String),
}

impl<Z, BQF> ClHSMqInstance<Z, BQF>
where
    Z: crate::z::Z + Debug + Clone,
    BQF: BinaryQuadraticForm<Z> + Clone,
{
    pub fn dlog_solve_F(&self, fm: &BQF) -> Z {
        let mut m = Z::from(0u64);
        if !fm.equals(&fm.identity()) {
            //let (u, _j) = Z::remove(&fm.b(), &self.q); /* j, u such that ft.b = q^j*u */
            /* tm = 1/u mod M=q */
            //println!("j={}", _j);
            // assert j = 1? (as k=1)? (and is that true?)
            m = fm.b().divide_exact(&self.q).invert_mod(&self.q).unwrap();
        }
        m
    }
}

fn compute_generator_H<
    Z: std::clone::Clone + std::fmt::Debug + z::Z,
    BQF: BinaryQuadraticForm<Z> + Clone,
>(
    q: &Z,
    discriminant_conductor_q: &Z,
) -> BQF {
    let mut l = Z::from(2);
    while discriminant_conductor_q.kronecker(&l) != 1 {
        l = l.next_prime();
    }
    BQF::prime_form(discriminant_conductor_q, &l)
        .double()
        .reduce()
        .pow(q)
        .reduce()
}

#[test]
fn test_encrypt_decrypt_identity() {
    let q = rug::Integer::from_str_radix(
        "52435875175126190479447740508185965837690552500527637822603658699938581184513",
        10,
    )
    .unwrap();
    let security_level = SecurityLevel::SecLvl112;
    let mut rng = &mut OsRng;

    println!("before new");
    let instance = ClHSMqInstance::new(q.clone(), security_level, &mut rng);

    println!("before keygen");

    let (public_key, private_key) = instance.keygen(&mut rng);

    println!("generated keys = pk={:?} sk={}", public_key, private_key);
    let cleartext = rug::Integer::from(13344545); // Sample cleartext
    let ciphertext: (BQF<Integer>, BQF<Integer>) =
        instance.encrypt(&public_key, &cleartext, &mut rng);

    println!("generated ciphertext {:?}", ciphertext);

    println!("before decrypt");

    let decrypted = instance.decrypt(&private_key, &ciphertext);

    assert_eq!(
        decrypted, cleartext,
        "Decrypted text does not match the original cleartext"
    );
}
