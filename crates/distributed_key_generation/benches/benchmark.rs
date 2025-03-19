use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use rand::rngs::OsRng;
use rug::Integer;
use blstrs::{G1Projective, Scalar};
use distributed_key_generation::{
    setup,
    generate_key_pair,
    verify_public_key,
    generate_dealing,
    verify_dealing,
    SecurityLevel,
};

fn bench_setup(c: &mut Criterion) {
    let mut group = c.benchmark_group("Setup");
    
    for &participants in [5, 10, 20, 50].iter() {
        let threshold = participants / 2;
        
        group.bench_with_input(
            BenchmarkId::new("setup", participants),
            &(participants, threshold),
            |b, &(n, t)| {
                b.iter(|| {
                    setup::<_, Scalar, rug::Integer>(
                        n,
                        t,
                        SecurityLevel::SecLvl112,
                        &mut OsRng,
                    )
                })
            },
        );
    }
    
    group.finish();
}

fn bench_key_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("KeyGeneration");
    
    for &participants in [5, 10, 20, 50].iter() {
        let threshold = participants / 2;
        let pp = setup::<_, Scalar, rug::Integer>(
            participants,
            threshold,
            SecurityLevel::SecLvl112,
            &mut OsRng,
        );
        
        group.bench_with_input(
            BenchmarkId::new("generate_key_pair", participants),
            &pp,
            |b, pp| {
                b.iter(|| generate_key_pair::<OsRng, Integer>(pp, &mut OsRng))
            },
        );
    }
    
    group.finish();
}

fn bench_verify_public_key(c: &mut Criterion) {
    let mut group = c.benchmark_group("VerifyPublicKey");
    
    for &participants in [5, 10, 20, 50].iter() {
        let threshold = participants / 2;
        let pp = setup::<_, Scalar, rug::Integer>(
            participants,
            threshold,
            SecurityLevel::SecLvl112,
            &mut OsRng,
        );
        let (pk, _, pop) = generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);
        
        group.bench_with_input(
            BenchmarkId::new("verify_public_key", participants),
            &(pp, pk, pop),
            |b, (pp, pk, pop)| {
                b.iter(|| verify_public_key(pp, pk, pop))
            },
        );
    }
    
    group.finish();
}

fn bench_generate_dealing(c: &mut Criterion) {
    let mut group = c.benchmark_group("GenerateDealing");
    
    for &participants in [5, 10, 20, 50].iter() {
        let threshold = participants / 2;
        let pp = setup::<_, Scalar, rug::Integer>(
            participants,
            threshold,
            SecurityLevel::SecLvl112,
            &mut OsRng,
        );
        
        // Generate keys for all participants
        let mut pks = vec![];
        for _ in 0..participants {
            let (pk, _, _) = generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);
            pks.push(pk);
        }
        
        group.bench_with_input(
            BenchmarkId::new("generate_dealing", participants),
            &(pp, pks.clone()),
            |b, (pp, pks)| {
                b.iter(|| {
                    generate_dealing::<
                        _,
                        _,
                        G1Projective,
                        Scalar,
                        Vec<Scalar>,
                    >(pp, pks, &mut OsRng)
                })
            },
        );
    }
    
    group.finish();
}

fn bench_verify_dealing(c: &mut Criterion) {
    let mut group = c.benchmark_group("VerifyDealing");
    
    for &participants in [5, 10, 20, 50].iter() {
        let threshold = participants / 2;
        let pp = setup::<_, Scalar, rug::Integer>(
            participants,
            threshold,
            SecurityLevel::SecLvl112,
            &mut OsRng,
        );
        
        // Generate keys for all participants
        let mut pks = vec![];
        for _ in 0..participants {
            let (pk, _, _) = generate_key_pair::<OsRng, Integer>(&pp, &mut OsRng);
            pks.push(pk);
        }
        
        // Generate one dealing
        let dealing = generate_dealing::<
            _,
            _,
            G1Projective,
            Scalar,
            Vec<Scalar>,
        >(&pp, &pks, &mut OsRng);
        
        group.bench_with_input(
            BenchmarkId::new("verify_dealing", participants),
            &(pp, pks.clone(), dealing.clone()),
            |b, (pp, pks, dealing)| {
                b.iter(|| {
                    verify_dealing::<
                        Integer,
                        G1Projective,
                        Scalar,
                        Vec<Scalar>,
                    >(pp, pks, dealing)
                })
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_setup,
    bench_key_generation,
    bench_verify_public_key,
    bench_generate_dealing,
    bench_verify_dealing,
);
criterion_main!(benches);
