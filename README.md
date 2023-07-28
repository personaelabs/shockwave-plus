# Shockwave+

## Motivation

Field-agnostic proof systems provide efficient proving by allowing *native-field arithmetic.* This property is significant for client-side programmable signatures, zkEVMs, and various other scenarios where proof of statements about primitives traditionally considered "zk-unfriendly" is required.

We aim to build a zero-knowledge proof system that is field-agnostic, efficient, and effortlessly composable with a SNARK with $O (1)$ proof size (e.g. Groth16, PLONK) to achieve minimal verification complexity.

## Overview

Shockwave is a variant of [Brakedown](https://eprint.iacr.org/2021/1043) that uses Reed-Solomon code instead of a linear-time encodable code. Brakedown has a linear-time prover and is field-agnostic (i.e. works over all finite fields), but its proofs are concretely larger than Shockwave’s. Shockwave provides shorter proofs and lower verification time but requires an FFT-friendly field to achieve $O (n\log{n})$ proving time. 

**Shockwave+** is an extension of Shockwave that works over all finite fields by using [ECFFT](https://arxiv.org/pdf/2107.08473.pdf) instead of FFT for low-degree extension of polynomial evaluations. It inherits the smaller proofs of Shockwave and is also field-agnostic. It uses the EXTEND operation from [ECFFT](https://arxiv.org/pdf/2107.08473.pdf) to run Reed-Solomon encoding in $O (n\log{n})$ time.

### Crates

[shockwave_plus](/shockwave_plus/) contains the prover/verifier for a zero-knowledge proof of R1CS satisfiability. It’s based on the PIOP from [Spartan](https://eprint.iacr.org/2019/550.pdf), and uses the multilinear polynomial commitment scheme implemented in [tensor_pcs](/tensor_pcs/).

The EXTEND operation is implemented in a separate crate [ecfft](https://github.com/DanTehrani/ecfft) and is used in [tensor_pcs](/tensor_pcs/).

### Zero-Knowledge

We use the zero-knowledge sum-check protocol from [Libra](https://eprint.iacr.org/2019/317.pdf) to transform the Spartan PIOP into a zero-knowledge PIOP. And use a technique from [BCG+17](https://eprint.iacr.org/2017/872.pdf) to make the polynomial commitment scheme zero-knowledge.


## Benchmarks

TBD

## Future work
- [ ]  Employ *self-recursion* techniques from [Vortex](https://eprint.iacr.org/2022/1633.pdf)/[Orion](https://eprint.iacr.org/2022/1010.pdf) to make the proofs smaller.
- [ ]  Support richer frontends (CCS, PLONKish).

## Run tests
```bash
cargo test
```
