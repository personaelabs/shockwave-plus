# Shockwave+

## Overview

Shockwave is a variant of [Brakedown](https://eprint.iacr.org/2021/1043) that uses Reed-Solomon code instead of a linear-time encodable code. **Shockwave+** is an extension of Shockwave that works over all finite fields by using [ECFFT](https://arxiv.org/pdf/2107.08473.pdf) instead of FFT for low-degree extension of polynomial evaluations.

Brakedown has a linear-time prover and is *field-agnostic* (i.e. works over all finite fields), but its proofs are concretely larger than Shockwave’s.

Shockwave provides shorter proofs and lower verification time but requires an FFT-friendly field to achieve $O (n\log{n})$ proving time. 

Shockwave+ inherits the smaller proofs of Shockwave and is also *field-agnostic*. It uses the EXTEND operation from [ECFFT](https://arxiv.org/pdf/2107.08473.pdf) to run Reed-Solomon encoding in $n\log{n}$ time.

**Crates**
[shockwave_plus](/shockwave_plus/) contains the prover/verifier for a zero-knowledge proof of R1CS satisfiability. It’s based on the PIOP from [Spartan](https://eprint.iacr.org/2019/550.pdf), and uses the multilinear polynomial commitment scheme implemented in [tensor_pcs](/tensor_pcs/).

**Zero-Knowledge**

We use the zero-knowledge sum-check protocol from Libra to transform the Spartan PIOP into a zero-knowledge PIOP. And use a technique from [BCG+17](https://eprint.iacr.org/2017/872.pdf) to make the polynomial commitment scheme zero-knowledge.



The EXTEND operation is implemented in a separate crate [ecfft](https://github.com/DanTehrani/ecfft) and is used in [tensor_pcs](/tensor_pcs/).

## Benchmarks

TBD

## Future work

- [ ]  Support richer frontends (CCS, PLONKish).
- [ ]  Employ *self-recursion* techniques from [Vortex](https://eprint.iacr.org/2022/1633.pdf)/[Orion](https://eprint.iacr.org/2022/1010.pdf) to make the proofs smaller.

## Run tests
```bash
cargo test
```