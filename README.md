# libctkem

**libctkem** is a lightweight and efficient C library implementing a Key Encapsulation Mechanism (KEM) based on the **Ring Learning With Errors (RLWE)** problem. The project is specifically optimized for ARM architectures and serves as a research platform for post-quantum cryptographic primitives.

## Core Objective
The library is designed to facilitate secure key exchange resistant to potential quantum computing threats. It implements the mathematical core of lattice-based algorithms (similar to ML-KEM/Kyber) with a strict focus on minimal code footprint and high execution speed on mobile and embedded systems.

## Key Features

### 1. Side-Channel Resistance (Constant-Time)
All security-critical operations are implemented in **constant-time**. This ensures that the execution flow is independent of secret key values, protecting the system against:
* Timing attacks.
* Microarchitectural data leakage.
* Branch-prediction exploitation (via branchless masks).

### 2. Optimized Polynomial Arithmetic
The library features highly tuned arithmetic for polynomial rings:
* **Barrett Reduction:** An optimized modular reduction algorithm that avoids expensive division instructions, tailored for the specific modulus used in RLWE.
* **Efficient Multiplication:** A streamlined implementation of polynomial multiplication for $N=256$, ensuring predictable instruction pipelining and cache efficiency.

### 3. ARM-Centric Design
The codebase is developed with **ARMv8 (Cortex-A/M)** architectures in mind. Algorithmic choices prioritize register-heavy operations and minimize memory overhead, making it ideal for resource-constrained environments.

### 4. Zero-Dependency Integration
`libctkem` consists of a minimal set of `.c` and `.h` files. It has no external dependencies, making it easy to port to "bare metal" firmware, RTOS, or mobile platforms.

## Technical Specifications
* **Primitive:** RLWE (Ring Learning With Errors).
* **PRNG:** Utilizes the fast `xoshiro128**` generator (modular design allows for easy replacement with SHAKE256/FIPS 202 for standardized compliance).
* **Scheme:** Full KEM cycle (Key Generation, Encapsulation, Decapsulation).

## Project Status
This is an **experimental library** designed to explore the performance boundaries of RLWE-based schemes on mobile hardware. It serves as a proof-of-concept for high-performance, side-channel-resistant cryptography on edge devices.

## Disclaimer
This project is intended for research and educational purposes. Production use should be preceded by a formal cryptographic audit, particularly regarding the entropy sources of the random number generator.
