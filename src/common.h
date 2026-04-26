/*
 * ============================================================
 * common.h  —  Global Parameters for Ring-LWE KEM
 *
 * Scheme: Kyber-512 compatible (N=256, Q=3329)
 * Target: ARM64 (NEON) and Portable C
 * ============================================================
 */

#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>

/* --- Algorithm Parameters --- */
#define RLWE_N            256      /* Degree of the polynomial ring */
#define RLWE_Q            3329     /* Prime modulus */
#define RLWE_K            2        /* Rank (Kyber-512 setting) */

/* --- Byte Lengths for API --- */
#define RLWE_SYMBYTES     32       /* 256-bit seed/entropy */
#define RLWE_SS_BYTES     32       /* Shared Secret length */
#define RLWE_PK_BYTES     800      /* Public Key: (k*12*N/8 + 32) */
#define RLWE_SK_BYTES     1632     /* Secret Key: (pk + sk + h(pk) + z) */
#define RLWE_CT_BYTES     768      /* Ciphertext: (k*10*N/8 + 4*N/8) */

/* --- Montgomery Constants (Q = 3329) --- */
#define MONT              2285     /* 2^16 mod Q */
#define QINV              3327     /* -Q^-1 mod 2^16 */

/* --- Compiler Optimization Macros --- */

/* Pointer aliasing hint for the compiler (Clang/GCC) */
#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

/* Inline hint for hot-path functions */
#define INLINE static inline

/* --- Data Structures --- */

/* * poly: Polynomial structure.
 * Aligned to 16 bytes for optimal NEON vector load/store (LDP/STP).
 */
typedef struct {
    int16_t coeffs[RLWE_N];
} __attribute__((aligned(16))) poly;

/* --- Extern Function Prototypes --- */

/* System entropy source (implemented in randombytes.c) */
void randombytes(uint8_t *out, size_t outlen);

#endif /* COMMON_H */
