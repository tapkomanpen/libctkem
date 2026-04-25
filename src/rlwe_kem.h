/*
 * ============================================================
 *  rlwe_kem.h  —  Public API for Compact Ring-LWE KEM
 *
 *  Scheme : Ring-LWE based Key Encapsulation (NewHope-style)
 *  Ring   : Z_q[x]/(x^n + 1),  n=256, q=7681
 *  Noise  : CBD(η=2), σ≈1
 *  Secret : 256-bit shared key
 *
 *  Security properties:
 *    - Constant-time encode/decode (no secret-dependent branches)
 *    - Constant-time comparison (ct_memcmp)
 *    - Zero dynamic allocation (stack-only)
 *    - 13-bit coefficient packing (20% smaller keys/ciphertexts)
 *
 *  Wire sizes (packed):
 *    RLWE_PK_BYTES   = 32 + 416  = 448  bytes
 *    RLWE_SK_BYTES   =      416  = 416  bytes
 *    RLWE_CT_BYTES   = 416 + 416 = 832  bytes
 *    RLWE_SS_BYTES   =       32  =  32  bytes
 *
 *  Compile: gcc -O2 -Wall -o rlwe_kem rlwe_kem.c
 * ============================================================
 */

#ifndef RLWE_KEM_H
#define RLWE_KEM_H

#include <stdint.h>
#include <stddef.h>

/* ── Ring parameters ──────────────────────────────────────── */
#define RLWE_N        256
#define RLWE_Q        7681
#define RLWE_ETA      2
#define RLWE_SS_BYTES 32

/*
 * Packed sizes:
 *   Each coefficient ∈ [0, 7681) needs ceil(log2(7681)) = 13 bits.
 *   256 coefficients × 13 bits = 3328 bits = 416 bytes exactly.
 */
#define RLWE_POLY_PACKED_BYTES  416   /* 256 * 13 / 8              */
#define RLWE_SEED_BYTES          32   /* seed for matrix A         */
#define RLWE_PK_BYTES  (RLWE_SEED_BYTES + RLWE_POLY_PACKED_BYTES)
#define RLWE_SK_BYTES  (RLWE_POLY_PACKED_BYTES)
#define RLWE_CT_BYTES  (RLWE_POLY_PACKED_BYTES + RLWE_POLY_PACKED_BYTES)

/*
 * ── kem_keypair ──────────────────────────────────────────────
 *
 * Generate a public/secret key pair.
 *
 *   pk  [out] RLWE_PK_BYTES  serialized public key
 *   sk  [out] RLWE_SK_BYTES  serialized secret key
 *
 * Returns 0 on success, -1 on entropy failure.
 */
int kem_keypair(uint8_t pk[RLWE_PK_BYTES],
                uint8_t sk[RLWE_SK_BYTES]);

/*
 * ── kem_encaps ───────────────────────────────────────────────
 *
 * Encapsulate: generate a shared secret and encrypt it.
 *
 *   ct  [out] RLWE_CT_BYTES  ciphertext
 *   ss  [out] RLWE_SS_BYTES  shared secret
 *   pk  [in]  RLWE_PK_BYTES  recipient's public key
 *
 * Returns 0 on success, -1 on entropy failure.
 */
int kem_encaps(uint8_t ct[RLWE_CT_BYTES],
               uint8_t ss[RLWE_SS_BYTES],
               const uint8_t pk[RLWE_PK_BYTES]);

/*
 * ── kem_decaps ───────────────────────────────────────────────
 *
 * Decapsulate: recover shared secret from ciphertext.
 * Runs in constant time (no secret-dependent branches).
 *
 *   ss  [out] RLWE_SS_BYTES  shared secret
 *   ct  [in]  RLWE_CT_BYTES  ciphertext
 *   sk  [in]  RLWE_SK_BYTES  secret key
 */
void kem_decaps(uint8_t ss[RLWE_SS_BYTES],
                const uint8_t ct[RLWE_CT_BYTES],
                const uint8_t sk[RLWE_SK_BYTES]);

#endif /* RLWE_KEM_H */
