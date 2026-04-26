/*
 * ============================================================
 *  rlwe_kem.h  —  Ring-LWE KEM  v3  (IND-CCA2 via FO Transform)
 *
 *  Scheme : Ring-LWE KEM with Fujisaki-Okamoto (FO) Transform
 *  Ring   : Z_q[x]/(x^n + 1),  n=256, q=3329  (ML-KEM/Kyber compatible)
 *  Noise  : CBD(η=2)
 *  XOF    : SHAKE256  (NIST FIPS 202)  — replaces xoshiro128**
 *  Hash   : SHA3-256  (NIST FIPS 202)
 *  Secret : 256-bit shared key
 *
 *  Security level:
 *    CPA-PKE  → Ring-LWE hardness assumption
 *    IND-CCA2 → Fujisaki-Okamoto transform (re-encapsulation check)
 *
 *  Side-channel protections:
 *    - Constant-time encode/decode  (no secret-dependent branches)
 *    - Constant-time ct_memcmp      (no early exit)
 *    - Constant-time ct_select      (mask-based, no branch)
 *    - Zero dynamic allocation      (stack only)
 *    - SHAKE256 for all PRF/XOF     (replaces non-standard PRNG)
 *
 *  FO Transform (Decaps):
 *    1. Unpack ct → (u, v)
 *    2. Decrypt m' = decode(v - u·s)
 *    3. Re-encapsulate: (ct', ss') = Encaps(pk, m')
 *    4. If ct' == ct  →  output SHA3-256(m' ‖ H(pk))
 *       If ct' != ct  →  output SHA3-256(fail_seed ‖ H(pk))
 *    Steps 3-4 are computed in constant time (no branch on result).
 *
 *  Wire sizes (packed, q=3329 → 12 bits/coeff):
 *    RLWE_PK_BYTES   = 32 + 384  = 416  bytes
 *    RLWE_SK_BYTES   = 384 + 416 + 32 + 32 = 864  bytes
 *                      (s ‖ pk ‖ H(pk) ‖ fail_seed)
 *    RLWE_CT_BYTES   = 384 + 384       = 768  bytes
 *    RLWE_SS_BYTES   = 32              bytes
 *
 *  Compile: gcc -O2 -Wall -Wextra -o rlwe_kem rlwe_kem.c randombytes.c
 * ============================================================
 */

#ifndef RLWE_KEM_H
#define RLWE_KEM_H

#include <stdint.h>
#include <stddef.h>

/* ── Ring parameters (ML-KEM / Kyber-512 compatible) ──────── */
#define RLWE_N        256
#define RLWE_Q        3329      /* ML-KEM standard q  */
#define RLWE_ETA      2
#define RLWE_SS_BYTES 32
#define RLWE_HASH_BYTES 32

/*
 * Packed sizes:
 *   q=3329 < 2^12=4096  →  12 bits per coefficient.
 *   256 coefficients × 12 bits = 3072 bits = 384 bytes.
 */
#define RLWE_POLY_PACKED_BYTES  384   /* 256 * 12 / 8              */
#define RLWE_SEED_BYTES          32   /* seed for matrix A (public) */

#define RLWE_PK_BYTES  (RLWE_SEED_BYTES + RLWE_POLY_PACKED_BYTES)  /* 416 */
#define RLWE_CT_BYTES  (RLWE_POLY_PACKED_BYTES * 2)                 /* 768 */

/*
 * Extended secret key (FO Transform):
 *   sk_inner  : pack(s)            — RLWE_POLY_PACKED_BYTES bytes
 *   pk        : full public key    — RLWE_PK_BYTES bytes
 *   H_pk      : SHA3-256(pk)       — RLWE_HASH_BYTES bytes
 *   fail_seed : random seed for    — RLWE_HASH_BYTES bytes
 *               implicit rejection
 */
#define RLWE_SK_BYTES  (RLWE_POLY_PACKED_BYTES + RLWE_PK_BYTES + \
                        RLWE_HASH_BYTES + RLWE_HASH_BYTES)         /* 864 */

/*
 * ── kem_keypair ──────────────────────────────────────────────
 *
 * Generate a public/secret key pair (IND-CCA2).
 * The secret key includes: s ‖ pk ‖ H(pk) ‖ fail_seed
 *
 *   pk  [out] RLWE_PK_BYTES  serialized public key
 *   sk  [out] RLWE_SK_BYTES  serialized extended secret key
 *
 * Returns 0 on success, -1 on entropy failure.
 */
int kem_keypair(uint8_t pk[RLWE_PK_BYTES],
                uint8_t sk[RLWE_SK_BYTES]);

/*
 * ── kem_encaps ───────────────────────────────────────────────
 *
 * Encapsulate: generate a shared secret and encrypt it.
 * Uses SHAKE256 for all randomness derivation.
 *
 *   ct  [out] RLWE_CT_BYTES  ciphertext
 *   ss  [out] RLWE_SS_BYTES  shared secret = SHA3-256(m ‖ H(pk))
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
 * Decapsulate: recover shared secret (IND-CCA2 via FO Transform).
 *
 * SECURITY: Performs re-encapsulation to detect ciphertext tampering.
 * Output is determined by ct_memcmp in constant time — no timing leak.
 * On invalid ciphertext: outputs a deterministic pseudorandom value
 * (implicit rejection), preventing adaptive attacks.
 *
 * Never fails: always writes RLWE_SS_BYTES to ss.
 *
 *   ss  [out] RLWE_SS_BYTES  shared secret (or rejection value)
 *   ct  [in]  RLWE_CT_BYTES  ciphertext
 *   sk  [in]  RLWE_SK_BYTES  extended secret key
 */
void kem_decaps(uint8_t ss[RLWE_SS_BYTES],
                const uint8_t ct[RLWE_CT_BYTES],
                const uint8_t sk[RLWE_SK_BYTES]);

#endif /* RLWE_KEM_H */
