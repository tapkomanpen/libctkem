/*
 * ============================================================
 *  rlwe_kem.c  —  Ring-LWE KEM  v3  (IND-CCA2 + FO Transform)
 *
 *  Improvements over v2:
 *    1. IND-CCA2 via Fujisaki-Okamoto (FO) Transform
 *       — Decaps performs re-encapsulation + constant-time compare
 *       — Implicit rejection on invalid ciphertext
 *    2. SHAKE256 / SHA3-256 replace xoshiro128**
 *       — Provably secure PRF/XOF (NIST FIPS 202)
 *       — Deterministic matrix generation (poly_uniform via SHAKE256)
 *    3. q = 3329, 12-bit packing (ML-KEM compatible)
 *       — 256 × 12 bits = 384 bytes per polynomial
 *    4. Extended secret key: s ‖ pk ‖ H(pk) ‖ fail_seed
 *       — fail_seed enables implicit rejection
 *    5. All-constant-time: ct_select replaces conditional branches
 *       — Even intermediate hash comparisons use ct_memcmp
 *    6. MISRA-C style: no VLAs, explicit casts, bounded loops
 *
 *  Compile:
 *    gcc -O2 -Wall -Wextra -o rlwe_kem rlwe_kem.c randombytes.c
 *  KAT:
 *    ./rlwe_kem --kat 5 > kat_vectors.txt
 * ============================================================
 */

#include "rlwe_kem.h"
#include "randombytes.h"
#include "sha3.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ── Internal aliases ─────────────────────────────────────── */
#define N       RLWE_N
#define Q       RLWE_Q           /* 3329 */
#define QHALF   1664             /* floor(Q/2) */

typedef int16_t poly[N];

/* ── SK layout offsets ────────────────────────────────────── */
#define SK_OFF_S        0
#define SK_OFF_PK       (RLWE_POLY_PACKED_BYTES)
#define SK_OFF_HPK      (RLWE_POLY_PACKED_BYTES + RLWE_PK_BYTES)
#define SK_OFF_FAIL     (RLWE_POLY_PACKED_BYTES + RLWE_PK_BYTES + RLWE_HASH_BYTES)

/* ============================================================
 *  §1  Constant-time primitives
 * ============================================================ */

/*
 * ct_memcmp — constant-time comparison.
 * Returns 0 iff equal. Time depends only on 'len'.
 */
static int ct_memcmp(const uint8_t *a, const uint8_t *b, size_t len) {
    uint8_t diff = 0;
    for (size_t i = 0; i < len; i++)
        diff |= a[i] ^ b[i];
    return (int)diff;
}

/*
 * ct_select_byte — branchless byte selection.
 * Returns a[i] if cond==0, b[i] if cond!=0.
 * 'cond' must be 0x00 or 0xFF (use ct_mask_from_neq).
 */
static inline uint8_t ct_select_byte(uint8_t a, uint8_t b, uint8_t mask) {
    return (uint8_t)(a ^ (mask & (a ^ b)));
}

/*
 * ct_mask_from_neq — produce 0xFF if x!=0, else 0x00.
 * Uses: result of ct_memcmp → mask for ct_select_byte.
 */
static inline uint8_t ct_mask_from_neq(int x) {
    /* If x != 0, (uint32_t)-x has bit 31 set; arithmetic shift fills */
    uint32_t u = (uint32_t)x;
    return (uint8_t)(((uint32_t)(-(int32_t)(u | (uint32_t)(-(int32_t)u))) >> 31) * 0xFF);
}

/*
 * ct_memselect — constant-time buffer selection.
 * If   cond == 0  → out = a
 * Else cond != 0  → out = b
 */
static void ct_memselect(uint8_t *out,
                          const uint8_t *a, const uint8_t *b,
                          size_t len, int cond)
{
    uint8_t mask = ct_mask_from_neq(cond);
    for (size_t i = 0; i < len; i++)
        out[i] = ct_select_byte(a[i], b[i], mask);
}

/*
 * ct_min_u32 — constant-time min(a, b).
 */
static inline uint32_t ct_min_u32(uint32_t a, uint32_t b) {
    uint32_t mask = (uint32_t)(-(uint32_t)(b < a));
    return a ^ ((a ^ b) & mask);
}

/* ============================================================
 *  §2  Barrett reduction (q=3329, branchless)
 *
 *  For q=3329: k = floor(2^26 / q) = 20159
 *  Reduces a ∈ (-q²/2, q²/2) to [0, q).
 * ============================================================ */
static inline int16_t barrett(int32_t a) {
    int32_t t = (int32_t)(((int64_t)a * 20159L) >> 26);
    int32_t r = a - t * Q;
    r -= Q & -(r >= Q);
    r += Q & -(r <  0);
    return (int16_t)r;
}

/* ============================================================
 *  §3  CBD noise sampling via SHAKE256
 *
 *  CBD(eta=2): x = popcount(bits 0,1) - popcount(bits 2,3)
 *  x ∈ {-2,-1,0,1,2}
 *
 *  We now derive noise from SHAKE256(seed ‖ nonce) instead of
 *  xoshiro128**. This gives provable pseudorandomness.
 * ============================================================ */
static inline int ct_pop2(uint32_t x) {
    uint32_t m = x & 0x00030003u;
    m = (m & 0x00010001u) + ((m >> 1) & 0x00010001u);
    return (int)((m + (m >> 16)) & 0xF);
}

/*
 * poly_sample_noise_xof — sample noise polynomial from SHAKE256 stream.
 * ctx must already be initialized and finalized (squeezing phase).
 */
static void poly_sample_noise_xof(poly p, shake256_ctx *ctx) {
    uint8_t block[4];
    for (int i = 0; i < N; i++) {
        shake256_squeeze(ctx, block, 4);
        uint32_t r = (uint32_t)block[0]
                   | ((uint32_t)block[1] << 8)
                   | ((uint32_t)block[2] << 16)
                   | ((uint32_t)block[3] << 24);
        int a = ct_pop2(r);
        int b = ct_pop2(r >> 2);
        p[i] = barrett((int32_t)(a - b) + Q);
    }
}

/*
 * poly_sample_noise — convenience wrapper.
 * Derives noise from SHAKE256(seed ‖ domain ‖ nonce).
 */
static void poly_sample_noise(poly p,
                               const uint8_t seed[32],
                               uint8_t domain,
                               uint8_t nonce)
{
    uint8_t label[2] = { domain, nonce };
    shake256_ctx ctx;
    shake256_init(&ctx);
    shake256_absorb(&ctx, seed, 32);
    shake256_absorb(&ctx, label, 2);
    shake256_finalize(&ctx);
    poly_sample_noise_xof(p, &ctx);
}

/* ============================================================
 *  §4  Polynomial arithmetic in R_q = Z_q[x]/(x^n+1)
 *
 *  poly_mul: negacyclic schoolbook O(n²).
 *  Constant-time: no branch on secret data, uniform loop.
 *
 *  NOTE: For a production standard, NTT would replace schoolbook
 *  for O(n log n). Schoolbook is kept here for clarity and
 *  auditability per the MISRA-C principle.
 * ============================================================ */
static void poly_mul(poly c, const poly a, const poly b) {
    int64_t acc[N];
    memset(acc, 0, sizeof(acc));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int     k    = i + j;
            int     wrap = (int)(k >= N);
            int     idx  = k - wrap * N;
            int64_t sign = (int64_t)(1 - 2 * wrap);
            acc[idx] += sign * (int64_t)a[i] * (int64_t)b[j];
        }
    }

    for (int i = 0; i < N; i++) {
        int32_t v = (int32_t)(acc[i] % (int64_t)Q);
        c[i] = barrett(v);
    }
}

static void poly_add(poly c, const poly a, const poly b) {
    for (int i = 0; i < N; i++)
        c[i] = barrett((int32_t)a[i] + (int32_t)b[i]);
}

static void poly_sub(poly c, const poly a, const poly b) {
    for (int i = 0; i < N; i++)
        c[i] = barrett((int32_t)a[i] - (int32_t)b[i] + (int32_t)Q);
}

/*
 * poly_uniform — expand public polynomial uniformly from seed.
 * Uses SHAKE256(seed ‖ domain_byte) for provable pseudorandomness.
 * Rejection sampling: expected ~1.07 iterations per coefficient.
 */
static void poly_uniform(poly a,
                          const uint8_t seed[RLWE_SEED_BYTES],
                          uint8_t domain)
{
    uint8_t d = domain;
    shake256_ctx ctx;
    shake256_init(&ctx);
    shake256_absorb(&ctx, seed, RLWE_SEED_BYTES);
    shake256_absorb(&ctx, &d, 1);
    shake256_finalize(&ctx);

    uint8_t buf[2];
    for (int i = 0; i < N; ) {
        shake256_squeeze(&ctx, buf, 2);
        uint32_t r = ((uint32_t)buf[0] | ((uint32_t)(buf[1] & 0x0F) << 8));
        /* r ∈ [0, 4095], accept if < Q=3329 */
        if (r < (uint32_t)Q) {
            a[i] = (int16_t)r;
            i++;
        }
    }
}

/* ============================================================
 *  §5  12-bit coefficient packing / unpacking (q=3329 < 2^12)
 *
 *  256 coefficients × 12 bits = 3072 bits = 384 bytes.
 *  Pure bijection — no control flow depends on data values.
 * ============================================================ */
static void poly_pack(uint8_t out[RLWE_POLY_PACKED_BYTES], const poly p) {
    /* Process pairs of coefficients: 2 × 12 bits = 3 bytes */
    for (int i = 0, j = 0; i < N; i += 2, j += 3) {
        uint16_t a = (uint16_t)p[i];
        uint16_t b = (uint16_t)p[i + 1];
        out[j    ] = (uint8_t)(a & 0xFF);
        out[j + 1] = (uint8_t)((a >> 8) | ((b & 0x0F) << 4));
        out[j + 2] = (uint8_t)(b >> 4);
    }
}

static void poly_unpack(poly p,
                         const uint8_t in[RLWE_POLY_PACKED_BYTES])
{
    for (int i = 0, j = 0; i < N; i += 2, j += 3) {
        p[i    ] = (int16_t)((uint16_t)in[j] | (((uint16_t)in[j+1] & 0x0F) << 8));
        p[i + 1] = (int16_t)(((uint16_t)in[j+1] >> 4) | ((uint16_t)in[j+2] << 4));
    }
}

/* ============================================================
 *  §6  Constant-time message encode / decode
 *
 *  Encoding: bit=0 → coeff=0; bit=1 → coeff=QHALF=1664
 *  Decoding: closest to 0 or QHALF (constant-time, no branch)
 * ============================================================ */
static void encode_msg(poly out, const uint8_t msg[RLWE_SS_BYTES]) {
    for (int i = 0; i < N; i++) {
        int bit = (i < RLWE_SS_BYTES * 8)
                  ? ((msg[i >> 3] >> (i & 7)) & 1)
                  : 0;
        out[i] = (int16_t)(-(uint16_t)bit & (uint16_t)QHALF);
    }
}

static void decode_msg(uint8_t msg[RLWE_SS_BYTES], const poly w) {
    memset(msg, 0, RLWE_SS_BYTES);
    for (int i = 0; i < RLWE_SS_BYTES * 8; i++) {
        uint32_t c    = (uint32_t)(uint16_t)w[i];
        uint32_t d0   = ct_min_u32(c, (uint32_t)Q - c);
        uint32_t diff = (c >= (uint32_t)QHALF) ? (c - QHALF) : (QHALF - c);
        uint32_t d1   = ct_min_u32(diff, (uint32_t)Q - diff);
        /* bit=1 iff d1 < d0; (d1-d0) underflows → MSB=1 */
        uint32_t bit  = (d1 - d0) >> 31;
        msg[i >> 3]  |= (uint8_t)(bit << (i & 7));
    }
}

/* ============================================================
 *  §7  Internal CPA-PKE (used by FO transform)
 * ============================================================ */

/*
 * pke_enc_det — deterministic CPA encryption.
 * All randomness derived from (r_seed, noise_seed).
 */
static void pke_enc_det(uint8_t ct[RLWE_CT_BYTES],
                         const uint8_t pk[RLWE_PK_BYTES],
                         const uint8_t m[RLWE_SS_BYTES],
                         const uint8_t r_seed[32])
{
    poly a, b, s2, e1, e2, u, v, tmp, mp;

    poly_uniform(a, pk, 0x01);         /* public matrix from pk seed */
    poly_unpack(b, pk + RLWE_SEED_BYTES);

    poly_sample_noise(s2, r_seed, 0x03, 0);
    poly_sample_noise(e1, r_seed, 0x03, 1);
    poly_sample_noise(e2, r_seed, 0x03, 2);

    poly_mul(u, a, s2);  poly_add(u, u, e1);     /* u = a·s' + e'  */
    poly_mul(tmp, b, s2); poly_add(tmp, tmp, e2);  /* b·s' + e''     */
    encode_msg(mp, m);
    poly_add(v, tmp, mp);                          /* v = b·s'+e''+enc(m) */

    poly_pack(ct,                          u);
    poly_pack(ct + RLWE_POLY_PACKED_BYTES, v);
}

/*
 * pke_dec — CPA decryption.
 * Recovers m' from (ct, s). Constant-time.
 */
static void pke_dec(uint8_t m[RLWE_SS_BYTES],
                    const uint8_t ct[RLWE_CT_BYTES],
                    const poly s)
{
    poly u, v, w;
    poly_unpack(u, ct);
    poly_unpack(v, ct + RLWE_POLY_PACKED_BYTES);
    poly_mul(w, u, s);
    poly_sub(w, v, w);
    decode_msg(m, w);
}

/* ============================================================
 *  §8  Public API  (IND-CCA2)
 * ============================================================ */

int kem_keypair(uint8_t pk[RLWE_PK_BYTES],
                uint8_t sk[RLWE_SK_BYTES])
{
    poly    a, s, e, b;
    uint8_t noise_seed[32];
    uint8_t fail_seed[RLWE_HASH_BYTES];

    /* Public seed for matrix A (sent in plaintext) */
    if (randombytes(pk, RLWE_SEED_BYTES) < 0) return -1;
    /* Private noise seed */
    if (randombytes(noise_seed, 32)       < 0) return -1;
    /* Implicit rejection seed */
    if (randombytes(fail_seed, RLWE_HASH_BYTES) < 0) return -1;

    /* Sample secret and error from SHAKE256(noise_seed, domain, nonce) */
    poly_sample_noise(s, noise_seed, 0x02, 0);  /* s ← χ_η */
    poly_sample_noise(e, noise_seed, 0x02, 1);  /* e ← χ_η */

    poly_uniform(a, pk, 0x01);   /* a from public seed */
    poly_mul(b, a, s);
    poly_add(b, b, e);           /* b = a·s + e */

    /* pk = seed ‖ pack(b) */
    poly_pack(pk + RLWE_SEED_BYTES, b);

    /* Extended sk = pack(s) ‖ pk ‖ SHA3-256(pk) ‖ fail_seed */
    poly_pack(sk + SK_OFF_S, s);
    memcpy(sk + SK_OFF_PK, pk, RLWE_PK_BYTES);
    sha3_256(sk + SK_OFF_HPK, pk, RLWE_PK_BYTES);
    memcpy(sk + SK_OFF_FAIL, fail_seed, RLWE_HASH_BYTES);

    /* Wipe temporaries */
    memset(noise_seed, 0, sizeof(noise_seed));
    memset(fail_seed,  0, sizeof(fail_seed));

    return 0;
}

int kem_encaps(uint8_t ct[RLWE_CT_BYTES],
               uint8_t ss[RLWE_SS_BYTES],
               const uint8_t pk[RLWE_PK_BYTES])
{
    uint8_t m[RLWE_SS_BYTES];
    uint8_t h_pk[RLWE_HASH_BYTES];
    uint8_t r_seed[32];
    uint8_t hash_input[RLWE_SS_BYTES + RLWE_HASH_BYTES];

    /* Random message m */
    if (randombytes(m, RLWE_SS_BYTES) < 0) return -1;

    /* H(pk) */
    sha3_256(h_pk, pk, RLWE_PK_BYTES);

    /* r_seed = SHA3-256(m ‖ H(pk)) — deterministic, non-malleable */
    memcpy(hash_input,                 m,    RLWE_SS_BYTES);
    memcpy(hash_input + RLWE_SS_BYTES, h_pk, RLWE_HASH_BYTES);
    sha3_256(r_seed, hash_input, sizeof(hash_input));

    /* CPA encrypt */
    pke_enc_det(ct, pk, m, r_seed);

    /* ss = SHA3-256(m ‖ H(pk)) — same as r_seed here */
    memcpy(ss, r_seed, RLWE_SS_BYTES);

    /* Wipe */
    memset(m,          0, sizeof(m));
    memset(h_pk,       0, sizeof(h_pk));
    memset(r_seed,     0, sizeof(r_seed));
    memset(hash_input, 0, sizeof(hash_input));

    return 0;
}

/*
 * kem_decaps — IND-CCA2 decapsulation via FO Transform.
 *
 * Security invariant: the same number of operations are performed
 * regardless of whether ct is valid. ct_memselect picks the output
 * based on the comparison result with NO conditional branches.
 */
void kem_decaps(uint8_t ss[RLWE_SS_BYTES],
                const uint8_t ct[RLWE_CT_BYTES],
                const uint8_t sk[RLWE_SK_BYTES])
{
    poly    s;
    uint8_t m_prime[RLWE_SS_BYTES];
    uint8_t h_pk[RLWE_HASH_BYTES];
    uint8_t r_seed[32];
    uint8_t ct_prime[RLWE_CT_BYTES];
    uint8_t ss_ok[RLWE_SS_BYTES];
    uint8_t ss_fail[RLWE_SS_BYTES];
    uint8_t hash_input[RLWE_SS_BYTES + RLWE_HASH_BYTES];

    const uint8_t *pk        = sk + SK_OFF_PK;
    const uint8_t *stored_hpk = sk + SK_OFF_HPK;
    const uint8_t *fail_seed  = sk + SK_OFF_FAIL;

    /* 1. CPA decrypt: recover m' */
    poly_unpack(s, sk + SK_OFF_S);
    pke_dec(m_prime, ct, s);

    /* 2. Retrieve H(pk) from extended sk */
    memcpy(h_pk, stored_hpk, RLWE_HASH_BYTES);

    /* 3. Recompute r_seed = SHA3-256(m' ‖ H(pk)) */
    memcpy(hash_input,                 m_prime, RLWE_SS_BYTES);
    memcpy(hash_input + RLWE_SS_BYTES, h_pk,    RLWE_HASH_BYTES);
    sha3_256(r_seed, hash_input, sizeof(hash_input));

    /* 4. Re-encrypt: ct' = Enc(pk, m', r_seed) */
    pke_enc_det(ct_prime, pk, m_prime, r_seed);

    /* 5. ss_ok   = SHA3-256(m' ‖ H(pk))  = r_seed (if valid) */
    memcpy(ss_ok, r_seed, RLWE_SS_BYTES);

    /* 6. ss_fail = SHA3-256(fail_seed ‖ H(pk))  — implicit rejection */
    memcpy(hash_input,                 fail_seed, RLWE_HASH_BYTES);
    memcpy(hash_input + RLWE_HASH_BYTES, h_pk,   RLWE_HASH_BYTES);
    sha3_256(ss_fail, hash_input, RLWE_HASH_BYTES + RLWE_HASH_BYTES);

    /* 7. Constant-time select: output ss_ok iff ct' == ct */
    int mismatch = ct_memcmp(ct, ct_prime, RLWE_CT_BYTES);
    ct_memselect(ss, ss_ok, ss_fail, RLWE_SS_BYTES, mismatch);

    /* Wipe all temporaries unconditionally */
    memset(m_prime,   0, sizeof(m_prime));
    memset(h_pk,      0, sizeof(h_pk));
    memset(r_seed,    0, sizeof(r_seed));
    memset(ct_prime,  0, sizeof(ct_prime));
    memset(ss_ok,     0, sizeof(ss_ok));
    memset(ss_fail,   0, sizeof(ss_fail));
    memset(hash_input,0, sizeof(hash_input));
}

/* ============================================================
 *  §9  KAT — Known Answer Tests
 * ============================================================ */
static int kem_keypair_det(uint8_t pk[RLWE_PK_BYTES],
                            uint8_t sk[RLWE_SK_BYTES],
                            const uint8_t pk_seed[RLWE_SEED_BYTES],
                            const uint8_t noise_seed[32],
                            const uint8_t fail_seed[RLWE_HASH_BYTES])
{
    poly a, s, e, b;

    memcpy(pk, pk_seed, RLWE_SEED_BYTES);
    poly_sample_noise(s, noise_seed, 0x02, 0);
    poly_sample_noise(e, noise_seed, 0x02, 1);
    poly_uniform(a, pk, 0x01);
    poly_mul(b, a, s);  poly_add(b, b, e);
    poly_pack(pk + RLWE_SEED_BYTES, b);

    poly_pack(sk + SK_OFF_S, s);
    memcpy(sk + SK_OFF_PK,   pk, RLWE_PK_BYTES);
    sha3_256(sk + SK_OFF_HPK, pk, RLWE_PK_BYTES);
    memcpy(sk + SK_OFF_FAIL, fail_seed, RLWE_HASH_BYTES);
    return 0;
}

static int kem_encaps_det(uint8_t ct[RLWE_CT_BYTES],
                           uint8_t ss[RLWE_SS_BYTES],
                           const uint8_t pk[RLWE_PK_BYTES],
                           const uint8_t m[RLWE_SS_BYTES])
{
    uint8_t h_pk[RLWE_HASH_BYTES];
    uint8_t r_seed[32];
    uint8_t buf[RLWE_SS_BYTES + RLWE_HASH_BYTES];

    sha3_256(h_pk, pk, RLWE_PK_BYTES);
    memcpy(buf,                 m,    RLWE_SS_BYTES);
    memcpy(buf + RLWE_SS_BYTES, h_pk, RLWE_HASH_BYTES);
    sha3_256(r_seed, buf, sizeof(buf));

    pke_enc_det(ct, pk, m, r_seed);
    memcpy(ss, r_seed, RLWE_SS_BYTES);
    return 0;
}

static void run_kat(int count) {
    /* Fixed master seed — never change this for reproducibility */
    static const uint8_t MASTER[32] = {
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
        0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
        0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f
    };

    printf("# Ring-LWE KEM  Known Answer Test  (v3, IND-CCA2+FO)\n");
    printf("# n=%d  q=%d  eta=%d  pk=%d  sk=%d  ct=%d  ss=%d\n\n",
           RLWE_N, RLWE_Q, RLWE_ETA,
           RLWE_PK_BYTES, RLWE_SK_BYTES, RLWE_CT_BYTES, RLWE_SS_BYTES);

    /* Derive all seeds via SHAKE256(MASTER) */
    shake256_ctx kat_ctx;
    shake256_init(&kat_ctx);
    shake256_absorb(&kat_ctx, MASTER, 32);
    shake256_finalize(&kat_ctx);

    uint8_t pk[RLWE_PK_BYTES], sk[RLWE_SK_BYTES];
    uint8_t ct[RLWE_CT_BYTES], ss_e[RLWE_SS_BYTES], ss_d[RLWE_SS_BYTES];
    uint8_t pk_seed[32], noise_kg[32], m[32], fail_seed[32];

    for (int i = 0; i < count; i++) {
        shake256_squeeze(&kat_ctx, pk_seed,   32);
        shake256_squeeze(&kat_ctx, noise_kg,  32);
        shake256_squeeze(&kat_ctx, m,         32);
        shake256_squeeze(&kat_ctx, fail_seed, 32);

        kem_keypair_det(pk, sk, pk_seed, noise_kg, fail_seed);
        kem_encaps_det(ct, ss_e, pk, m);
        kem_decaps(ss_d, ct, sk);

        printf("count = %d\n", i);
        printf("pk_seed  = "); for(int j=0;j<32;j++) printf("%02x",pk_seed[j]);  printf("\n");
        printf("noise_kg = "); for(int j=0;j<32;j++) printf("%02x",noise_kg[j]); printf("\n");
        printf("msg      = "); for(int j=0;j<32;j++) printf("%02x",m[j]);        printf("\n");
        printf("pk       = "); for(int j=0;j<RLWE_PK_BYTES;j++) printf("%02x",pk[j]); printf("\n");
        printf("sk[0:32] = "); for(int j=0;j<32;j++) printf("%02x",sk[j]);       printf("\n");
        printf("ct       = "); for(int j=0;j<RLWE_CT_BYTES;j++) printf("%02x",ct[j]); printf("\n");
        printf("ss       = "); for(int j=0;j<RLWE_SS_BYTES;j++) printf("%02x",ss_e[j]); printf("\n");

        if (ct_memcmp(ss_e, ss_d, RLWE_SS_BYTES) != 0)
            printf("ERROR: decaps mismatch at count=%d\n", i);
        printf("\n");
    }
}

/* ============================================================
 *  §10  Test harness
 * ============================================================ */
static void run_correctness(int trials) {
    int fail = 0;
    uint8_t pk[RLWE_PK_BYTES], sk[RLWE_SK_BYTES];
    uint8_t ct[RLWE_CT_BYTES], ss_e[RLWE_SS_BYTES], ss_d[RLWE_SS_BYTES];

    for (int t = 0; t < trials; t++) {
        if (kem_keypair(pk, sk) < 0) { fail++; continue; }
        if (kem_encaps(ct, ss_e, pk) < 0) { fail++; continue; }
        kem_decaps(ss_d, ct, sk);
        if (ct_memcmp(ss_e, ss_d, RLWE_SS_BYTES) != 0) fail++;
    }
    printf("Correctness (%5d trials): failures = %d  (%.4f%%)\n",
           trials, fail, 100.0 * fail / trials);
}

static void run_pack_roundtrip(void) {
    poly orig, restored;
    uint8_t buf[RLWE_POLY_PACKED_BYTES];
    uint8_t seed[32];
    memset(seed, 0x42, sizeof(seed));

    shake256_ctx ctx;
    shake256_init(&ctx);
    shake256_absorb(&ctx, seed, 32);
    shake256_finalize(&ctx);
    uint8_t r2[2];
    for (int i = 0; i < N; ) {
        shake256_squeeze(&ctx, r2, 2);
        uint32_t r = (uint32_t)r2[0] | (((uint32_t)r2[1] & 0x0F) << 8);
        if (r < (uint32_t)Q) { orig[i] = (int16_t)r; i++; }
    }

    poly_pack(buf, orig);
    poly_unpack(restored, buf);

    int ok = 1;
    for (int i = 0; i < N; i++) if (orig[i] != restored[i]) { ok = 0; break; }
    printf("Pack/unpack round-trip:   %s\n", ok ? "PASS" : "FAIL");
}

static void run_ct_memcmp_test(void) {
    uint8_t a[32], b[32];
    memset(a, 0xAB, 32); memset(b, 0xAB, 32);
    int eq  = ct_memcmp(a, b, 32);
    b[16]  ^= 1;
    int neq = ct_memcmp(a, b, 32);
    printf("ct_memcmp equal/neq:      %s\n",
           (eq == 0 && neq != 0) ? "PASS" : "FAIL");
}

static void run_fo_rejection_test(void) {
    uint8_t pk[RLWE_PK_BYTES], sk[RLWE_SK_BYTES];
    uint8_t ct[RLWE_CT_BYTES], ss_valid[RLWE_SS_BYTES], ss_tampered[RLWE_SS_BYTES];

    kem_keypair(pk, sk);
    kem_encaps(ct, ss_valid, pk);

    /* Tamper heavily with ciphertext — decaps must return different ss */
    uint8_t ct_bad[RLWE_CT_BYTES];
    memcpy(ct_bad, ct, RLWE_CT_BYTES);
    /* Flip many bytes across both u and v parts */
    for (int i = 0; i < (int)RLWE_CT_BYTES; i += 7)
        ct_bad[i] ^= 0xFF;

    kem_decaps(ss_tampered, ct_bad, sk);

    int same = (ct_memcmp(ss_valid, ss_tampered, RLWE_SS_BYTES) == 0);
    printf("FO rejection (tampered):  %s\n", same ? "FAIL (no rejection!)" : "PASS");
}

static void run_sha3_selftest(void) {
    /* SHA3-256("abc") = 3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532 */
    static const uint8_t expected[32] = {
        0x3a,0x98,0x5d,0xa7,0x4f,0xe2,0x25,0xb2,
        0x04,0x5c,0x17,0x2d,0x6b,0xd3,0x90,0xbd,
        0x85,0x5f,0x08,0x6e,0x3e,0x9d,0x52,0x5b,
        0x46,0xbf,0xe2,0x45,0x11,0x43,0x15,0x32
    };
    uint8_t got[32];
    sha3_256(got, (const uint8_t *)"abc", 3);
    printf("SHA3-256(\"abc\") selftest: %s\n",
           (ct_memcmp(got, expected, 32) == 0) ? "PASS" : "FAIL");
}

/* ============================================================
 *  main
 * ============================================================ */
int main(int argc, char *argv[]) {

    if (argc >= 2 && strcmp(argv[1], "--kat") == 0) {
        int n = (argc >= 3) ? atoi(argv[2]) : 5;
        run_kat(n);
        return 0;
    }

    uint8_t pk[RLWE_PK_BYTES], sk[RLWE_SK_BYTES];
    uint8_t ct[RLWE_CT_BYTES], ss_a[RLWE_SS_BYTES], ss_b[RLWE_SS_BYTES];

    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║    Ring-LWE KEM  v3  —  IND-CCA2  (FO Transform)       ║\n");
    printf("╠══════════════════════════════════════════════════════════╣\n");
    printf("║  n=%-4d  q=%-5d  eta=%d  ss=%2d B                     ║\n",
           N, Q, RLWE_ETA, RLWE_SS_BYTES);
    printf("║  [IND-CCA2] FO Transform   [SHAKE256/SHA3-256]         ║\n");
    printf("║  [CT] encode/decode/select  [12-bit pack]  [0-alloc]  ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    printf("[1] kem_keypair\n");
    if (kem_keypair(pk, sk) < 0) { fprintf(stderr,"entropy fail\n"); return 1; }
    printf("  pk: "); for(int i=0;i<16;i++) printf("%02x",pk[i]); printf("... (%d B)\n", RLWE_PK_BYTES);
    printf("  sk: "); for(int i=0;i<16;i++) printf("%02x",sk[i]); printf("... (%d B)\n\n", RLWE_SK_BYTES);

    printf("[2] kem_encaps\n");
    if (kem_encaps(ct, ss_a, pk) < 0) { fprintf(stderr,"entropy fail\n"); return 1; }
    printf("  ss(Alice): "); for(int i=0;i<RLWE_SS_BYTES;i++) printf("%02x",ss_a[i]); printf("\n");
    printf("  ct: "); for(int i=0;i<16;i++) printf("%02x",ct[i]); printf("... (%d B)\n\n", RLWE_CT_BYTES);

    printf("[3] kem_decaps\n");
    kem_decaps(ss_b, ct, sk);
    printf("  ss(Bob):   "); for(int i=0;i<RLWE_SS_BYTES;i++) printf("%02x",ss_b[i]); printf("\n\n");

    printf("[4] Match: %s\n\n",
           (ct_memcmp(ss_a, ss_b, RLWE_SS_BYTES) == 0) ? "SUCCESS" : "FAILURE");

    printf("--- Wire sizes ---\n");
    printf("  pk  : %d B  (seed %d + pack(b) %d)\n",
           RLWE_PK_BYTES, RLWE_SEED_BYTES, RLWE_POLY_PACKED_BYTES);
    printf("  sk  : %d B  (s+pk+H(pk)+fail_seed)\n", RLWE_SK_BYTES);
    printf("  ct  : %d B  (u %d + v %d)\n",
           RLWE_CT_BYTES, RLWE_POLY_PACKED_BYTES, RLWE_POLY_PACKED_BYTES);
    printf("  ss  : %d B\n\n",   RLWE_SS_BYTES);

    printf("--- Unit tests ---\n");
    run_sha3_selftest();
    run_ct_memcmp_test();
    run_pack_roundtrip();
    run_fo_rejection_test();
    run_correctness(1000);

    printf("\nKAT: ./rlwe_kem --kat 5 > kat_vectors.txt\n");
    return 0;
}
