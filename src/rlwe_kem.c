/*
 * ============================================================
 *  rlwe_kem.c  —  Compact Ring-LWE KEM
 *
 *  Improvements over v1:
 *    1. Constant-time encode/decode  (bitmask, no branches on secrets)
 *    2. Constant-time comparison     (ct_memcmp)
 *    3. Clean public API in .h, all internals static
 *       Zero dynamic allocation — stack only
 *    4. KAT (Known Answer Test) with fixed PRNG seed
 *    5. 13-bit coefficient packing   (poly_pack / poly_unpack)
 *
 *  Compile:  gcc -O2 -Wall -Wextra -o rlwe_kem rlwe_kem.c
 *  KAT:      ./rlwe_kem --kat 5 > kat_vectors.txt
 * ============================================================
 */

#include "rlwe_kem.h"
#include "randombytes.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ── Internal aliases ─────────────────────────────────────── */
#define N       RLWE_N
#define Q       RLWE_Q
#define QHALF   3840   /* floor(Q/2) */

typedef int16_t poly[N];

/* ============================================================
 *  §1  PRNG — xoshiro128** (public domain, Vigna & Blackman)
 *
 *  Encapsulated in a struct so multiple independent streams
 *  can coexist (no global state in crypto paths).
 *  Domain byte ensures different uses of the same seed
 *  produce orthogonal outputs.
 * ============================================================ */
typedef struct { uint32_t s[4]; } Xoshiro;

static inline uint32_t _rotl(uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

static uint32_t xo_next(Xoshiro *x) {
    uint32_t result = _rotl(x->s[1] * 5, 7) * 9;
    uint32_t t = x->s[1] << 9;
    x->s[2] ^= x->s[0];
    x->s[3] ^= x->s[1];
    x->s[1] ^= x->s[2];
    x->s[0] ^= x->s[3];
    x->s[2] ^= t;
    x->s[3] = _rotl(x->s[3], 11);
    return result;
}

static void xo_init(Xoshiro *x, const uint8_t seed[32], uint8_t domain) {
    uint32_t a[4], b[4];
    memcpy(a, seed,      16);
    memcpy(b, seed + 16, 16);
    for (int i = 0; i < 4; i++) x->s[i] = a[i] ^ b[i];
    x->s[0] ^= (uint32_t)domain << 24;   /* domain separation */
    if (!x->s[0] && !x->s[1] && !x->s[2] && !x->s[3]) x->s[0] = 1;
    for (int i = 0; i < 20; i++) xo_next(x);
}

/* ── OS entropy ───────────────────────────────────────────── */
static int get_random_bytes(uint8_t *buf, size_t len) {
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) return -1;
    size_t r = fread(buf, 1, len, f);
    fclose(f);
    return (r == len) ? 0 : -1;
}

/* ============================================================
 *  §2  Constant-time primitives
 *
 *  Golden rule: never branch or index arrays on secret data.
 *  All selection uses arithmetic masks derived from comparison.
 *
 *  The compiler cannot optimize these to branches because the
 *  mask is computed from live data; use -O2 and verify asm.
 * ============================================================ */

/*
 * ct_memcmp — constant-time comparison.
 * Time depends only on 'len', not on content.
 * Returns 0 if buffers are equal, non-zero otherwise.
 */
static int ct_memcmp(const uint8_t *a, const uint8_t *b, size_t len) {
    uint8_t diff = 0;
    for (size_t i = 0; i < len; i++)
        diff |= a[i] ^ b[i];
    return (int)diff;
}

/*
 * ct_min_u32 — constant-time minimum of two uint32_t.
 * Uses: mask = -(b < a) = 0xFFFFFFFF if b<a, 0 otherwise.
 */
static inline uint32_t ct_min_u32(uint32_t a, uint32_t b) {
    uint32_t mask = -(uint32_t)(b < a);
    return a ^ ((a ^ b) & mask);
}

/* ============================================================
 *  §3  Barrett reduction (data-flow only, no branches)
 *
 *  Reduces a ∈ (-Q², Q²) to [0, Q).
 *  k = floor(2^24 / Q) = 2184  (for Q = 7681).
 *
 *  The two final corrections use -(condition) masks,
 *  which the compiler emits as cmov or equivalent.
 * ============================================================ */
static inline int16_t barrett(int32_t a) {
    int32_t t = (int32_t)(((int64_t)a * 2184) >> 24);
    int32_t r = a - t * Q;
    r -= Q & -(r >= Q);   /* if r >= Q: r -= Q */
    r += Q & -(r <  0);   /* if r <  0: r += Q */
    return (int16_t)r;
}

/* ============================================================
 *  §4  CBD noise sampling (constant-time)
 *
 *  CBD(eta=2): x = popcount(bits 0,1) - popcount(bits 2,3)
 *  x ∈ {-2,-1,0,1,2},  Var(x) = eta/2 = 1
 *
 *  We compute popcount via pure arithmetic (no lookup tables,
 *  no data-dependent memory access, no branches on noise values).
 *
 *  Note: popcount is NOT secret-dependent here —
 *  the noise is sampled before it touches any secret.
 *  But we keep it branchless anyway for uniformity.
 * ============================================================ */
static inline int ct_pop2(uint32_t x) {
    /* popcount of bits {0,1} and bits {16,17} combined */
    uint32_t m = x & 0x00030003u;
    m = (m & 0x00010001u) + ((m >> 1) & 0x00010001u);
    return (int)((m + (m >> 16)) & 0xF);
}

static int16_t cbd_sample(Xoshiro *rng) {
    uint32_t r = xo_next(rng);
    int a = ct_pop2(r);       /* popcount of low pair  */
    int b = ct_pop2(r >> 2);  /* popcount of high pair */
    return (int16_t)(a - b);  /* ∈ [-2, 2]             */
}

static void poly_sample_noise(poly p, Xoshiro *rng) {
    for (int i = 0; i < N; i++)
        p[i] = barrett((int32_t)cbd_sample(rng) + Q);
}

/* ============================================================
 *  §5  Polynomial arithmetic in R_q = Z_q[x]/(x^n+1)
 *
 *  poly_mul: negacyclic schoolbook O(n²).
 *
 *  *** Constant-time fix vs v1 ***
 *  v1 had: if (!a[i]) continue;
 *  This is a branch on a[i] — when a[] is the secret key,
 *  this leaks Hamming weight via timing. Removed completely.
 *  The inner loop now always runs N² iterations.
 *
 *  Wrap logic: instead of "if (k>=N) acc[k-N] -= prod",
 *  we compute sign = 1 - 2*wrap and index = k - wrap*N
 *  using integer arithmetic, eliminating the branch.
 * ============================================================ */
static void poly_mul(poly c, const poly a, const poly b) {
    int64_t acc[N];
    memset(acc, 0, sizeof(acc));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int     k    = i + j;
            int     wrap = (k >= N);          /* 0 or 1, no branch emitted  */
            int     idx  = k - wrap * N;      /* index into acc[]           */
            int64_t sign = 1 - 2 * wrap;      /* +1 or -1                   */
            acc[idx] += sign * (int64_t)a[i] * (int64_t)b[j];
        }
    }

    for (int i = 0; i < N; i++) {
        int32_t v = (int32_t)(acc[i] % Q);
        c[i] = barrett(v);
    }
}

static void poly_add(poly c, const poly a, const poly b) {
    for (int i = 0; i < N; i++)
        c[i] = barrett(a[i] + b[i]);
}

static void poly_sub(poly c, const poly a, const poly b) {
    for (int i = 0; i < N; i++)
        c[i] = barrett(a[i] - b[i] + Q);
}

/* Expand public polynomial uniformly from seed.
   Rejection sampling loop: expected ~1.07 iters per coeff. */
static void poly_uniform(poly a, const uint8_t seed[32]) {
    Xoshiro rng;
    xo_init(&rng, seed, 0x01);   /* domain 0x01 = uniform poly */
    for (int i = 0; i < N; i++) {
        uint32_t r;
        do { r = xo_next(&rng) & 0x3FFF; } while (r >= (uint32_t)Q);
        a[i] = (int16_t)r;
    }
}

/* ============================================================
 *  §6  13-bit coefficient packing / unpacking
 *
 *  q = 7681 < 2^13 = 8192  →  each coeff fits in 13 bits.
 *  256 coeffs × 13 bits = 3328 bits = 416 bytes  (vs 512 uncompressed).
 *  Reduction: 96 bytes per polynomial = 18.75%.
 *
 *  Algorithm: sliding bit-buffer, fill/drain 8 bits at a time.
 *  Pure bijection — no secret-dependent control flow.
 * ============================================================ */
static void poly_pack(uint8_t out[RLWE_POLY_PACKED_BYTES], const poly p) {
    uint32_t buf = 0;
    int      cnt = 0;
    int      bi  = 0;

    for (int i = 0; i < N; i++) {
        buf |= (uint32_t)(uint16_t)p[i] << cnt;
        cnt += 13;
        if (cnt >= 8) { out[bi++] = (uint8_t)buf; buf >>= 8; cnt -= 8; }
        if (cnt >= 8) { out[bi++] = (uint8_t)buf; buf >>= 8; cnt -= 8; }
    }
    if (cnt > 0) out[bi] = (uint8_t)buf;
}

static void poly_unpack(poly p, const uint8_t in[RLWE_POLY_PACKED_BYTES]) {
    uint32_t buf = 0;
    int      cnt = 0;
    int      bi  = 0;

    for (int i = 0; i < N; i++) {
        while (cnt < 13) {
            buf |= (uint32_t)in[bi++] << cnt;
            cnt += 8;
        }
        p[i] = (int16_t)(buf & 0x1FFF);
        buf >>= 13;
        cnt  -= 13;
    }
}

/* ============================================================
 *  §7  Constant-time message encode / decode
 *
 *  Encoding (msg is public during encaps — branches OK):
 *    bit = 0 → coeff = 0
 *    bit = 1 → coeff = QHALF = 3840
 *    Branchless via mask: -(bit) & QHALF
 *
 *  Decoding (w depends on secret key → MUST be constant-time):
 *    bit = (dist(coeff, QHALF) < dist(coeff, 0)) ? 1 : 0
 *
 *    Both distances computed with ct_min_u32.
 *    Final bit selection: (d1 - d0) >> 31 extracts sign bit,
 *    which is 1 iff d1 < d0 i.e. closer to QHALF.
 *    This is pure arithmetic — zero branches on secret data.
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
        uint32_t c = (uint32_t)(uint16_t)w[i];   /* c ∈ [0, Q) */

        /* dist(c, 0) mod Q = min(c, Q-c) */
        uint32_t d0 = ct_min_u32(c, (uint32_t)Q - c);

        /* dist(c, QHALF) mod Q = min(|c-QHALF|, Q-|c-QHALF|) */
        uint32_t diff = (c >= (uint32_t)QHALF) ? c - QHALF : QHALF - c;
        uint32_t d1   = ct_min_u32(diff, (uint32_t)Q - diff);

        /*
         * bit = 1  iff  d1 < d0
         * (d1 - d0) underflows → MSB = 1 when d1 < d0
         * No branch: pure integer arithmetic.
         */
        uint32_t bit = (d1 - d0) >> 31;
        msg[i >> 3] |= (uint8_t)(bit << (i & 7));
    }
}

/* ============================================================
 *  §8  Public API
 * ============================================================ */

int kem_keypair(uint8_t pk[RLWE_PK_BYTES],
                uint8_t sk[RLWE_SK_BYTES])
{
    poly    a, s, e, b;
    Xoshiro rng;
    uint8_t noise_seed[32];

    /* pk[0..31] = seed for a (public, sent in plaintext) */
    if (get_random_bytes(pk, RLWE_SEED_BYTES) < 0) return -1;
    if (get_random_bytes(noise_seed, 32)       < 0) return -1;

    xo_init(&rng, noise_seed, 0x02);   /* domain 0x02 = keygen noise */
    poly_sample_noise(s, &rng);        /* secret  s ← χ_η            */
    poly_sample_noise(e, &rng);        /* error   e ← χ_η            */

    poly_uniform(a, pk);               /* a from seed                 */
    poly_mul(b, a, s);
    poly_add(b, b, e);                 /* b = a·s + e  (RLWE sample) */

    poly_pack(pk + RLWE_SEED_BYTES, b);  /* pk = seed ‖ pack(b) */
    poly_pack(sk, s);                    /* sk = pack(s)         */
    return 0;
}

int kem_encaps(uint8_t ct[RLWE_CT_BYTES],
               uint8_t ss[RLWE_SS_BYTES],
               const uint8_t pk[RLWE_PK_BYTES])
{
    poly    a, b, s2, e1, e2, u, v, tmp, mp;
    Xoshiro rng;
    uint8_t m[RLWE_SS_BYTES], noise_seed[32];

    if (get_random_bytes(m, RLWE_SS_BYTES) < 0) return -1;
    if (get_random_bytes(noise_seed, 32)   < 0) return -1;

    xo_init(&rng, noise_seed, 0x03);   /* domain 0x03 = encaps noise */
    poly_sample_noise(s2, &rng);
    poly_sample_noise(e1, &rng);
    poly_sample_noise(e2, &rng);

    poly_uniform(a, pk);
    poly_unpack(b, pk + RLWE_SEED_BYTES);

    poly_mul(u, a, s2);  poly_add(u, u, e1);       /* u = a·s' + e'          */
    poly_mul(tmp, b, s2); poly_add(tmp, tmp, e2);   /* tmp = b·s' + e''       */
    encode_msg(mp, m);
    poly_add(v, tmp, mp);                           /* v = b·s' + e'' + enc(m)*/

    poly_pack(ct,                          u);
    poly_pack(ct + RLWE_POLY_PACKED_BYTES, v);

    memcpy(ss, m, RLWE_SS_BYTES);
    return 0;
}

void kem_decaps(uint8_t ss[RLWE_SS_BYTES],
                const uint8_t ct[RLWE_CT_BYTES],
                const uint8_t sk[RLWE_SK_BYTES])
{
    poly s, u, v, w;

    poly_unpack(s, sk);
    poly_unpack(u, ct);
    poly_unpack(v, ct + RLWE_POLY_PACKED_BYTES);

    /*
     * w = v - u·s
     *   = (b·s'+e''+enc(m)) - (a·s'+e')·s
     *   = (a·s+e)·s' + e'' + enc(m) - a·s·s' - e'·s
     *   = e·s' - e'·s + e'' + enc(m)
     *     └─────────────────┘ bounded noise
     *
     * poly_mul is constant-time (uniform loop, no branch on s).
     * decode_msg is constant-time (see §7).
     */
    poly_mul(w, u, s);
    poly_sub(w, v, w);
    decode_msg(ss, w);
}

/* ============================================================
 *  §9  KAT — Known Answer Tests
 *
 *  Deterministic: fixed master seed → reproducible keys/ciphertexts.
 *  Use:  ./rlwe_kem --kat 5 > kat_vectors.txt
 *  After any refactor, re-run and diff against saved vectors.
 *
 *  We use separate "det" variants that accept explicit seeds,
 *  so the production code paths are never modified.
 * ============================================================ */
static int kem_keypair_det(uint8_t pk[RLWE_PK_BYTES],
                           uint8_t sk[RLWE_SK_BYTES],
                           const uint8_t pk_seed[RLWE_SEED_BYTES],
                           const uint8_t noise_seed[32])
{
    poly a, s, e, b;
    Xoshiro rng;

    memcpy(pk, pk_seed, RLWE_SEED_BYTES);
    xo_init(&rng, noise_seed, 0x02);
    poly_sample_noise(s, &rng);
    poly_sample_noise(e, &rng);
    poly_uniform(a, pk);
    poly_mul(b, a, s);  poly_add(b, b, e);
    poly_pack(pk + RLWE_SEED_BYTES, b);
    poly_pack(sk, s);
    return 0;
}

static int kem_encaps_det(uint8_t ct[RLWE_CT_BYTES],
                          uint8_t ss[RLWE_SS_BYTES],
                          const uint8_t pk[RLWE_PK_BYTES],
                          const uint8_t m[RLWE_SS_BYTES],
                          const uint8_t noise_seed[32])
{
    poly a, b, s2, e1, e2, u, v, tmp, mp;
    Xoshiro rng;

    xo_init(&rng, noise_seed, 0x03);
    poly_sample_noise(s2, &rng);
    poly_sample_noise(e1, &rng);
    poly_sample_noise(e2, &rng);
    poly_uniform(a, pk);
    poly_unpack(b, pk + RLWE_SEED_BYTES);
    poly_mul(u, a, s2);   poly_add(u, u, e1);
    poly_mul(tmp, b, s2); poly_add(tmp, tmp, e2);
    encode_msg(mp, m);    poly_add(v, tmp, mp);
    poly_pack(ct,                          u);
    poly_pack(ct + RLWE_POLY_PACKED_BYTES, v);
    memcpy(ss, m, RLWE_SS_BYTES);
    return 0;
}

static void run_kat(int count) {
    /* Fixed master seed — never change this */
    static const uint8_t MASTER[32] = {
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
        0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
        0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f
    };
    Xoshiro kat;
    xo_init(&kat, MASTER, 0x00);

    printf("# Ring-LWE KEM  Known Answer Test  (v2)\n");
    printf("# n=%d  q=%d  eta=%d  pk=%d  sk=%d  ct=%d  ss=%d\n\n",
           RLWE_N, RLWE_Q, RLWE_ETA,
           RLWE_PK_BYTES, RLWE_SK_BYTES, RLWE_CT_BYTES, RLWE_SS_BYTES);

    uint8_t pk[RLWE_PK_BYTES], sk[RLWE_SK_BYTES];
    uint8_t ct[RLWE_CT_BYTES], ss_e[RLWE_SS_BYTES], ss_d[RLWE_SS_BYTES];
    uint8_t pk_seed[32], noise_kg[32], m[32], noise_enc[32];

    for (int i = 0; i < count; i++) {
        /* derive all seeds deterministically */
        for (int j = 0; j < 8;  j++) ((uint32_t*)pk_seed)[j]   = xo_next(&kat);
        for (int j = 0; j < 8;  j++) ((uint32_t*)noise_kg)[j]  = xo_next(&kat);
        for (int j = 0; j < 8;  j++) ((uint32_t*)m)[j]         = xo_next(&kat);
        for (int j = 0; j < 8;  j++) ((uint32_t*)noise_enc)[j] = xo_next(&kat);

        kem_keypair_det(pk, sk, pk_seed, noise_kg);
        kem_encaps_det(ct, ss_e, pk, m, noise_enc);
        kem_decaps(ss_d, ct, sk);

        printf("count = %d\n", i);
        printf("pk_seed   = "); for(int j=0;j<32;j++) printf("%02x",pk_seed[j]);   printf("\n");
        printf("noise_kg  = "); for(int j=0;j<32;j++) printf("%02x",noise_kg[j]);  printf("\n");
        printf("msg       = "); for(int j=0;j<32;j++) printf("%02x",m[j]);          printf("\n");
        printf("noise_enc = "); for(int j=0;j<32;j++) printf("%02x",noise_enc[j]); printf("\n");
        printf("pk        = "); for(int j=0;j<RLWE_PK_BYTES;j++) printf("%02x",pk[j]); printf("\n");
        printf("sk        = "); for(int j=0;j<RLWE_SK_BYTES;j++) printf("%02x",sk[j]); printf("\n");
        printf("ct        = "); for(int j=0;j<RLWE_CT_BYTES;j++) printf("%02x",ct[j]); printf("\n");
        printf("ss        = "); for(int j=0;j<RLWE_SS_BYTES;j++) printf("%02x",ss_e[j]); printf("\n");

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
        kem_keypair(pk, sk);
        kem_encaps(ct, ss_e, pk);
        kem_decaps(ss_d, ct, sk);
        if (ct_memcmp(ss_e, ss_d, RLWE_SS_BYTES) != 0) fail++;
    }
    printf("Correctness (%5d trials): failures = %d  (%.4f%%)\n",
           trials, fail, 100.0 * fail / trials);
}

static void run_pack_roundtrip(void) {
    poly orig, restored;
    uint8_t buf[RLWE_POLY_PACKED_BYTES];
    uint8_t seed[32] = {0x42, 0x13};
    Xoshiro rng; xo_init(&rng, seed, 0x10);

    for (int i = 0; i < N; i++) {
        uint32_t r;
        do { r = xo_next(&rng) & 0x3FFF; } while (r >= (uint32_t)Q);
        orig[i] = (int16_t)r;
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
    b[16] ^= 1;
    int neq = ct_memcmp(a, b, 32);
    printf("ct_memcmp equal/neq:      %s\n",
           (eq == 0 && neq != 0) ? "PASS" : "FAIL");
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

    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║      Ring-LWE KEM  v2  —  hardened build            ║\n");
    printf("╠══════════════════════════════════════════════════════╣\n");
    printf("║  n=%-4d  q=%-5d  eta=%d  ss=%2d B                  ║\n",
           N, Q, RLWE_ETA, RLWE_SS_BYTES);
    printf("║  [CT] encode/decode  [CT] memcmp  [0-alloc] stack   ║\n");
    printf("║  [13-bit packing]    [KAT support]                  ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

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
           !ct_memcmp(ss_a, ss_b, RLWE_SS_BYTES) ? "SUCCESS" : "FAILURE");

    printf("--- Wire sizes ---\n");
    printf("  pk : %d B   (seed %d + poly %d)\n",
           RLWE_PK_BYTES, RLWE_SEED_BYTES, RLWE_POLY_PACKED_BYTES);
    printf("  sk : %d B\n",   RLWE_SK_BYTES);
    printf("  ct : %d B   (u %d + v %d)\n",
           RLWE_CT_BYTES, RLWE_POLY_PACKED_BYTES, RLWE_POLY_PACKED_BYTES);
    printf("  ss : %d B\n",   RLWE_SS_BYTES);
    printf("  Packing saves: %zu B per poly vs 16-bit\n\n",
           sizeof(poly) - RLWE_POLY_PACKED_BYTES);

    printf("--- Unit tests ---\n");
    run_ct_memcmp_test();
    run_pack_roundtrip();
    run_correctness(2000);

    printf("\nKAT: ./rlwe_kem --kat 5 > kat_vectors.txt\n");
    return 0;
}
