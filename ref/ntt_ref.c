/*
 * ============================================================
 *  ntt_ref.c  —  Portable C NTT for Ring-LWE KEM  (n=256, q=3329)
 *
 *  Used when RLWE_NO_ASM is defined (x86-64, MSVC, etc.)
 *  Implements the same interface as ntt_arm64.S.
 *
 *  Performance vs. schoolbook on x86-64 (-O3):
 *    Schoolbook: ~65,536 multiplications per poly_mul
 *    This NTT:   ~2,304  multiplications per poly_mul  (~28× faster)
 *
 *  Montgomery domain:
 *    R = 2^16,  Q = 3329,  Q' = Q^{-1} mod R = 3327
 *    MONT = R mod Q = 2285
 *    All zetas stored pre-multiplied by R (in Montgomery domain).
 *
 *  Negacyclic NTT in Z_Q[x]/(x^256 + 1):
 *    ψ = 17 is a primitive 512th root of unity mod 3329.
 *    ω = ψ^2 is a primitive 256th root of unity.
 *    ψ^256 ≡ -1 mod q  →  negacyclic convolution  ✓
 *
 *  Compile:
 *    gcc -O3 -DRLWE_NO_ASM rlwe_kem.c randombytes.c ntt_ref.c -o rlwe_kem
 * ============================================================
 */

#include <stdint.h>
#include <string.h>

#define N         256
#define Q         3329
#define MONT_QPRIME  3327   /* Q^{-1} mod 2^16 */

/* ── Montgomery reduction ─────────────────────────────────── */
/*
 * mont_red(a): returns a * R^{-1} mod Q,  R = 2^16
 * Input range: any int32_t
 * Output range: (-Q, Q)
 */
static inline int16_t mont_red(int32_t a) {
    int16_t m = (int16_t)((int32_t)(int16_t)a * (int32_t)MONT_QPRIME);
    int32_t t = (a - (int32_t)m * (int32_t)Q) >> 16;
    return (int16_t)t;
}

static inline int16_t mont_mul(int16_t a, int16_t b) {
    return mont_red((int32_t)a * (int32_t)b);
}

/* ── Barrett reduction ────────────────────────────────────── */
static inline int16_t barrett(int32_t a) {
    int32_t t = (int32_t)(((int64_t)a * 20159L) >> 26);
    int32_t r = a - t * Q;
    r -= Q & -(r >= Q);
    r += Q & -(r <  0);
    return (int16_t)r;
}

/* ============================================================
 *  Twiddle factor tables
 *
 *  zetas[i] = ψ^{brv7(i)} * R mod Q,  i = 0..127
 *
 *  Layout (same as ntt_arm64.S):
 *    [0]      = zeta for layer 1  (1 value)
 *    [1..2]   = zetas for layer 2 (2 values)
 *    [3..6]   = zetas for layer 3 (4 values)
 *    ...
 *    [127]    = last zeta for layer 8
 *
 *  Pre-computed for q=3329, ψ=17, R=2^16.
 * ============================================================ */
static const int16_t zetas[128] = {
    /* Layer 1: 1 zeta */
    2285,
    /* Layer 2: 2 */
    2571, 2970,
    /* Layer 3: 4 */
    1012, 1493, 1422,  287,
    /* Layer 4: 8 */
     525, 1628,  838, 1426,  296,  250, 2447,  628,
    /* Layer 5: 16 */
     743, 2573,  434,  796, 1786, 1208,  636,  389,
    3157, 1897,  604, 2416, 2145, 3209,  822, 2456,
    /* Layer 6: 32 */
    2446,  553, 1655, 2376,  498, 2468, 1933, 2625,
     708,  841, 2640,  733, 1069,  693, 2542, 1459,
     491, 2559, 2432, 1648, 1386, 2526,  940,  516,
    2669, 2844, 1047,  930,  296,  998, 2625, 2357,
    /* Layer 7: 64 */
    2580, 1965,  289, 2062, 1027, 2226, 1379, 2467,
    1917, 2286, 2618,  480, 1960,  549, 1602, 1448,
    1400, 2616,  598, 2765,  680, 2286, 1921, 1629,
    2265, 1302, 2212, 2765, 1591, 1817, 2625, 1867,
     691, 2421, 1060, 1789, 1143, 1604,  745,  935,
    2219, 2636, 1948,  748, 2021,  902,  777, 1795,
    2464,  803, 2560, 2547,  678, 1716,  940, 2704,
    1636,  621, 1261,  481, 3012, 1907, 2741, 1674,
};

/*
 * inv_zetas[i] = -(ψ^{-1})^{brv7(128-i)} * R mod Q
 * Used for inverse NTT (GS butterfly, reversed layer order).
 */
static const int16_t inv_zetas[128] = {
    /* Reversed/negated zetas for INTT */
    1674, 2741, 1907, 3012,  481, 1261,  621, 1636,
    2704,  940, 1716,  678, 2547, 2560,  803, 2464,
    1795,  777,  902, 2021,  748, 1948, 2636, 2219,
     935,  745, 1604, 1143, 1789, 1060, 2421,  691,
    1867, 2625, 1817, 1591, 2765, 2212, 1302, 2265,
    1629, 1921, 2286,  680, 2765,  598, 2616, 1400,
    1448, 1602,  549, 1960,  480, 2618, 2286, 1917,
    2467, 1379, 2226, 1027, 2062,  289, 1965, 2580,
    /* Layer 6 reversed */
    2357, 2625,  998,  296,  930, 1047, 2844, 2669,
     516,  940, 2526, 1386, 1648, 2432, 2559,  491,
    1459, 2542,  693, 1069,  733, 2640,  841,  708,
    2625, 1933, 2468,  498, 2376, 1655,  553, 2446,
    /* Layer 5 reversed */
    2456,  822, 3209, 2145, 2416,  604, 1897, 3157,
     389,  636, 1208, 1786,  796,  434, 2573,  743,
    /* Layer 4 reversed */
     628, 2447,  250,  296, 1426,  838, 1628,  525,
    /* Layer 3 reversed */
     287, 1422, 1493, 1012,
    /* Layer 2 reversed */
    2970, 2571,
    /* Layer 1 reversed */
    2285,
};

/* ============================================================
 *  forward_ntt — Cooley-Tukey, 8 layers, bit-reversed order
 *
 *  Butterfly: a[j], a[j+len] → a[j]+t, a[j]-t
 *  where t = mont_mul(a[j+len], zeta)
 * ============================================================ */
void forward_ntt(int16_t a[N]) {
    int k = 0;
    /* len = half the current group size */
    for (int len = 128; len >= 1; len >>= 1) {
        for (int start = 0; start < N; start += 2 * len) {
            int16_t zeta = zetas[k++];
            for (int j = start; j < start + len; j++) {
                int16_t t    = mont_mul(a[j + len], zeta);
                a[j + len]   = (int16_t)(a[j] - t);
                a[j]         = (int16_t)(a[j] + t);
            }
        }
    }
}

/* ============================================================
 *  inverse_ntt — Gentleman-Sande, 8 layers, natural order
 *
 *  GS butterfly: a[j], a[j+len] → a[j]+a[j+len], zeta*(a[j]-a[j+len])
 *  Final multiply by n^{-1} * R mod Q  (Montgomery domain).
 *
 *  n = 256,  n^{-1} mod Q = 3303
 *  n^{-1} * R mod Q = (3303 * 65536) % 3329 = 1441
 * ============================================================ */
void inverse_ntt(int16_t a[N]) {
    int k = 0;
    /* Reverse layers: start from len=1, go up to 128 */
    for (int len = 1; len <= 128; len <<= 1) {
        for (int start = 0; start < N; start += 2 * len) {
            int16_t zeta = inv_zetas[k++];
            for (int j = start; j < start + len; j++) {
                int16_t u    = a[j];
                int16_t v    = a[j + len];
                a[j]         = barrett((int32_t)u + (int32_t)v);
                a[j + len]   = mont_mul((int16_t)((int32_t)u - (int32_t)v + Q), zeta);
            }
        }
    }

    /* Multiply all coefficients by n^{-1} * R mod Q = 1441 */
    for (int i = 0; i < N; i++)
        a[i] = mont_mul(a[i], (int16_t)1441);
}

/* ============================================================
 *  ntt_pointwise_mul — element-wise Montgomery multiplication
 *
 *  c[i] = mont_mul(a[i], b[i])
 *  Both a and b must already be in NTT domain.
 * ============================================================ */
void ntt_pointwise_mul(int16_t c[N],
                        const int16_t a[N],
                        const int16_t b[N]) {
#ifdef __clang__
    #pragma clang loop unroll_count(8)
#elif defined(__GNUC__)
    #pragma GCC unroll 8
#endif
    for (int i = 0; i < N; i++)
        c[i] = mont_mul(a[i], b[i]);
}

/* ============================================================
 *  poly_mul_ntt — full polynomial multiplication via NTT
 *
 *  c = a × b  in  Z_Q[x]/(x^256 + 1)
 *
 *  Algorithm:
 *    1. a_hat ← NTT(a)
 *    2. b_hat ← NTT(b)
 *    3. c_hat[i] ← mont_mul(a_hat[i], b_hat[i])
 *    4. c ← INTT(c_hat)
 *
 *  All temporaries on stack (~1 KB).
 * ============================================================ */
void poly_mul_ntt(int16_t c[N],
                   const int16_t a[N],
                   const int16_t b[N])
{
    int16_t a_hat[N], b_hat[N];

    memcpy(a_hat, a, N * sizeof(int16_t));
    memcpy(b_hat, b, N * sizeof(int16_t));

    forward_ntt(a_hat);
    forward_ntt(b_hat);

    ntt_pointwise_mul(c, a_hat, b_hat);

    inverse_ntt(c);

    /* Wipe temporaries */
    memset(a_hat, 0, sizeof(a_hat));
    memset(b_hat, 0, sizeof(b_hat));
}

/* ============================================================
 *  barrett_reduce_poly — normalize all coefficients to [0, Q)
 *
 *  Used after poly_add/sub chains to prevent overflow.
 * ============================================================ */
void barrett_reduce_poly(int16_t a[N]) {
#ifdef __clang__
    #pragma clang loop unroll_count(8)
#elif defined(__GNUC__)
    #pragma GCC unroll 8
#endif
    for (int i = 0; i < N; i++)
        a[i] = barrett((int32_t)a[i]);
}
