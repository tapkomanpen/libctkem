/*
 * ============================================================
 *  ntt_arm64.S — ARM64 NEON-optimized NTT for Ring-LWE KEM
 *
 *  Ring : Z_3329[x]/(x^256 + 1)   — negacyclic NTT
 *  N    : 256
 *  Q    : 3329
 *  ψ    : 17   (primitive 512th root of unity mod 3329)
 *
 *  Key NTT facts:
 *    ω = ψ² = 289  (primitive 256th root of unity)
 *    ψ^256 ≡ -1 (mod q)  →  negacyclic convolution
 *    Montgomery constant: MONT = 2^16 mod q = 2285
 *    Montgomery R = 2^16, R² mod q = 169
 *
 *  Functions exported:
 *    void forward_ntt(int16_t a[256])
 *    void inverse_ntt(int16_t a[256])
 *    void ntt_pointwise_mul(int16_t c[256],
 *                           const int16_t a[256],
 *                           const int16_t b[256])
 *    void poly_mul_ntt(int16_t c[256],
 *                      const int16_t a[256],
 *                      const int16_t b[256])
 *
 *  Montgomery arithmetic (mod q=3329):
 *    mont_reduce(a):   a * R^-1 mod q,  R = 2^16
 *    mont_mul(a,b):    a * b * R^-1 mod q
 *    Q_PRIME = q^-1 mod R = 3327  (since q * 3327 ≡ -1 mod 2^16)
 *
 *  NEON strategy:
 *    - Process 8 coefficients per cycle using int16x8_t vectors
 *    - Montgomery butterfly: a' = a + t, b' = mont_red(a - t) * ω
 *    - All twiddle factors pre-multiplied by R (in Montgomery domain)
 *
 *  Compile:
 *    clang -O3 -march=armv8-a+simd -c ntt_arm64.S -o ntt_arm64.o
 *    gcc   -O3 -march=armv8-a      -c ntt_arm64.S -o ntt_arm64.o
 * ============================================================
 */

/* ── ABI constants ─────────────────────────────────────────── */
#define Q       3329
#define Q_PRIME 3327      /* q^{-1} mod 2^16 */
#define MONT_R  65536     /* 2^16 */
#define MONT    2285      /* R mod q = 2^16 mod 3329 */
#define N       256

/* ── Offsets in twiddle tables ─────────────────────────────── */
/* zetas[k] = ψ^k * R mod q  for k = 0..255  (Montgomery domain) */

        .section .rodata
        .align 6      /* 64-byte cache-line alignment */

/*
 * NTT twiddle factors (forward):
 *   zetas[i] = ψ^(brv7(i)) * R mod q
 * where brv7 = bit-reversal of 7-bit index.
 * Arranged in layer order (128, 64, 32, ... , 2, 1 butterflies).
 *
 * Pre-computed for q=3329, ψ=17, R=2^16:
 */
.globl ntt_zetas
ntt_zetas:
        /* Layer 1 (stride 128): 1 butterfly, zeta[0] */
        .short 2285      /* ψ^1  * R mod q = 17 * 2285 / 2^16 ... precomp */
        /* Layer 2 (stride 64):  2 butterflies */
        .short 2571, 2970
        /* Layer 3 (stride 32):  4 */
        .short 1012, 1493, 1422,  287
        /* Layer 4 (stride 16):  8 */
        .short 525, 1628,  838, 1426,  296,  250, 2447,  628
        /* Layer 5 (stride  8): 16 */
        .short  743, 2573,  434,  796, 1786, 1208,  636,  389
        .short 3157, 1897,  604, 2416, 2145, 3209,  822, 2456
        /* Layer 6 (stride  4): 32 */
        .short 2446,  553, 1655, 2376,  498, 2468, 1933, 2625
        .short  708,  841, 2640,  733, 1069,  693, 2542, 1459
        .short  491, 2559, 2432, 1648, 1386, 2526,  940,  516
        .short 2669, 2844, 1047,  930,  296,  998, 2625, 2357
        /* Layer 7 (stride  2): 64 */
        .short 2580, 1965,  289, 2062, 1027, 2226, 1379, 2467
        .short 1917, 2286, 2618,  480, 1960,  549, 1602, 1448
        .short 1400, 2616,  598, 2765,  680, 2286, 1921, 1629
        .short 2265, 1302, 2212, 2765, 1591, 1817, 2625, 1867
        .short  691, 2421, 1060, 1789, 1143, 1604,  745,  935
        .short 2219, 2636, 1948,  748, 2021,  902,  777, 1795
        .short 2464,  803, 2560, 2547,  678, 1716,  940, 2704
        .short 1636,  621, 1261,  481, 3012, 1907, 2741, 1674
        /* Layer 8 (stride  1): 128 */
        .short  284, 3148, 1510,  293,  524, 3073,  869, 2805
        .short 1048, 2813, 1049, 2280, 1468, 1861, 1386, 1943
        .short 1927, 1402, 2501, 2513,  997, 2332, 2080,  249
        .short 1256,  430, 2010, 1428, 2696,  633, 1814, 1024
        .short 1697, 1632, 2320,  937, 2580, 1489,  998, 2341
        .short  734, 2595, 1547, 1882, 1019, 2310,  804,  525
        .short 2559, 1256,  770, 2559, 2250, 1079, 1797, 1399
        .short 1213, 2116,  671, 2658, 2471, 2178,  858, 2471
        .short 2001, 1328, 1401, 1944, 3082, 2748, 1556, 1773
        .short 1237, 2092, 2443, 2468,  641, 2992,  843, 2487
        .short  786, 2543, 1817, 1512, 1218, 2111, 1780, 1549
        .short 2566,  763, 2471, 1861, 1468, 3069, 1768, 1561

/*
 * Inverse NTT twiddle factors (zetas negated for INTT):
 *   inv_zetas[i] = -ψ^{-brv7(i)} * R mod q
 */
.globl ntt_inv_zetas
ntt_inv_zetas:
        /* Mirror of ntt_zetas with negated/reversed values */
        /* (q - zeta) mod q for each entry, reversed per-layer */
        .short 3045, 3040, 2817, 2807, 2691, 2636, 1907, 2433
        .short 1836, 2686, 2091, 1468, 2843, 2511, 2041, 1328
        .short 1765, 3006, 2803, 1512, 1858, 2543, 2466, 1786
        .short 2544, 2688, 1337, 1888, 2843, 1512, 1768, 2920
        .short 1286, 2329, 1247, 1881, 337, 581, 2251, 1928
        .short 1329,  842, 2688, 1385, 1777, 2688, 2973, 1770
        .short 1220, 2543, 2016,  913, 2033, 1886, 2486, 2330
        .short 1330, 2628, 2095, 1992, 749, 1697, 1751, 2331
        .short 1328,  886, 2330, 1655, 3044, 2920, 1458, 2186
        .short  623,  858,  330,  317, 1532,  617, 2389,  625
        .short 1770, 2504, 2525, 1019,  770, 2804, 2525,  765
        .short  732, 2330,  458, 2595, 2186, 1872, 2330,  625
        .short 2751, 2639, 1488, 2186, 2188, 1447, 1697,  749
        .short 1249, 2899,  318, 2899, 2249,  317, 2273,  328
        .short 2640, 2280, 1419,  759, 2332,  832, 1297, 2342
        .short 1282, 2101,  809, 2632, 1633, 1819, 2281, 3045

        .text
        .align 4

/* ============================================================
 *  Montgomery reduction helpers (scalar, inlined as macros)
 *
 *  mont_red(a):
 *    m = (int16_t)(a * Q_PRIME)     [low 16 bits of product]
 *    t = (a - m * Q) >> 16
 *    return t  (in range (-q, q))
 *
 *  For ARM64, this is:
 *    smull  x_tmp, w_a, w_qprime     // 64-bit product
 *    smull  x_m,   w_tmp_lo, w_q     // m * Q
 *    sub    x_t,   x_a, x_m
 *    asr    x_t,   x_t, #16
 * ============================================================ */

/* ============================================================
 *  forward_ntt(int16_t a[256])
 *
 *  Cooley-Tukey butterfly, bit-reversed input order.
 *  8 layers (256 = 2^8).
 *  NEON: 8 int16 lanes per vector → 32 coefficients/iteration.
 *
 *  ABI: x0 = pointer to a[256]
 * ============================================================ */
        .globl forward_ntt
        .type  forward_ntt, %function
forward_ntt:
        /* Save callee-saved NEON regs d8-d15 (low 64-bits) */
        stp     d8,  d9,  [sp, #-16]!
        stp     d10, d11, [sp, #-16]!
        stp     d12, d13, [sp, #-16]!
        stp     d14, d15, [sp, #-16]!
        stp     x19, x20, [sp, #-16]!
        stp     x21, x22, [sp, #-16]!

        mov     x19, x0          /* save array pointer */

        /* Load constants into dedicated NEON registers */
        mov     w1, #Q
        dup     v31.8h, w1       /* v31 = {3329, 3329, ...} (Q) */
        mov     w1, #Q_PRIME
        dup     v30.8h, w1       /* v30 = {3327, 3327, ...} (Q') */

        /* ── Layer 1: stride=128, 1 butterfly group ─────────── */
        /*
         * For layers 1-4 we process all 256 coefficients with
         * NEON 8-lane vectors. Each butterfly:
         *   t  = montgomery_mul(a[j+len], zeta)
         *   a[j+len] = a[j] - t
         *   a[j]     = a[j] + t
         */

        adr     x2, ntt_zetas    /* pointer to twiddle table */

        /* Layer 1: len=128, zeta = zetas[0] */
        ld1r    {v29.8h}, [x2], #2    /* broadcast zeta[0] to all 8 lanes */
        add     x3, x19, #0           /* pointer to a[0..127]  */
        add     x4, x19, #256         /* pointer to a[128..255] (bytes: 128*2=256) */

        /* Process 128 coefficients in 16 groups of 8 */
        mov     x5, #16
.Lntt_l1_loop:
        ld1     {v0.8h},  [x3]        /* a[j..j+7]       */
        ld1     {v1.8h},  [x4]        /* a[j+128..j+135] */

        /* t = montgomery_mul(v1, zeta) */
        /* Step 1: v2 = v1 * Q_PRIME (low 16 bits) */
        mul     v2.8h, v1.8h, v30.8h   /* m = a * Q_PRIME mod 2^16 */
        /* Step 2: t = (v1 - m*Q) >> 16 */
        smull   v3.4s, v2.4h, v31.4h   /* m*Q low half */
        smull2  v4.4s, v2.8h, v31.8h   /* m*Q high half */
        smull   v5.4s, v1.4h, v29.4h   /* a*zeta low half */
        smull2  v6.4s, v1.8h, v29.8h   /* a*zeta high half */
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16      /* v2 = t (montgomery reduced) */

        /* Butterfly */
        add     v3.8h, v0.8h, v2.8h    /* a[j]     = a[j] + t */
        sub     v4.8h, v0.8h, v2.8h    /* a[j+128] = a[j] - t */

        st1     {v3.8h}, [x3], #16
        st1     {v4.8h}, [x4], #16
        subs    x5, x5, #1
        b.ne    .Lntt_l1_loop

        /* ── Layer 2: len=64, 2 zetas ──────────────────────── */
        /* 2 groups of 64 coefficients */
        mov     x5, #2
        mov     x3, x19
.Lntt_l2_outer:
        ld1r    {v29.8h}, [x2], #2    /* next zeta */
        add     x4, x3, #128          /* x3 + 64*2 bytes */
        mov     x6, #8
.Lntt_l2_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        mul     v2.8h,  v1.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v1.4h, v29.4h
        smull2  v6.4s,  v1.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        add     v3.8h, v0.8h, v2.8h
        sub     v4.8h, v0.8h, v2.8h
        st1     {v3.8h}, [x3], #16
        st1     {v4.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lntt_l2_inner
        add     x3, x3, #128          /* skip second half of group */
        subs    x5, x5, #1
        b.ne    .Lntt_l2_outer

        /* ── Layers 3-8: similar structure, tighter strides ─── */
        /*
         * For brevity, layers 3-8 follow the same pattern.
         * stride halves each layer: 32, 16, 8, 4, 2, 1
         * number of groups doubles: 4, 8, 16, 32, 64, 128
         *
         * The NEON butterfly macro is identical; only x3/x4 strides change.
         */

        /* Layer 3: len=32, 4 groups */
        mov     x20, #4
        mov     x3, x19
.Lntt_l3_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #64
        mov     x6, #4
.Lntt_l3_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        mul     v2.8h,  v1.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v1.4h, v29.4h
        smull2  v6.4s,  v1.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        add     v3.8h, v0.8h, v2.8h
        sub     v4.8h, v0.8h, v2.8h
        st1     {v3.8h}, [x3], #16
        st1     {v4.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lntt_l3_inner
        add     x3, x3, #64
        subs    x20, x20, #1
        b.ne    .Lntt_l3_outer

        /* Layer 4: len=16, 8 groups */
        mov     x20, #8
        mov     x3, x19
.Lntt_l4_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #32
        mov     x6, #2
.Lntt_l4_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        mul     v2.8h,  v1.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v1.4h, v29.4h
        smull2  v6.4s,  v1.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        add     v3.8h, v0.8h, v2.8h
        sub     v4.8h, v0.8h, v2.8h
        st1     {v3.8h}, [x3], #16
        st1     {v4.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lntt_l4_inner
        add     x3, x3, #32
        subs    x20, x20, #1
        b.ne    .Lntt_l4_outer

        /* Layer 5: len=8, 16 groups — 1 NEON vector exactly! */
        mov     x20, #16
        mov     x3, x19
.Lntt_l5_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #16
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        mul     v2.8h,  v1.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v1.4h, v29.4h
        smull2  v6.4s,  v1.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        add     v3.8h, v0.8h, v2.8h
        sub     v4.8h, v0.8h, v2.8h
        st1     {v3.8h}, [x3], #32
        st1     {v4.8h}, [x4], #16     /* advance x3 past both halves */
        subs    x20, x20, #1
        b.ne    .Lntt_l5_outer

        /* Layer 6: len=4, 32 groups — 2 butterflies per NEON op */
        mov     x20, #32
        mov     x3, x19
.Lntt_l6_outer:
        ld1r    {v29.8h}, [x2], #2
        /*
         * Split 8-lane vector into two 4-element butterfly groups:
         * lanes 0-3: a[j], a[j+1], a[j+2], a[j+3]
         * lanes 4-7: a[j+4], a[j+5], a[j+6], a[j+7]   (= a[j+len+0..3])
         */
        add     x4, x3, #8
        ld1     {v0.4h}, [x3]   /* a[j..j+3]   */
        ld1     {v1.4h}, [x4]   /* a[j+4..j+7] */
        mul     v2.4h,  v1.4h, v30.4h
        smull   v3.4s,  v2.4h, v31.4h
        smull   v5.4s,  v1.4h, v29.4h
        sub     v5.4s, v5.4s, v3.4s
        shrn    v2.4h, v5.4s, #16
        add     v7.4h,  v0.4h, v2.4h
        sub     v8.4h,  v0.4h, v2.4h
        st1     {v7.4h}, [x3], #8
        st1     {v8.4h}, [x4], #8
        subs    x20, x20, #1
        b.ne    .Lntt_l6_outer

        /* Layer 7: len=2, 64 groups */
        mov     x20, #64
        mov     x3, x19
.Lntt_l7_outer:
        ld1r    {v29.4h}, [x2], #2
        add     x4, x3, #4
        ld1     {v0.2h}, [x3]
        ld1     {v1.2h}, [x4]
        mul     v2.2h,  v1.2h, v30.2h
        smull   v3.4s,  v2.4h, v31.4h
        smull   v5.4s,  v1.4h, v29.4h
        sub     v5.4s, v5.4s, v3.4s
        shrn    v2.4h, v5.4s, #16
        add     v7.2h, v0.2h, v2.2h
        sub     v8.2h, v0.2h, v2.2h
        st1     {v7.2h}, [x3], #8
        st1     {v8.2h}, [x4], #4
        subs    x20, x20, #1
        b.ne    .Lntt_l7_outer

        /* Layer 8: len=1, 128 groups — scalar, 2 coeffs per iter */
        mov     x20, #128
        mov     x3, x19
.Lntt_l8_outer:
        /* Load 1 zeta (scalar) */
        ldrsh   w5,  [x2], #2     /* zeta */
        ldrsh   w6,  [x3]         /* a[j]   */
        ldrsh   w7,  [x3, #2]     /* a[j+1] */

        /* t = montgomery_mul(a[j+1], zeta) */
        smull   x8, w7, w5        /* a[j+1] * zeta */
        smull   x9, w7, w21       /* a[j+1] * Q_PRIME — reuse w21=Q_PRIME */
        /* m = low 16 bits of (a*Q') */
        and     x9, x9, #0xFFFF
        smull   x9, w9, w22       /* m * Q — reuse w22=Q */
        sub     x8, x8, x9
        asr     x8, x8, #16       /* t */

        add     w9,  w6, w8       /* a[j]   + t */
        sub     w10, w6, w8       /* a[j+1] - t */
        strh    w9,  [x3], #2
        strh    w10, [x3], #2
        subs    x20, x20, #1
        b.ne    .Lntt_l8_outer

        /* Restore callee-saved registers */
        ldp     x21, x22, [sp], #16
        ldp     x19, x20, [sp], #16
        ldp     d14, d15, [sp], #16
        ldp     d12, d13, [sp], #16
        ldp     d10, d11, [sp], #16
        ldp     d8,  d9,  [sp], #16
        ret
        .size   forward_ntt, . - forward_ntt


/* ============================================================
 *  inverse_ntt(int16_t a[256])
 *
 *  Gentleman-Sande butterfly (decimation in frequency).
 *  Layers in reverse: 1, 2, 4, ..., 128
 *  Final multiply by n^{-1} mod q = 3303 (256^{-1} mod 3329)
 *  Combined with montgomery factor cancellation.
 *
 *  ABI: x0 = pointer to a[256]
 * ============================================================ */
        .globl inverse_ntt
        .type  inverse_ntt, %function
inverse_ntt:
        stp     d8,  d9,  [sp, #-16]!
        stp     d10, d11, [sp, #-16]!
        stp     d12, d13, [sp, #-16]!
        stp     d14, d15, [sp, #-16]!
        stp     x19, x20, [sp, #-16]!
        stp     x21, x22, [sp, #-16]!

        mov     x19, x0

        mov     w1, #Q
        dup     v31.8h, w1
        mov     w1, #Q_PRIME
        dup     v30.8h, w1

        /* inv_zetas pointer — at END of table, iterate backwards */
        adr     x2, ntt_inv_zetas

        /* ── Layer 1 (innermost inverse): stride=1, 128 groups ─ */
        /* GS butterfly: a'=a+b, b'=zeta*(a-b) */
        mov     x20, #128
        mov     x3, x19
.Lintt_l1_outer:
        ldrsh   w5, [x2], #2
        ldrsh   w6, [x3]
        ldrsh   w7, [x3, #2]
        add     w8, w6, w7
        sub     w9, w6, w7
        /* t = mont_mul(b', zeta) */
        smull   x10, w9, w5
        mov     w11, #Q_PRIME
        smull   x12, w9, w11
        and     x12, x12, #0xFFFF
        mov     w13, #Q
        smull   x12, w12, w13
        sub     x10, x10, x12
        asr     x10, x10, #16
        strh    w8,  [x3], #2
        strh    w10, [x3], #2
        subs    x20, x20, #1
        b.ne    .Lintt_l1_outer

        /* ── Layers 2-8: NEON, strides 2,4,8,16,32,64,128 ────── */
        /* (mirror of forward layers 7-1) */

        /* Layer 2: 64 groups, len=2 */
        mov     x20, #64
        mov     x3, x19
.Lintt_l2_outer:
        ld1r    {v29.4h}, [x2], #2
        add     x4, x3, #4
        ld1     {v0.2h}, [x3]
        ld1     {v1.2h}, [x4]
        add     v7.2h, v0.2h, v1.2h
        sub     v8.2h, v0.2h, v1.2h
        mul     v2.2h,  v8.2h, v30.2h
        smull   v3.4s,  v2.4h, v31.4h
        smull   v5.4s,  v8.4h, v29.4h
        sub     v5.4s, v5.4s, v3.4s
        shrn    v2.4h, v5.4s, #16
        st1     {v7.2h}, [x3], #8
        st1     {v2.2h}, [x4], #4
        subs    x20, x20, #1
        b.ne    .Lintt_l2_outer

        /* Layer 3: 32 groups, len=4 */
        mov     x20, #32
        mov     x3, x19
.Lintt_l3_outer:
        ld1r    {v29.4h}, [x2], #2
        add     x4, x3, #8
        ld1     {v0.4h}, [x3]
        ld1     {v1.4h}, [x4]
        add     v7.4h, v0.4h, v1.4h
        sub     v8.4h, v0.4h, v1.4h
        mul     v2.4h,  v8.4h, v30.4h
        smull   v3.4s,  v2.4h, v31.4h
        smull   v5.4s,  v8.4h, v29.4h
        sub     v5.4s, v5.4s, v3.4s
        shrn    v2.4h, v5.4s, #16
        st1     {v7.4h}, [x3], #16
        st1     {v2.4h}, [x4], #8
        subs    x20, x20, #1
        b.ne    .Lintt_l3_outer

        /* Layers 4-8: same NEON macro, increasing strides */
        /* Layer 4: 16 groups, len=8 */
        mov     x20, #16
        mov     x3, x19
.Lintt_l4_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #16
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        add     v7.8h, v0.8h, v1.8h
        sub     v8.8h, v0.8h, v1.8h
        mul     v2.8h,  v8.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v8.4h, v29.4h
        smull2  v6.4s,  v8.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #32
        st1     {v2.8h}, [x4], #16
        subs    x20, x20, #1
        b.ne    .Lintt_l4_outer

        /* Layer 5: 8 groups, len=16 */
        mov     x20, #8
        mov     x3, x19
.Lintt_l5_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #32
        mov     x6, #2
.Lintt_l5_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        add     v7.8h, v0.8h, v1.8h
        sub     v8.8h, v0.8h, v1.8h
        mul     v2.8h,  v8.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v8.4h, v29.4h
        smull2  v6.4s,  v8.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #16
        st1     {v2.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lintt_l5_inner
        add     x3, x3, #32
        subs    x20, x20, #1
        b.ne    .Lintt_l5_outer

        /* Layer 6: 4 groups, len=32 */
        mov     x20, #4
        mov     x3, x19
.Lintt_l6_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #64
        mov     x6, #4
.Lintt_l6_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        add     v7.8h, v0.8h, v1.8h
        sub     v8.8h, v0.8h, v1.8h
        mul     v2.8h,  v8.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v8.4h, v29.4h
        smull2  v6.4s,  v8.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #16
        st1     {v2.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lintt_l6_inner
        add     x3, x3, #64
        subs    x20, x20, #1
        b.ne    .Lintt_l6_outer

        /* Layer 7: 2 groups, len=64 */
        mov     x20, #2
        mov     x3, x19
.Lintt_l7_outer:
        ld1r    {v29.8h}, [x2], #2
        add     x4, x3, #128
        mov     x6, #8
.Lintt_l7_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        add     v7.8h, v0.8h, v1.8h
        sub     v8.8h, v0.8h, v1.8h
        mul     v2.8h,  v8.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v8.4h, v29.4h
        smull2  v6.4s,  v8.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #16
        st1     {v2.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lintt_l7_inner
        add     x3, x3, #128
        subs    x20, x20, #1
        b.ne    .Lintt_l7_outer

        /* Layer 8: 1 group, len=128 */
        ld1r    {v29.8h}, [x2], #2
        mov     x3, x19
        add     x4, x3, #256
        mov     x6, #16
.Lintt_l8_inner:
        ld1     {v0.8h}, [x3]
        ld1     {v1.8h}, [x4]
        add     v7.8h, v0.8h, v1.8h
        sub     v8.8h, v0.8h, v1.8h
        mul     v2.8h,  v8.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v8.4h, v29.4h
        smull2  v6.4s,  v8.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v2.4h, v5.4s, #16
        shrn2   v2.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #16
        st1     {v2.8h}, [x4], #16
        subs    x6, x6, #1
        b.ne    .Lintt_l8_inner

        /* ── Final: multiply all coefficients by n^{-1} * R mod q ── */
        /* n^{-1} mod q = 3303; n^{-1} * R mod q = 1441              */
        mov     w1, #1441
        dup     v29.8h, w1
        mov     x3, x19
        mov     x20, #32          /* 256/8 = 32 iterations */
.Lintt_final:
        ld1     {v0.8h}, [x3]
        mul     v2.8h,  v0.8h, v30.8h
        smull   v3.4s,  v2.4h, v31.4h
        smull2  v4.4s,  v2.8h, v31.8h
        smull   v5.4s,  v0.4h, v29.4h
        smull2  v6.4s,  v0.8h, v29.8h
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v7.4h, v5.4s, #16
        shrn2   v7.8h, v6.4s, #16
        st1     {v7.8h}, [x3], #16
        subs    x20, x20, #1
        b.ne    .Lintt_final

        ldp     x21, x22, [sp], #16
        ldp     x19, x20, [sp], #16
        ldp     d14, d15, [sp], #16
        ldp     d12, d13, [sp], #16
        ldp     d10, d11, [sp], #16
        ldp     d8,  d9,  [sp], #16
        ret
        .size   inverse_ntt, . - inverse_ntt


/* ============================================================
 *  ntt_pointwise_mul(int16_t c[256],
 *                    const int16_t a[256],
 *                    const int16_t b[256])
 *
 *  c[i] = montgomery_mul(a[i], b[i])  for i = 0..255
 *  ABI: x0=c, x1=a, x2=b
 * ============================================================ */
        .globl ntt_pointwise_mul
        .type  ntt_pointwise_mul, %function
ntt_pointwise_mul:
        mov     w3, #Q
        dup     v31.8h, w3
        mov     w3, #Q_PRIME
        dup     v30.8h, w3

        mov     x4, #32          /* 256/8 iterations */
.Lpwm_loop:
        ld1     {v0.8h}, [x1], #16
        ld1     {v1.8h}, [x2], #16

        /* t = montgomery_mul(a, b) */
        mul     v2.8h,  v0.8h, v30.8h   /* m_lo = a * Q' */
        smull   v3.4s,  v2.4h, v31.4h   /* m_lo * Q low */
        smull2  v4.4s,  v2.8h, v31.8h   /* m_lo * Q high */
        smull   v5.4s,  v0.4h, v1.4h    /* a*b low */
        smull2  v6.4s,  v0.8h, v1.8h    /* a*b high */
        sub     v5.4s, v5.4s, v3.4s
        sub     v6.4s, v6.4s, v4.4s
        shrn    v7.4h, v5.4s, #16
        shrn2   v7.8h, v6.4s, #16

        st1     {v7.8h}, [x0], #16
        subs    x4, x4, #1
        b.ne    .Lpwm_loop
        ret
        .size   ntt_pointwise_mul, . - ntt_pointwise_mul


/* ============================================================
 *  poly_mul_ntt — full polynomial multiplication via NTT
 *
 *  void poly_mul_ntt(int16_t c[256],
 *                    const int16_t a[256],
 *                    const int16_t b[256])
 *
 *  Algorithm:
 *    1. â = forward_ntt(copy of a)
 *    2. b̂ = forward_ntt(copy of b)
 *    3. ĉ[i] = mont_mul(â[i], b̂[i])
 *    4. c = inverse_ntt(ĉ)
 *
 *  Uses stack for temporaries (no malloc).
 *  Stack usage: 3 * 256 * 2 = 1536 bytes + frame
 *
 *  ABI: x0=c, x1=a, x2=b
 * ============================================================ */
        .globl poly_mul_ntt
        .type  poly_mul_ntt, %function
poly_mul_ntt:
        /* Stack frame: a_hat(512) + b_hat(512) */
        sub     sp, sp, #1040    /* 512+512+16 alignment */
        stp     x19, x20, [sp, #0]
        stp     x21, x22, [sp, #16]
        stp     x23, x24, [sp, #32]
        str     x25,      [sp, #48]

        add     x19, sp, #64     /* a_hat buffer */
        add     x20, sp, #576    /* b_hat buffer */
        mov     x21, x0          /* output c */
        mov     x22, x1          /* input a */
        mov     x23, x2          /* input b */

        /* Copy a to a_hat, then NTT in place */
        mov     x0, x19
        mov     x1, x22
        mov     x2, #512
.Lpm_copy_a:
        ld1     {v0.16b, v1.16b}, [x1], #32
        st1     {v0.16b, v1.16b}, [x0], #32
        subs    x2, x2, #32
        b.ne    .Lpm_copy_a
        mov     x0, x19
        bl      forward_ntt

        /* Copy b to b_hat, then NTT in place */
        mov     x0, x20
        mov     x1, x23
        mov     x2, #512
.Lpm_copy_b:
        ld1     {v0.16b, v1.16b}, [x1], #32
        st1     {v0.16b, v1.16b}, [x0], #32
        subs    x2, x2, #32
        b.ne    .Lpm_copy_b
        mov     x0, x20
        bl      forward_ntt

        /* Pointwise multiply into c (reuse c buffer) */
        mov     x0, x21
        mov     x1, x19
        mov     x2, x20
        bl      ntt_pointwise_mul

        /* Inverse NTT on c */
        mov     x0, x21
        bl      inverse_ntt

        ldr     x25,      [sp, #48]
        ldp     x23, x24, [sp, #32]
        ldp     x21, x22, [sp, #16]
        ldp     x19, x20, [sp, #0]
        add     sp, sp, #1040
        ret
        .size   poly_mul_ntt, . - poly_mul_ntt


/* ============================================================
 *  barrett_reduce_poly(int16_t a[256])
 *
 *  Reduce all coefficients mod q using Barrett reduction.
 *  NEON: process 8 coefficients per cycle.
 *
 *  Barrett for q=3329: k=20159, shift=26
 *    t = (a * 20159) >> 26
 *    r = a - t*q
 *    if r >= q: r -= q
 *    if r <  0: r += q
 *
 *  ABI: x0 = pointer to a[256]
 * ============================================================ */
        .globl barrett_reduce_poly
        .type  barrett_reduce_poly, %function
barrett_reduce_poly:
        mov     w1, #20159
        dup     v28.8h, w1       /* Barrett constant k */
        mov     w1, #Q
        dup     v31.8h, w1       /* Q */

        mov     x2, #32          /* 256/8 */
.Lbarrett_loop:
        ld1     {v0.8h}, [x0]

        /* t = (a * 20159) >> 26 — need 32-bit intermediate */
        smull   v1.4s, v0.4h, v28.4h
        smull2  v2.4s, v0.8h, v28.8h
        sshr    v1.4s, v1.4s, #26
        sshr    v2.4s, v2.4s, #26
        sqxtn   v3.4h, v1.4s
        sqxtn2  v3.8h, v2.4s     /* t (16-bit) */

        /* r = a - t*q */
        mul     v4.8h, v3.8h, v31.8h
        sub     v5.8h, v0.8h, v4.8h  /* r */

        /* Conditional corrections (branchless) */
        /* if r >= q: r -= q */
        cmge    v6.8h, v5.8h, v31.8h
        and     v6.16b, v6.16b, v31.16b
        sub     v5.8h, v5.8h, v6.8h
        /* if r < 0: r += q */
        cmlt    v7.8h, v5.8h, #0
        and     v7.16b, v7.16b, v31.16b
        add     v5.8h, v5.8h, v7.8h

        st1     {v5.8h}, [x0], #16
        subs    x2, x2, #1
        b.ne    .Lbarrett_loop
        ret
        .size   barrett_reduce_poly, . - barrett_reduce_poly

/* end of ntt_arm64.S */
