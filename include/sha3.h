/*
 * ============================================================
 *  sha3.h  —  SHA-3 / SHAKE256 (Keccak-f[1600])
 *
 *  Implements:
 *    - SHA3-256  (32-byte hash)
 *    - SHAKE256  (eXtendable Output Function, XOF)
 *
 *  Properties:
 *    - Zero dynamic allocation (stack only)
 *    - Constant-time (no secret-dependent branches)
 *    - Self-contained — no external dependencies
 *    - NIST FIPS 202 compliant
 *
 *  Used by rlwe_kem.c for:
 *    - poly_uniform()  →  SHAKE256(seed, domain)
 *    - kem_keypair()   →  SHA3-256(noise_seed) as PRF
 *    - FO Transform    →  SHA3-256(m ‖ pk) → shared secret
 * ============================================================
 */

#ifndef SHA3_H
#define SHA3_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* ── Keccak constants ──────────────────────────────────────── */
#define KECCAK_ROUNDS 24

static const uint64_t KECCAK_RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL,
};

static const int KECCAK_RHO[24] = {
     1,  3,  6, 10, 15, 21, 28, 36,
    45, 55,  2, 14, 27, 41, 56,  8,
    25, 43, 62, 18, 39, 61, 20, 44,
};

static const int KECCAK_PI[24] = {
    10,  7, 11, 17, 18,  3,  5, 16,
     8, 21, 24,  4, 15, 23, 19, 13,
    12,  2, 20, 14, 22,  9,  6,  1,
};

/* ── Keccak-f[1600] permutation ───────────────────────────── */
static inline uint64_t rotl64(uint64_t x, int n) {
    return (x << n) | (x >> (64 - n));
}

static void keccak_f1600(uint64_t st[25]) {
    uint64_t bc[5], tmp;

    for (int r = 0; r < KECCAK_ROUNDS; r++) {
        /* θ */
        for (int i = 0; i < 5; i++)
            bc[i] = st[i] ^ st[i+5] ^ st[i+10] ^ st[i+15] ^ st[i+20];
        for (int i = 0; i < 5; i++) {
            tmp = bc[(i+4)%5] ^ rotl64(bc[(i+1)%5], 1);
            for (int j = 0; j < 25; j += 5) st[j+i] ^= tmp;
        }
        /* ρ and π */
        tmp = st[1];
        for (int i = 0; i < 24; i++) {
            int j = KECCAK_PI[i];
            uint64_t t = st[j];
            st[j] = rotl64(tmp, KECCAK_RHO[i]);
            tmp = t;
        }
        /* χ */
        for (int j = 0; j < 25; j += 5) {
            uint64_t t[5];
            for (int i = 0; i < 5; i++) t[i] = st[j+i];
            for (int i = 0; i < 5; i++)
                st[j+i] ^= (~t[(i+1)%5]) & t[(i+2)%5];
        }
        /* ι */
        st[0] ^= KECCAK_RC[r];
    }
}

/* ============================================================
 *  SHAKE256 XOF  (rate = 136 bytes = 1088 bits)
 * ============================================================ */
#define SHAKE256_RATE 136

typedef struct {
    uint64_t st[25];
    uint8_t  buf[SHAKE256_RATE];
    int      absorb_pos; /* current byte position within rate block (absorb) */
    int      buf_pos;    /* squeeze position */
    int      absorbed;   /* 1 after pad+permute */
} shake256_ctx;

static void shake256_init(shake256_ctx *ctx) {
    memset(ctx->st,  0, sizeof(ctx->st));
    memset(ctx->buf, 0, sizeof(ctx->buf));
    ctx->absorb_pos = 0;
    ctx->buf_pos    = SHAKE256_RATE; /* force refill on first squeeze */
    ctx->absorbed   = 0;
}

static void shake256_absorb(shake256_ctx *ctx,
                             const uint8_t *in, size_t inlen)
{
    uint8_t *s   = (uint8_t *)ctx->st;
    size_t   pos = 0;
    while (pos < inlen) {
        s[ctx->absorb_pos] ^= in[pos];
        ctx->absorb_pos++;
        pos++;
        if (ctx->absorb_pos == SHAKE256_RATE) {
            keccak_f1600(ctx->st);
            ctx->absorb_pos = 0;
        }
    }
}

static void shake256_finalize(shake256_ctx *ctx) {
    /* SHAKE256 domain: 0x1F at current position, 0x80 at rate-1 */
    uint8_t *s = (uint8_t *)ctx->st;
    s[ctx->absorb_pos] ^= 0x1F;
    s[SHAKE256_RATE - 1] ^= 0x80;
    keccak_f1600(ctx->st);
    /* copy state to buf for squeezing */
    memcpy(ctx->buf, ctx->st, SHAKE256_RATE);
    ctx->buf_pos  = 0;
    ctx->absorbed = 1;
}

static void shake256_squeeze(shake256_ctx *ctx,
                              uint8_t *out, size_t outlen)
{
    size_t pos = 0;
    while (pos < outlen) {
        if (ctx->buf_pos >= SHAKE256_RATE) {
            keccak_f1600(ctx->st);
            memcpy(ctx->buf, ctx->st, SHAKE256_RATE);
            ctx->buf_pos = 0;
        }
        size_t avail = (size_t)(SHAKE256_RATE - ctx->buf_pos);
        size_t need  = outlen - pos;
        size_t take  = need < avail ? need : avail;
        memcpy(out + pos, ctx->buf + ctx->buf_pos, take);
        ctx->buf_pos += (int)take;
        pos          += take;
    }
}

/* ── One-shot SHAKE256 ─────────────────────────────────────── */
static void shake256(uint8_t *out, size_t outlen,
                     const uint8_t *in, size_t inlen)
{
    shake256_ctx ctx;
    shake256_init(&ctx);
    shake256_absorb(&ctx, in, inlen);
    shake256_finalize(&ctx);
    shake256_squeeze(&ctx, out, outlen);
}

/* ============================================================
 *  SHA3-256  (rate = 136 bytes, output = 32 bytes)
 *  Domain suffix = 0x06 (differs from SHAKE)
 * ============================================================ */
static void sha3_256(uint8_t out[32],
                     const uint8_t *in, size_t inlen)
{
    uint64_t st[25];
    memset(st, 0, sizeof(st));
    uint8_t *s   = (uint8_t *)st;
    int      apos = 0;

    for (size_t i = 0; i < inlen; i++) {
        s[apos++] ^= in[i];
        if (apos == SHAKE256_RATE) { keccak_f1600(st); apos = 0; }
    }
    /* SHA3-256 domain = 0x06 */
    s[apos]              ^= 0x06;
    s[SHAKE256_RATE - 1] ^= 0x80;
    keccak_f1600(st);
    memcpy(out, st, 32);
}

/* ── Convenience: SHAKE256 with domain byte ────────────────── */
static void shake256_domain(uint8_t *out, size_t outlen,
                             const uint8_t *seed, size_t seedlen,
                             uint8_t domain)
{
    shake256_ctx ctx;
    shake256_init(&ctx);
    shake256_absorb(&ctx, seed, seedlen);
    shake256_absorb(&ctx, &domain, 1);
    shake256_finalize(&ctx);
    shake256_squeeze(&ctx, out, outlen);
}

#endif /* SHA3_H */
