/*
 * ============================================================
 *  randombytes.h  —  Entropy Source Interface  (v2)
 *
 *  Platform-agnostic CSPRNG interface.
 *  Implementations: Linux/macOS (/dev/urandom),
 *                   Windows (BCryptGenRandom),
 *                   Embedded (HW TRNG register).
 *
 *  Contract:
 *    - Returns 0 on success, -1 on failure.
 *    - Output is cryptographically strong (CSPRNG quality).
 *    - Must NOT be used for deterministic key derivation —
 *      use SHAKE256 for that.
 * ============================================================
 */

#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

#include <stdint.h>
#include <stddef.h>

/**
 * Fill buffer with cryptographically strong random bytes.
 *
 * @param out    [out] Output buffer
 * @param outlen [in]  Number of bytes to generate
 * @return 0 on success, -1 on failure.
 */
int randombytes(uint8_t *out, size_t outlen);

#endif /* RANDOMBYTES_H */
