/*
 * ============================================================
 * randombytes.h  —  Entropy Source Interface
 *
 * This header defines the interface for the system's
 * cryptographically secure random number generator (CSPRNG).
 * The implementation should be provided based on the target 
 * platform (e.g., /dev/urandom, BCryptGenRandom, or HW TRNG).
 * ============================================================
 */

#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

#include <stdint.h>

/**
 * Fill a buffer with cryptographically strong random bytes.
 * * @param out     [out] Pointer to the output buffer
 * @param outlen  [in]  Number of bytes to generate
 * * @return 0 on success, -1 on failure.
 */
int randombytes(uint8_t *out, size_t outlen);

#endif /* RANDOMBYTES_H */
