/*
 * ============================================================
 *  randombytes.c  —  OS Entropy Source  (v2)
 *
 *  Platform detection:
 *    - Windows  : BCryptGenRandom (no /dev/urandom)
 *    - Linux    : getrandom(2) syscall  [preferred, no fd leak]
 *    - macOS    : getentropy(2)
 *    - Fallback : /dev/urandom (POSIX)
 *
 *  Security properties:
 *    - Constant-time: timing does not depend on output bytes
 *    - No global state
 *    - Zero dynamic allocation
 * ============================================================
 */

#include "randombytes.h"
#include <string.h>

/* ── Platform detection ───────────────────────────────────── */
#if defined(_WIN32) || defined(_WIN64)
  #define PLATFORM_WINDOWS
#elif defined(__linux__)
  #define PLATFORM_LINUX
#elif defined(__APPLE__)
  #define PLATFORM_APPLE
#else
  #define PLATFORM_POSIX
#endif

/* ── Windows ──────────────────────────────────────────────── */
#ifdef PLATFORM_WINDOWS
#include <windows.h>
#include <bcrypt.h>
#pragma comment(lib, "bcrypt.lib")

int randombytes(uint8_t *out, size_t outlen) {
    if (outlen == 0) return 0;
    NTSTATUS status = BCryptGenRandom(
        NULL, (PUCHAR)out, (ULONG)outlen,
        BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    return BCRYPT_SUCCESS(status) ? 0 : -1;
}

/* ── Linux: getrandom(2) ──────────────────────────────────── */
#elif defined(PLATFORM_LINUX)
#include <sys/random.h>
#include <errno.h>

int randombytes(uint8_t *out, size_t outlen) {
    size_t done = 0;
    while (done < outlen) {
        ssize_t r = getrandom(out + done, outlen - done, 0);
        if (r < 0) {
            if (errno == EINTR) continue; /* retry on signal */
            return -1;
        }
        done += (size_t)r;
    }
    return 0;
}

/* ── macOS: getentropy(2) ─────────────────────────────────── */
#elif defined(PLATFORM_APPLE)
#include <sys/random.h>

int randombytes(uint8_t *out, size_t outlen) {
    /* getentropy is limited to 256 bytes per call */
    size_t done = 0;
    while (done < outlen) {
        size_t chunk = outlen - done;
        if (chunk > 256) chunk = 256;
        if (getentropy(out + done, chunk) != 0) return -1;
        done += chunk;
    }
    return 0;
}

/* ── POSIX fallback: /dev/urandom ─────────────────────────── */
#else
#include <stdio.h>

int randombytes(uint8_t *out, size_t outlen) {
    if (outlen == 0) return 0;
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) return -1;
    size_t r = fread(out, 1, outlen, f);
    fclose(f);
    /* Wipe any partial read to avoid silent short reads */
    if (r != outlen) {
        memset(out, 0, outlen); /* don't use partial entropy */
        return -1;
    }
    return 0;
}
#endif
