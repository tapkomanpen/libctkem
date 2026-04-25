#include "randombytes.h"
#include <stdio.h>

int randombytes(uint8_t *out, size_t outlen) {
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) return -1;
    size_t read = fread(out, 1, outlen, f);
    fclose(f);
    return (read == outlen) ? 0 : -1;
}
