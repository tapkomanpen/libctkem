import ctypes
import os

class CTKEM:
    PK_SIZE = 448
    SK_SIZE = 416
    CT_SIZE = 832
    SS_SIZE = 32

    def __init__(self, lib_path="./libctkem.so"):
        self.lib = ctypes.CDLL(os.path.abspath(lib_path))

        self.lib.kem_keypair.argtypes = [
            ctypes.POINTER(ctypes.c_uint8 * self.PK_SIZE),
            ctypes.POINTER(ctypes.c_uint8 * self.SK_SIZE)
        ]
        self.lib.kem_keypair.restype = ctypes.c_int

        self.lib.kem_encaps.argtypes = [
            ctypes.POINTER(ctypes.c_uint8 * self.CT_SIZE),
            ctypes.POINTER(ctypes.c_uint8 * self.SS_SIZE),
            ctypes.POINTER(ctypes.c_uint8 * self.PK_SIZE)
        ]
        self.lib.kem_encaps.restype = ctypes.c_int

        self.lib.kem_decaps.argtypes = [
            ctypes.POINTER(ctypes.c_uint8 * self.SS_SIZE),
            ctypes.POINTER(ctypes.c_uint8 * self.CT_SIZE),
            ctypes.POINTER(ctypes.c_uint8 * self.SK_SIZE)
        ]
        self.lib.kem_decaps.restype = None

    def generate_keypair(self):
        pk = (ctypes.c_uint8 * self.PK_SIZE)()
        sk = (ctypes.c_uint8 * self.SK_SIZE)()
        status = self.lib.kem_keypair(ctypes.byref(pk), ctypes.byref(sk))
        if status != 0:
            raise RuntimeError("Keypair generation failed")
        return bytes(pk), bytes(sk)

    def encapsulate(self, public_key):
        pk_buf = (ctypes.c_uint8 * self.PK_SIZE).from_buffer_copy(public_key)
        ct = (ctypes.c_uint8 * self.CT_SIZE)()
        ss = (ctypes.c_uint8 * self.SS_SIZE)()
        status = self.lib.kem_encaps(ctypes.byref(ct), ctypes.byref(ss), ctypes.byref(pk_buf))
        if status != 0:
            raise RuntimeError("Encapsulation failed")
        return bytes(ct), bytes(ss)

    def decapsulate(self, ciphertext, secret_key):
        ct_buf = (ctypes.c_uint8 * self.CT_SIZE).from_buffer_copy(ciphertext)
        sk_buf = (ctypes.c_uint8 * self.SK_SIZE).from_buffer_copy(secret_key)
        ss = (ctypes.c_uint8 * self.SS_SIZE)()
        self.lib.kem_decaps(ctypes.byref(ss), ctypes.byref(ct_buf), ctypes.byref(sk_buf))
        return bytes(ss)

if __name__ == "__main__":
    kem = CTKEM()
    
    pk, sk = kem.generate_keypair()
    print(f"Public Key:  {pk[:16].hex()}...")
    
    ct, ss_e = kem.encapsulate(pk)
    print(f"Ciphertext:  {ct[:16].hex()}...")
    print(f"Shared Alice: {ss_e.hex()}")
    
    ss_d = kem.decapsulate(ct, sk)
    print(f"Shared Bob:   {ss_d.hex()}")
    
    if ss_e == ss_d:
        print("Status: Success")
    else:
        print("Status: Failure")
