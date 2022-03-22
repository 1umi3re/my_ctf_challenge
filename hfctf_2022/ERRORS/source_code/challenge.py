from ntru import NTRUCipher
from hashlib import sha3_384, sha256
from os import urandom
from Crypto.Util.number import long_to_bytes, bytes_to_long
from Crypto.Util.strxor import strxor as xor
from ast import literal_eval
import signal
from random import choices
import string
from secret import flag

H = lambda x: sha3_384(x).digest()


def proof_of_work():
    alphabet = string.ascii_letters + string.digits
    nonce = "".join(choices(alphabet, k=8))
    print(f'SHA256("{nonce}" + ?) starts with "000000"')
    message = (nonce + input().strip()).encode()
    return sha256(message).digest().hex().startswith("000000")


class Challenge:
    def __init__(self, *params):
        self.cipher = NTRUCipher(*params)

    def encrypt(self, pt: bytes):
        assert len(pt) == 48
        ss = H(pt)
        encoded_msg = self.cipher.encode(bytes_to_long(pt) ^ bytes_to_long(ss))
        ct = self.cipher.encrypt(encoded_msg)
        return ct.list(), ss

    def decrypt(self, ct: list, ss: bytes):
        assert len(ct) == self.cipher.N
        encoded_msg = self.cipher.decrypt(self.cipher.poly_from_list(ct))
        msg = long_to_bytes(self.cipher.decode(encoded_msg)).rjust(48, b'\x00')
        if len(msg) > 48:
            return False, None
        pt = xor(msg, ss)
        ss_ = H(pt)
        if ss_ != ss:
            return False, None
        else:
            return True, pt


# NTRU parameters (N, p, q, d)
PARAMS = (256, 3, 2048, 75)

BANNER = '''I have recently learnt a public key cryptosystem NTRU, which is
more resistant to attacks from quantum computers. Now we can securely share
our secrets with each other. However, since there may be errors in the 
decoding process, I add a hash for verification. So would you like to tell
me some little secrets?
'''


MENU = '''[+] 1. Get my public key
[+] 2. Share your secret
[+] 3. Share my secret 
[+] 4. Exit'''

if __name__ == "__main__":
    signal.alarm(60)
    if not proof_of_work():
        exit(0)
    signal.alarm(0)

    print(BANNER)
    chal = Challenge(*PARAMS)
    
    for _ in range(1024):
        print(MENU)
        c = input("Enter your command: ")
        if c.startswith("1"):
            pk = chal.cipher.h.list()
            print(f"Here is my public key: {pk}")
        elif c.startswith("2"):
            try:
                ct = input("Please input your encrypted secret: ")
                ss = input("Please input the hash of the secret (in hex): ")
                ct = literal_eval(ct)
                ss = bytes.fromhex(ss)
                if len(ss) != 48 or any([type(x) != int for x in ct]):
                    print("Sorry, the format of the ciphertext is not correct.")
                else:
                    verified, msg_ = chal.decrypt(ct, ss)
                    if verified:
                        print("The message is verified. I promise I won't tell anybody your secret.")
                    else:
                        print("There is something wrong in your ciphertext, try again!")
            except:
                print('Malformed ciphertext.')
                pass
        elif c.startswith("3"):
            tmp_flag = flag + urandom(48 - len(flag))
            enc_flag, ss = chal.encrypt(tmp_flag)
            print(f"Here is my little secret: {enc_flag}")
            print(f"Here is the hash of my little secret: {ss.hex()}")
        elif c.startswith("4"):
            print("BYE! HAVE A BEAUTIFUL TIME!")
            exit()
        else:
            print("Unsupported command.")