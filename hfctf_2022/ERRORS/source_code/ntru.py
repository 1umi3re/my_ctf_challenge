from sage.all import *
from random import shuffle, getrandbits

Zx = PolynomialRing(ZZ, 'x')
x = Zx.gen()


def convolution(f, g, R):
    return (f * g) % R


def balancedmod(f, q, R):
    g = list(map(lambda x: ((x + q//2) % q) - q//2, f.list()))
    return Zx(g) % R


def random_poly(n, d1, d2):
    assert d1 + d2 <= n
    result = d1 * [1] + d2 * [-1] + (n - d1 - d2) * [0]
    shuffle(result)
    return Zx(result)


def invert_poly_mod_prime(f, R, p):
    T = Zx.change_ring(Integers(p)).quotient(R)
    return Zx(lift(1 / T(f)))


def invert_poly_mod_powerof2(f, R, q):
    g = invert_poly_mod_prime(f, R, 2)
    e = log(q, 2)
    for i in range(1, e):
        g = ((2 * g - f * g ** 2) % R) % q
    return g


class NTRUCipher:
    def __init__(self, N, p, q, d):
        self.N = N
        self.p = p
        self.q = q
        self.d = d
        self.R = x ** N - 1

        # key generation
        self.g = random_poly(self.N, d, d)
        while True:
            try:
                self.f = random_poly(self.N, d + 1, d)
                self.fp = invert_poly_mod_prime(self.f, self.R, self.p)
                self.fq = invert_poly_mod_powerof2(self.f, self.R, self.q)
                break
            except:
                pass

        self.h = balancedmod(self.p * convolution(self.fq, self.g, self.R), self.q, self.R)

    def encrypt(self, m):
        r = random_poly(self.N, self.d, self.d)
        return balancedmod(convolution(self.h, r, self.R) + m, self.q, self.R)

    def decrypt(self, c):
        a = balancedmod(convolution(c, self.f, self.R), self.q, self.R)
        return balancedmod(convolution(a, self.fp, self.R), self.p, self.R)

    def encode(self, val):
        poly = 0
        for i in range(self.N):
            poly += ((val % self.p) - self.p // 2) * (x ** i)
            val //= self.p
        return poly

    def decode(self, poly):
        result = 0
        ll = poly.list()
        for idx, val in enumerate(ll):
            result += (val + self.p // 2) * (self.p ** idx)
        return result

    def poly_from_list(self, l: list):
        return Zx(l)


if __name__ == '__main__':
    N = 256
    d = 75
    p = 3
    q = 2048

    cipher = NTRUCipher(N, p, q, d)
    msg = getrandbits(384)
    encode_msg = cipher.encode(msg)
    c = cipher.encrypt(encode_msg)
    mm = cipher.decrypt(c)
    decode_msg = cipher.decode(mm)

    assert decode_msg == msg
