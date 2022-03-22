from gmpy2 import next_prime
from Crypto.Util.number import getPrime, inverse, GCD, bytes_to_long
from random import getrandbits
from secret import flag
from sage.all import *

def seq(r, k, m):
    v = vector(Zmod(m), [r, 2])
    if k >= 2:
        M = Matrix(Zmod(m), [[r, -1], [1, 0]])
        v = (M**(k-1)) * v
    ret = v[0] if k != 0 else v[1]
    return int(ret)

def encrypt(m, e, n):
    while True:
        r = randint(1, n - 1)
        if r != 2 and r != n - 2 and GCD(r, n) == 1:
            break
    v = seq(r, e, n**2)
    print(r)
    return ((1 + m*n) * v) % n**2

def legendre_symbol(a, p):
    """ Compute the Legendre symbol a|p using
        Euler's criterion. p is a prime, a is
        relatively prime to p (if p divides
        a, then a|p = 0)
        Returns 1 if a has a square root modulo
        p, -1 otherwise.
    """
    ls = pow(2,(p-1)//2,p)
    return -1 if ls == p - 1 else ls

def decrypt(c, e, p, q):
    d_p = {1: int(pow(e, -1, p-1)), -1: int(pow(e, -1, p+1))}
    d_q = {1: int(pow(e, -1, q-1)), -1: int(pow(e, -1, q+1))}

    inv_q = int(pow(p, -1, q))
    inv_p = int(pow(q, -1, p))

    i_p = legendre_symbol(c**2-4, p)
    i_q = legendre_symbol(c**2-4, q)
    r_p = seq(c, d_p[i_p], p)
    r_q = seq(c, d_q[i_q], q)

    r = CRT([r_p, r_q], [p, q])
    v_rp = seq(r, e, p**2)
    t_p = int((c * pow(v_rp, -1, p**2)) % p**2)
    s_p = (t_p - 1) // p

    v_rq = seq(r, e, q**2)
    t_q = int((c * pow(v_rq, -1, q**2)) % q**2)
    s_q = (t_q - 1) // q

    m_p = (s_p * inv_p) % p
    m_q = (s_q * inv_q) % q

    m = CRT([m_p, m_q], [p, q])

    return m


beta = 0.44
nbits = 1024
delta = 0.63

p = getPrime(nbits // 2)
# This is a super dummy generation of prime and it should be modified like below
# q = next_prime(p + (1 << int(nbits * beta))) 
q = next_prime(p + getrandbits(int(nbits*data)))

phi = (p**2 - 1) * (q**2 - 1)

while True:
    d = getPrime(int(nbits * delta))
    if GCD(d, phi) == 1:
        break
e = inverse(d, phi)
alpha = int(e).bit_length() / nbits
print(alpha)

n = p * q
m = bytes_to_long(flag)
while True:
    c = int(encrypt(m, e, n))
    m_ = decrypt(c, e, p, q)
    if m == m_:
        break

print('------------------------------------------')
print(f"n = {n}")
print(f"e = {e}")
print(f"c = {c}")
print(f"d = {d}")
print(f"p = {p}")
print(f"q = {q}")
print('------------------------------------------')