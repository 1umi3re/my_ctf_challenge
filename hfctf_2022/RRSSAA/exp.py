from gmpy2 import next_prime, iroot
from Crypto.Util.number import getPrime, inverse, GCD, bytes_to_long, long_to_bytes
from sage.all import *

def attack2(N, e, m, t, X, Y):
    PR = PolynomialRing(QQ, 'x,y', 2, order='lex')
    x, y = PR.gens()
    A = -(N-1)**2
    F = x * y**2 + A * x + 1

    G_polys = []
    # G_{k,i_1,i_2}(x,y) = x^{i_1-k}y_{i_2-2k}f(x,y)^{k}e^{m-k} 
    for k in range(m + 1):
        for i_1 in range(k, m+1):
            for i_2 in [2*k, 2*k + 1]:
                G_polys.append(x**(i_1-k) * y**(i_2-2*k) * F**k * e**(m-k))

    H_polys = []
    # y_shift H_{k,i_1,i_2}(x,y) = y^{i_2-2k} f(x,y)^k e^{m-k}
    for k in range(m + 1):
        for i_2 in range(2*k+2, 2*k+t+1):
            H_polys.append(y**(i_2-2*k) * F**k * e**(m-k))

    polys = G_polys + H_polys
    monomials = []
    for poly in polys:
        monomials.append(poly.lm())
    
    dims1 = len(polys)
    dims2 = len(monomials)
    MM = matrix(QQ, dims1, dims2)
    for idx, poly in enumerate(polys):
        for idx_, monomial in enumerate(monomials):
            if monomial in poly.monomials():
                MM[idx, idx_] = poly.monomial_coefficient(monomial) * monomial(X, Y)
    B = MM.LLL()

    found_polynomials = False

    for pol1_idx in range(B.nrows()):
        for pol2_idx in range(pol1_idx + 1, B.nrows()):
            P = PolynomialRing(QQ, 'a,b', 2)
            a, b = P.gens()
            pol1 = pol2 = 0
            for idx_, monomial in enumerate(monomials):
                pol1 += monomial(a,b) * B[pol1_idx, idx_] / monomial(X, Y)
                pol2 += monomial(a,b) * B[pol2_idx, idx_] / monomial(X, Y)

            # resultant
            rr = pol1.resultant(pol2)
            # are these good polynomials?
            if rr.is_zero() or rr.monomials() == [1]:
                continue
            else:
                print(f"found them, using vectors {pol1_idx}, {pol2_idx}")
                found_polynomials = True
                break
        if found_polynomials:
            break

    if not found_polynomials:
        print("no independant vectors could be found. This should very rarely happen...")


    PRq = PolynomialRing(QQ, 'z')
    z = PRq.gen()
    rr = rr(z, z)
    soly = rr.roots()[0][0]

    ppol = pol1(z, soly)
    solx = ppol.roots()[0][0]
    return solx, soly


def seq(r, k, m):
    v = vector(Zmod(m), [r, 2])
    if k >= 2:
        M = Matrix(Zmod(m), [[r, -1], [1, 0]])
        v = (M**(k-1)) * v
    ret = v[0] if k != 0 else v[1]
    return int(ret)


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

if __name__ == '__main__':
    n = 59969098213446598961510550233718258878862148298191323654672950330070587404726715299685997489142290693126366408044603303463518341243526241117556011994804902686998166238333549719269703453450958140262475942580009981324936992976252832887660977703209225426388975233018602730303262439218292062822981478737257836581
    e = 970698965238639683403205181589498135440069660016843488485401994654202837058754446853559143754852628922125327583411039117445415303888796067576548626904070971514824878024057391507617988385537930417136322298476467215300995795105008488692961624917433064070351961856959734368784774555385603000155569897078026670993484466622344106374637350023474339105113172687604783395923403613555236693496567851779400707953027457705617050061193750124237055690801725151098972239120476113241310088089420901051617493693842562637896252448161948655455277146925913049354086353328749354876619287042077221173795354616472050669799421983520421287
    c = 2757297249371055260112176788534868300821961060153993508569437878576838431569949051806118959108641317578931985550844206475198216543139472405873345269094341570473142756599117266569746703013099627523306340748466413993624965897996985230542275127290795414763432332819334757831671028121489964563214463689614865416498886490980692515184662350519034273510244222407505570929178897273048405431658365659592815446583970229985655015539079874797518564867199632672678818617933927005198847206019475149998468493858071672920824599672525667187482558622701227716212254925837398813278836428805193481064316937182435285668656233017810444672

    alpha = ZZ(e).nbits() / ZZ(n).nbits()
    beta = 0.44
    nbits = 1024
    delta = 0.63

    X = 2 ** int(nbits*(alpha+delta-2)+3)
    Y = 2 ** int(nbits*beta+3)

    x, y = map(int, attack2(n, e, 8, 12, X, Y))
    p_minus_q = y
    p_plus_q = iroot(p_minus_q**2 + 4 * n, 2)[0]

    p = (p_minus_q + p_plus_q) // 2
    q = n // p
    assert p * q == n
    phi = (p**2 - 1) * (q**2 - 1)
    d = inverse(e, phi)
    m = decrypt(c, e, p, q)
    print(long_to_bytes(m))