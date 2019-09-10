from math import fmod, floor, ceil, exp, log
from gmpy2 import invert
from array import array
import logging
from sympy.ntheory.residue_ntheory import jacobi_symbol


class QuadraticSieve():

    def __init__(self, max_numero):
        self.max_numero = max_numero
        self.primos = array('I')
        self._inicializa_criba()

    def _inicializa_criba(self):
        # XXX: https://stackoverflow.com/questions/521674/initializing-a-list-to-a-known-number-of-elements-in-python
        bandera_primos = array('B', (1,) * (self.max_numero + 1))
        for i in range(2, self.max_numero + 1):
            if bandera_primos[i]:
                self.primos.append(i)
            
            for primo in self.primos:
                compuesto = primo * i
                if compuesto > self.max_numero:
                    break
                bandera_primos[compuesto] = 0
                if not (i % primo):
                    break
                
    def simbolo_legendre(self, a, p):
        s = (a ** ((p - 1) >> 1)) % p
        return s

    # XXX: https://www.johndcook.com/blog/2019/02/12/computing-jacobi-symbols/
    def simbolo_jacobi(self, a, n):
#        logger.debug("a {} n {}".format(a, n))
        assert(n % 2 == 1)
        t = 1
        if a < 0 or a > n:
            a = a % n
            if a < 0:
                a = -a
                if n % 4 == 3:
                    t = -t
            
        while a != 0:
            while a % 2 == 0:
                a /= 2
                r = n % 8
                if r == 3 or r == 5:
                    t = -t
            a, n = n, a
            if a % 4 == n % 4 == 3:
                t = -t
            a %= n
        if n == 1:
            return t
        else:
            return 0
    
    def _encuentra_base_de_factores(self, n, b):
        rc = array("I")
        for p in self.primos:
            if p > b:
                break
            if self.simbolo_legendre(n, p) == 1:
                rc.append(p)
        return rc
    
    def _calcula_limite_b_suave(self, n):
        b = ceil(exp(((log(n) * log(log(n))) / 2) / 2))
        return int(b)
    
    # XXX: https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python
    def _calcula_z_shanks_tonelli(self, p):
        for z in range(2, p):
            if p - 1 == self.simbolo_legendre(z, p):
                break
        return z

    def calcula_conguencia_residuo_cuadratico(self, n, p):
        logger.debug("calculando n {} p {}".format(n, p))
        p_menos_1 = p - 1
        S = 0
        while not (p_menos_1 & 1):
            p_menos_1 >>= 1
            S += 1
        logger.debug("S {} Q {}".format(S, p_menos_1))
        Z = self._calcula_z_shanks_tonelli(p)
        logger.debug("Z {}".format(Z))
        Q = p_menos_1
        c = (Z ** Q) % p
        R = (n ** ((Q + 1) >> 1)) % p
        t = (n ** Q) % p
        M = S
        logger.debug("t {} R {}".format(t, R))
        while (t % p) != 1:
            i = 0
            logger.debug("M {}".format(M))
            for i in range(1, M):
#                 logger.debug("t mod {}".format((t ** (i << 1)) % p))
#                 if ((t ** (i << 1)) % p) == 1:
                if ((t ** (2 ** i)) % p) == 1:
                    break
            b = (c ** (1 << (M - i - 1))) % p
            logger.debug("b {} R {} t {} c {} M {} i {}".format(b, R, t, c, M, i))
            R = (R * b) % p
            # XXX: https://eli.thegreenplace.net/2009/03/07/computing-modular-square-roots-in-python
            t = (t * (b ** 2)) % p
            c = (b ** 2) % p
            M = i
        return R, p - R


def pow_mod(a, n, m):
    r = 1
    pot = a
    aa = 1
    while(n):
        if(n & 1):
#            r=((r%m)*(pot%m))%m
#            print("fmod r {} fmod pot {}".format(fmod(r,m),fmod(pot,m)))
            r = fmod(fmod(r, m) * fmod(pot, m), m)
        n >>= 1
#        pot=((pot%m)*(pot%m))%m
#        print("pot {} de pot {} con a {} a la n {} m {}".format(fmod(pot,m),pot,a,aa,m))
        pot = fmod(fmod(pot, m) * fmod(pot, m), m)
        aa <<= 1
    return r


logging.basicConfig(format='%(asctime)s  %(levelname)-10s %(processName)s [%(filename)s:%(lineno)s - %(funcName)20s() ] %(name)s %(message)s')
logger = logging.getLogger('main')
logger.setLevel(logging.DEBUG)
# logger.setLevel(logging.INFO)
# logger.setLevel(logging.ERROR)
PHI = 1.61803398874989484820458683436563811772
PHI1 = -0.618033988749894902525738871190696954727
SQRT5 = 2.236067977499789696409173668731276235441
maxn = 50
qs = QuadraticSieve(maxn)
for m in range(5, maxn):
    logger.debug("m {}".format(m))
    if not(m&1):
        continue
#    ls = qs.simbolo_jacobi(5, m)
    ls = jacobi_symbol(5,m)
    logger.debug("m {} ls {}".format(m, ls))
    if(ls == -1):
        continue
    if(ls == 1):
        x1, x2 = qs.calcula_conguencia_residuo_cuadratico(5, m)
        logger.debug("m {} x1 {} x2 {}".format(m, x1, x2))
        x = min(x1, x2)
    else:
        x = 0
    if not x:
        continue
    for n in range(5, maxn):
        # fn_mod_m=floor(pow_mod(PHI, n, SQRT5 * m) / SQRT5)
        # fn_mod_m=(floor(pow_mod(PHI, n, m)) * invert(4,m))%m
        # fn_mod_m=(fmod((PHI**3)/SQRT5,m)*pow_mod(PHI,n-3,m))
        invert_2 = invert(2, m)
        Phi = (1 + x) * invert_2
        phi = (1 - x) * invert_2
        x_inv = invert(x, m) if x else 0
        phis = pow(Phi, n, m) - pow(phi, n, m)
        logger.debug("para n {} m {} phis {}".format(n, m, phis))
        fn_mod_m = (phis * x_inv) % m
        aaa = floor((PHI ** n) / SQRT5 - (PHI1 ** n / SQRT5)) % m
        assert fn_mod_m == aaa, "calc {} real {} para F({})mod{}".format(fn_mod_m, aaa, n, m)
