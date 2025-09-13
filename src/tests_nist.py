import math
import numpy as np
from itertools import product
from typing import List, Tuple
from scipy.special import erfc, gammaincc

def monobit(seq: List[int]) -> float:
    n = len(seq)
    s = sum(1 if x == 1 else -1 for x in seq)
    return erfc(abs(s) / math.sqrt(2*n))

def block_frequency(seq: List[int], M: int = 128) -> float:
    n = len(seq)
    N = n // M
    if N == 0:
        return 0.0
    x = np.array(seq[:N*M], dtype=int).reshape(N, M)
    prop = np.mean(x, axis=1)
    chi = 4*M*np.sum((prop-0.5)**2)
    return math.exp(-chi/2)

def runs_test(seq: List[int]) -> float:
    n = len(seq)
    pi = sum(seq)/n
    if abs(pi-0.5) >= (2/math.sqrt(n)):
        return 0.0
    v = 1 + sum(1 for i in range(n-1) if seq[i] != seq[i+1])
    num = abs(v - 2*n*pi*(1-pi))
    den = 2*math.sqrt(2*n)*pi*(1-pi)
    return erfc(num/den)

def cusum_test(seq: List[int], forward: bool = True) -> float:
    n = len(seq)
    x = np.array([1 if b == 1 else -1 for b in seq], dtype=int)
    if not forward:
        x = x[::-1]
    s = np.cumsum(x)
    z = max(abs(s))
    k = math.floor((n/z - 1)/4)
    s1 = sum(math.exp(-((4*j+1)**2)*(z**2)/(2*n)) for j in range(-k, k+1))
    k = math.floor((n/z - 3)/4)
    s2 = sum(math.exp(-((4*j+3)**2)*(z**2)/(2*n)) for j in range(-k, k+1))
    return 1.0 - s1 + s2

def approx_entropy(seq: List[int], m: int = 2) -> float:
    n = len(seq)
    def phi(mm: int) -> float:
        v = seq + seq[:mm-1]
        counts = {}
        for i in range(n):
            key = tuple(v[i:i+mm])
            counts[key] = counts.get(key, 0) + 1
        return sum((c/n)*math.log(c/n) for c in counts.values())
    ap_en = phi(m) - phi(m+1)
    chi = 2*n*(math.log(2) - ap_en)
    return gammaincc(2**(m-1), chi/2)

def maurer_test(seq: List[int], L: int = 6) -> float:
    n = len(seq)
    Q = 10*(2**L)
    K = n//L - Q
    if K <= 0:
        return 0.0
    def to_int(b):
        v = 0
        for bit in b:
            v = (v<<1) | bit
        return v
    T = [0]*(2**L)
    v = seq[:Q*L]
    for i in range(Q):
        blk = v[i*L:(i+1)*L]
        T[to_int(blk)] = i+1
    s = 0.0
    w = seq[Q*L:(Q+K)*L]
    for i in range(K):
        blk = w[i*L:(i+1)*L]
        j = to_int(blk)
        d = i+Q+1 - T[j]
        T[j] = i+Q+1
        s += math.log2(d)
    fn = s/K
    exp_var = {
        6:(5.2177052,2.954),7:(6.1962507,3.125),8:(7.1836656,3.238),
        9:(8.1764248,3.311),10:(9.1723243,3.356),11:(10.170032,3.384),
        12:(11.168765,3.401),13:(12.168070,3.410),14:(13.167693,3.416),
        15:(14.167488,3.419),16:(15.167379,3.421)
    }
    if L not in exp_var:
        return 0.0
    ev, var = exp_var[L]
    c = 0.7 - 0.8/L + (4 + 32/L)*(K**(-3/L))/15
    sigma = c*math.sqrt(var/K)
    return erfc(abs(fn-ev)/(math.sqrt(2)*sigma))

def serial_test(seq: List[int], m: int = 3) -> Tuple[float, float]:
    n = len(seq)
    def count_blocks(mm: int) -> float:
        v = seq + seq[:mm-1]
        counts = [0]*(2**mm)
        val = 0
        mask = (1<<mm)-1
        for i in range(mm):
            val = ((val<<1) | v[i]) & mask
        counts[val] += 1
        for i in range(1, n):
            val = ((val<<1) | v[i+mm-1]) & mask
            counts[val] += 1
        return (2**mm)/n * sum(c*c for c in counts) - n
    psi_m = count_blocks(m)
    psi_m1 = count_blocks(m-1)
    psi_m2 = count_blocks(m-2) if m >= 2 else 0.0
    d1 = psi_m - psi_m1
    d2 = psi_m - 2*psi_m1 + psi_m2
    p1 = gammaincc(2**(m-2), d1/2)
    p2 = gammaincc(2**(m-3), d2/2) if m >= 3 else 0.0
    return p1, p2
