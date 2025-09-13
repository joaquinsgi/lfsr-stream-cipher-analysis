from typing import List, Tuple

def berlekamp_massey(s: List[int]) -> Tuple[int, List[int]]:
    c = [1]
    b = [1]
    l = 0
    m = -1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, l+1):
            d ^= c[i] & s[n-i]
        if d == 1:
            t = c[:]
            pad = [0]*(n-m)
            c = c + [0]*max(0, len(b)+len(pad)-len(c))
            for i in range(len(b)):
                c[i+len(pad)] ^= b[i]
            if 2*l <= n:
                l = n+1-l
                b = t
                m = n
    return l, c
