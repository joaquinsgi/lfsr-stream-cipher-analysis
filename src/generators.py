from typing import List

def shrinking(a: List[int], b: List[int]) -> List[int]:
    n = min(len(a), len(b))
    s = []
    for i in range(n):
        if a[i] == 1:
            s.append(b[i])
    return s

def geffe(a: List[int], b: List[int], c: List[int]) -> List[int]:
    n = min(len(a), len(b), len(c))
    s = []
    for i in range(n):
        if a[i] == 1:
            s.append(b[i])
        else:
            s.append(c[i])
    return s

def majority(*seqs: List[int]) -> List[int]:
    n = min(len(s) for s in seqs)
    s = []
    m = len(seqs)
    for i in range(n):
        t = sum(seqs[j][i] for j in range(m))
        s.append(1 if t >= (m+1)//2 else 0)
    return s
