import numpy as np
from typing import List

class LFSR:
    def __init__(self, taps: List[int], state: List[int]):
        self.taps = np.array([-t for t in taps], dtype=int)
        s = np.array(state, dtype=np.uint8)
        if s.ndim != 1:
            raise ValueError("state must be 1-D")
        if s.size <= 0 or s.size != max(taps):
            raise ValueError("state size must equal max(taps)")
        if not np.any(s):
            raise ValueError("state must not be all zeros")
        self.state = s

    def next_bit(self) -> int:
        out = int(self.state[-1])
        fb = int(np.bitwise_xor.reduce(self.state[self.taps]))
        self.state[1:] = self.state[:-1]
        self.state[0] = fb
        return out

    def run_k(self, k: int) -> list:
        out = np.empty(k, dtype=np.uint8)
        st = self.state
        taps = self.taps
        for i in range(k):
            out[i] = st[-1]
            fb = int(np.bitwise_xor.reduce(st[taps]))
            st[1:] = st[:-1]
            st[0] = fb
        return out.tolist()
