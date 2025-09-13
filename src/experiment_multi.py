import os, csv, math
import numpy as np
from typing import List, Dict, Any
from src.lfsr import LFSR
from src.generators import shrinking, geffe, majority
from src.bm import berlekamp_massey
from src.tests_nist import (
    monobit, block_frequency, runs_test, cusum_test,
    approx_entropy, maurer_test, serial_test
)

def eval_one(seq: List[int]) -> Dict[str, float]:
    print(f"    > Berlekamp–Massey en {len(seq)} bits (usando 20000 prefix)...")
    L, _ = berlekamp_massey(seq[:20000])
    print("    > Ejecutando tests NIST...")
    s1, s2 = serial_test(seq, 3)
    return {
        "len": float(len(seq)),
        "lin_complexity": float(L),
        "monobit": float(monobit(seq)),
        "blockfreq_128": float(block_frequency(seq, 128)),
        "runs": float(runs_test(seq)),
        "cusum_fwd": float(cusum_test(seq, True)),
        "cusum_rev": float(cusum_test(seq, False)),
        "approx_entropy_m2": float(approx_entropy(seq, 2)),
        "maurer_L6": float(maurer_test(seq, 6)),
        "serial_m3_p1": float(s1),
        "serial_m3_p2": float(s2),
    }

def pass_flags(m: Dict[str, float], alpha: float = 0.01) -> Dict[str, int]:
    keys = ["monobit","blockfreq_128","runs","cusum_fwd","cusum_rev",
            "approx_entropy_m2","maurer_L6","serial_m3_p1","serial_m3_p2"]
    return {k: int(m[k] >= alpha) for k in keys}

def aggregate(rows: List[Dict[str, float]]) -> Dict[str, Any]:
    keys = ["len","lin_complexity","monobit","blockfreq_128","runs","cusum_fwd",
            "cusum_rev","approx_entropy_m2","maurer_L6","serial_m3_p1","serial_m3_p2"]
    out = {}
    for k in keys:
        v = np.array([r[k] for r in rows], dtype=float)
        out[f"{k}_mean"] = float(v.mean())
        out[f"{k}_std"] = float(v.std(ddof=1)) if v.size > 1 else 0.0
    flags = [pass_flags(r) for r in rows]
    for k in ["monobit","blockfreq_128","runs","cusum_fwd","cusum_rev",
              "approx_entropy_m2","maurer_L6","serial_m3_p1","serial_m3_p2"]:
        out[f"{k}_pass_rate"] = float(np.mean([f[k] for f in flags]))
    return out

def random_state(rng: np.random.Generator, n: int) -> List[int]:
    s = rng.integers(0, 2, size=n, dtype=np.uint8)
    while not np.any(s):
        s = rng.integers(0, 2, size=n, dtype=np.uint8)
    return s.tolist()

def main():
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(12345)
    nbits = 100000
    reps = 5

    taps1 = [19,5,1]
    taps2 = [19,6,2]
    taps3 = [19,9,1]
    deg = 19

    results = {}
    for r in range(reps):
        print(f"\n=== Repetición {r+1}/{reps} ===")
        s1 = random_state(rng, deg)
        s2 = random_state(rng, deg)
        s3 = random_state(rng, deg)

        print("  Generando secuencias LFSR...")
        L1 = LFSR(taps=taps1, state=s1)
        L2 = LFSR(taps=taps2, state=s2)
        L3 = LFSR(taps=taps3, state=s3)
        a = L1.run_k(nbits)
        b = L2.run_k(nbits)
        c = L3.run_k(nbits)

        cases = {
            "L1": a,
            "L2": b,
            "L3": c,
            "shrink_L1_L2": shrinking(a, b),
            "geffe_L1_L2_L3": geffe(a, b, c),
            "maj_L1_L2_L3": majority(a, b, c),
        }

        for name, seq in cases.items():
            print(f"\n  Evaluando {name}...")
            res = eval_one(seq)
            results.setdefault(name, []).append(res)

    cols = [
        "name",
        "len_mean","len_std",
        "lin_complexity_mean","lin_complexity_std",
        "monobit_mean","monobit_std","monobit_pass_rate",
        "blockfreq_128_mean","blockfreq_128_std","blockfreq_128_pass_rate",
        "runs_mean","runs_std","runs_pass_rate",
        "cusum_fwd_mean","cusum_fwd_std","cusum_fwd_pass_rate",
        "cusum_rev_mean","cusum_rev_std","cusum_rev_pass_rate",
        "approx_entropy_m2_mean","approx_entropy_m2_std","approx_entropy_m2_pass_rate",
        "maurer_L6_mean","maurer_L6_std","maurer_L6_pass_rate",
        "serial_m3_p1_mean","serial_m3_p1_std","serial_m3_p1_pass_rate",
        "serial_m3_p2_mean","serial_m3_p2_std","serial_m3_p2_pass_rate",
    ]

    out_path = "results/summary_multi.csv"
    print(f"\nGuardando resultados en {out_path}...")
    with open(out_path,"w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for name, lst in results.items():
            agg = aggregate(lst)
            row = {"name": name}
            for k in cols[1:]:
                row[k] = agg.get(k, "")
            w.writerow(row)

    print("\n Finalizado. Resultados guardados en results/summary_multi.csv")

if __name__ == "__main__":
    main()
