from typing import List, Dict, Any
import os, csv
from lfsr import LFSR
from generators import shrinking, geffe, majority
from bm import berlekamp_massey
from tests_nist import (
    monobit, block_frequency, runs_test, cusum_test,
    approx_entropy, maurer_test, serial_test
)

def evaluate(seq: List[int]) -> Dict[str, Any]:
    print(f"  > Calculando Berlekamp–Massey para {len(seq)} bits...")
    L, _ = berlekamp_massey(seq)

    print("  > Ejecutando tests NIST...")
    s1, s2 = serial_test(seq, 3)

    return {
        "lin_complexity": L,
        "monobit": monobit(seq),
        "blockfreq_128": block_frequency(seq, 128),
        "runs": runs_test(seq),
        "cusum_fwd": cusum_test(seq, True),
        "cusum_rev": cusum_test(seq, False),
        "approx_entropy_m2": approx_entropy(seq, 2),
        "maurer_L6": maurer_test(seq, 6),
        "serial_m3_p1": s1,
        "serial_m3_p2": s2,
    }

def main():
    os.makedirs("results", exist_ok=True)
    n = 100000   # empieza con menos para debug

    print("Generando secuencias con LFSR...")
    L1 = LFSR(taps=[19,5,1], state=[1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,0,1,0,1])
    L2 = LFSR(taps=[19,6,2], state=[1,1,0,1,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1])
    L3 = LFSR(taps=[19,9,1], state=[1,0,1,1,1,0,0,0,1,1,0,1,0,1,0,1,1,1,0])

    a = L1.run_k(n)
    b = L2.run_k(n)
    c = L3.run_k(n)

    cases = [
        ("L1", a),
        ("L2", b),
        ("L3", c),
        ("shrink_L1_L2", shrinking(a, b)),
        ("geffe_L1_L2_L3", geffe(a, b, c)),
        ("maj_L1_L2_L3", majority(a, b, c)),
    ]

    cols = ["name","lin_complexity","monobit","blockfreq_128","runs",
            "cusum_fwd","cusum_rev","approx_entropy_m2","maurer_L6",
            "serial_m3_p1","serial_m3_p2"]

    print("Evaluando casos...")
    with open("results/summary.csv","w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for name, seq in cases:
            print(f"\nEvaluando {name}...")
            res = evaluate(seq)
            res["name"] = name
            w.writerow({k: res.get(k,"") for k in cols})

    print("\n Resultados guardados en results/summary.csv")

if __name__ == "__main__":
    main()
