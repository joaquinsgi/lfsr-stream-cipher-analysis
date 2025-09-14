# Thesis – Study and Optimization of Pseudo-Random Sequences Based on LFSR

This repository contains the Bachelor’s Thesis of **Joaquín Sergio García Ibáñez**, titled:

> *Study and Optimization of Pseudo-Random Bit Sequences Based on LFSR for Stream Ciphers*

The project combines **theoretical analysis in LaTeX** and **practical experiments in Python** to evaluate pseudo-random sequence generators based on Linear Feedback Shift Registers (LFSR) and classical combiners such as Shrinking, Geffe, and Majority.
It studies their statistical properties, linear complexity, and resistance against NIST randomness tests.

---

## 📂 Repository Structure

```
.
├── src/                # Python source code
│   ├── lfsr.py             # LFSR implementation
│   ├── generators.py       # Shrinking, Geffe, and Majority generators
│   ├── bm.py               # Berlekamp–Massey algorithm
│   ├── tests_nist.py       # Selected NIST statistical tests
│   ├── experiment.py       # Single experiment execution
│   └── experiment_multi.py # Multiple experiments execution
│
├── results/            # Results in CSV format
│   ├── summary.csv
│   └── summary_multi.csv
│
├── TFG/                # LaTeX source of the thesis
│   ├── chapters/           # Thesis chapters in .tex files
│   ├── figures/            # Figures used in the document
│   ├── main.tex            # Main LaTeX file
│   └── referencias.bib     # Bibliography
│
├── data/               # (Optional) auxiliary data for experiments
├── requirements.txt    # Python dependencies
└── README.md           # This file
```

---

## ▶️ How to Use

### 1. Create a virtual environment and install dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Run a single experiment

```bash
python src/experiment.py
```

### 3. Run multiple experiments and save results

```bash
python src/experiment_multi.py
```

---

## 📊 LaTeX Document

The full thesis document is available inside the `TFG/` directory.
To compile the PDF:

```bash
cd TFG
pdflatex main.tex
biber main
pdflatex main.tex
pdflatex main.tex
```

The compiled file will be available as `TFG/main.pdf`.

---

## 📜 License

This work is published under the MIT License for the code and with academic purposes for the LaTeX document.