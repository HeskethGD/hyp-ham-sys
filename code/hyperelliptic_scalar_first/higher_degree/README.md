# Higher-degree polynomial corrections and determinants

Open [the marimo notebook](notebook.py) or its [executed HTML](notebook.html).
The [derivation](proof.md) gives arbitrary-genus differential reduction and
the longer determinant interpolation argument; [results.json](results.json)
contains every numerical and symbolic check.

```sh
PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.higher_degree.verification
.venv/bin/python -m marimo edit code/hyperelliptic_scalar_first/higher_degree/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/higher_degree/notebook.py -o code/hyperelliptic_scalar_first/higher_degree/notebook.html --force
```

The genus-two example tests quartic restricted q, degree-five corrections,
and a seven-row sigma determinant containing xy/2. The genus-three example
tests degree-eight corrections while keeping restricted q constant.
Independent comparisons use direct quadrature and the full Hamiltonian
field and Jacobian. Both moving sheets are tested.

At 40 working digits, all 629 exact and numerical checks pass. The largest
scaled numerical errors are 2.19e-35 (genus two) and 1.70e-34 (genus three),
against a threshold of 1e-26. The executed HTML reruns the verification.

The new endpoint and second-kind corrections cancel in the complete mode
products. The longer determinant needs separate prefactors for its two
fixed-sheet choices. Moving second-kind primitives are evaluated from a
finite initial place so their lifts agree on both sheets.

`model.py` contains the exact differential reduction and scalar models;
`verification.py` runs the algebraic and sigma checks. The package reuses
the main sigma-mode evaluator and the shared `abel_lifts.py` helper, with no
dependency on exploratory files. Even-degree curves, singular curves,
confluent determinants and global continuation remain outside these checks.
