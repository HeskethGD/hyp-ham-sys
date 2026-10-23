# Simple branch fixed points

Symbolic and numerical checks supporting the main paper's shifted sigma
mode formulas at simple branch fixed points, including continuation through
a turning point and the determinant identities.

From the repository root:

```sh
PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.branch_points.verification
.venv/bin/python -m marimo edit code/hyperelliptic_scalar_first/branch_points/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/branch_points/notebook.py -o code/hyperelliptic_scalar_first/branch_points/notebook.html --force
PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.branch_points.determinants
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/branch_points/determinant_notebook.py -o code/hyperelliptic_scalar_first/branch_points/determinant_notebook.html --force
```

The executable notebook runs the checks and saves `results.json`. The CLI
supports `--genus 2`, `--genus 3`, and `--dps`. The models reuse the main
package's sigma jets and canonical algebra. `model.py` supplies an endpoint-path
subclass that regularizes branch endpoints; `verification.py` checks the modes
and `determinants.py` checks the corresponding determinant identities. No Generapy change is needed.

Coverage includes compatible first- and second-kind endpoint integrals,
fixed sigma normalizers and stratum identities, independent u/v sigma ratios,
canonical equations, polynomial potentials, a non-special Q on the same leaf,
and signed continuation through a turning point, including the zero itself.
Exact checks cover the repeated-weight singularity obstruction.

## Arbitrary-genus justification and determinant extension

Read [proof.md](proof.md) for the arbitrary-genus branch-point continuation
argument and the extension of the two existing determinant families. The
proof cites the precise Onishi and prime-form results it uses and separates
its scope from the numerical evidence.

The [determinant notebook](determinant_notebook.py), [executed HTML](determinant_notebook.html)
and [saved results](determinant_results.json) check the full sigma determinant,
auxiliary residual splitting, period balancing and invariance of the corrected
modes. Numerical sigma checks cover genera two and three; exact mixed-row
interpolation checks cover genera two through six and both column families.

Global monodromy, a second independent solution, singular curves and
confluent determinant extensions remain outside these checks. Longer column
lists are tested separately in [higher_degree](../higher_degree/README.md).
