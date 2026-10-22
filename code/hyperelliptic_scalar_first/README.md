# Computational material

See the repository-root README for installation, execution and paper builds.

Start with `genus_two_example/notebook.py` or `genus_three_example/notebook.py` for the worked examples, and `onishi_products/notebook.py` for the determinant argument. Each has a static HTML walkthrough and saved results.

For the special-polynomial results, use `special_q/notebook.py` for the
weighted classification, Wronskians, polynomial potentials and independent
ODE comparisons. Continue with `sigma_modes/notebook.py` for reconstruction
of the same three examples using sigma functions and analytic first and
second derivatives. Both follow the same organization: `model.py`,
`verification.py`, marimo source, an executed HTML export and `results.json`.
Their READMEs give the corresponding package commands.

The top-level `notebook.py` is a supporting derivation notebook. It contains the fixed-point sigma formulas, regularized second-kind integrals, reference-point coefficient extraction and short regular trajectories on both sheets. Its `verification.py` uses `model.py` and `sigma_jets.py` and reads this notebook's source for formatting checks. These files therefore remain together as one package. The longer worked examples have their own initial states and validation reports.

`paper_plots.py` produces the combined publication figures from `plot_data.json` in the example folders. Those saved samples allow layout changes without rerunning the numerical calculations.
