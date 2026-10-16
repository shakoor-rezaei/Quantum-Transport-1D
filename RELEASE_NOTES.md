# Release notes — numerical research portfolio, v1.0

This repository combines three original MSc model families and two extensions
in a single numerical engine. Full model conventions, compilation and example
commands are in `README.md`.

Validation performed on this codebase:

- Clean chain and one-impurity analytic tests: pass.
- Binary, exponential, zeta and both generalized weighting rules: pass unit tests.
- Seeded outputs invariant across one vs two threads: pass.
- Four independent full-length scans: N=10,000, 1,000 disorder realizations per
  energy, 77 energies each. Model parameters match the published standard cases.
- Recorded two-column MSc results and newly generated summaries are compared
  under `validation/`; all individual numeric pairings are archived as CSV.

The median absolute differences outside specified problematic energies are
2.90%, 3.00%, 1.21%, and 1.30% respectively for binary (2,5), binary (3,6),
exponential and standard zeta. The discrepancy near nominal E=±1 in binary (3,6) is largely explained by the finite-chain fitting length: historical 5000-point fits and new N=5000 fits agree within 0.15%. The N=10000 fits differ, so the asymptotic localization length is not established; see validation/resonance_binary_3_6/README.md.
 Binary (2,5) resonance lengths have extremely
small fit R² and should not be taken as robust finite values.

These checks establish meaningful tested numerical agreement, not identity to
the original stochastic samples or proof across all possible model parameters.
The thesis figures and article's fourth-order perturbation calculation are
not recreated directly in this release. The included analytic comparison script
is second order only. The generalized cases are model extensions and are not
claimed to be research findings of the original article.
