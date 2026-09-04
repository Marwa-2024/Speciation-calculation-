# Dolomite surface titration — speciation and surface charge

Reproduces the batch surface titration of dolomite in NaCl solutions following
**Pokrovsky et al. (1999)**. The aqueous speciation is solved with Davies activity
coefficients and Plummer–Busenberg / Nordstrom formation constants, and the net
surface charge is the charge-sum difference of Charlet et al. (1990) — no surface
complexation model is fitted or assumed.

## Files

| File | Purpose |
|------|---------|
| `DolomiteSurfaceTitration.ipynb` | Google Colab / Jupyter notebook (Python). The executed, verified reference. |
| `DolomiteSurfaceTitration.nb` | Mathematica notebook, mirrors the Python. |
| `DolomiteSurfaceTitration.wl` | The same Mathematica code as a plain package/script. |

### Running in Google Colab

Open [colab.research.google.com](https://colab.research.google.com), choose **File ▸ Upload
notebook**, pick `DolomiteSurfaceTitration.ipynb`, then **Runtime ▸ Run all**. It uses only
numpy, pandas and matplotlib, which ship with Colab. The Mathematica files mirror the same
method; the Python notebook is the version that has been executed and verified.

## Method

**Surface charge (Pokrovsky Eqn 1).** The net surface charge is the charge-sum difference
between the conservative A+B mix (the zero-charge reference) and the actual reactor C:

```
σ_T = ( Σ z_k [k]_0  −  Σ z_k [k]_f ) / S
```

`[k]_0` is the closed 1:1 mix of reactors A and B; `[k]_f` is reactor C, both speciated at the
**measured** pH; `S` is the dolomite surface area per litre of C (0.84 m²/g × 30 g/L = 25.2 m²/L).
Free Na⁺ and Cl⁻ are omitted from the charge sum — they are conserved on mixing and cancel in the
difference.

**Solved step by step.** Each stage is a separate, labeled block in the notebook:

1. **Activity coefficients (Davies).** `log γ_z = −A z² ( √I/(1+√I) − 0.3 I )`, iterated on ionic strength.
2. **Law of mass action, γ-substituted.** Every complex from its formation constant in activities,
   then converted to concentration, e.g. `[CaCO₃] = 10^K · γ₂² · [Ca²⁺][CO₃²⁻]`.
3. **Free ions from the mass balances.** `[Ca²⁺] = Ca_T / d_Ca` with `d_Ca` the sum over all Ca species.
4. **Charge balance.** `Σ z_k [k]`, grouped as `q_H + q_DIC + q_Ca + q_Mg`.
5. **Surface charge.** Eqn 1 above.

**Carbonate switch.** `CARBONATE = 'closed'` (default) uses the measured alkalinity, as in the paper,
and keeps carbonate near the measured ~1–2 mM. `CARBONATE = 'open'` fixes CO₂(aq) by Henry's law at
atmospheric pCO₂. Open pCO₂ predicts hundreds of mM of carbonate above pH 9 and the base-run surface
charge blows up to ±10–25 mmol/m², so `'closed'` is the default; set `'open'` to reproduce the breakdown.

**Measured pH.** Pokrovsky computes `[k]_f` from the measured pH, alkalinity and total Ca/Mg, so the
speciation uses the electrode pH directly — it is the observable that encodes the surface reaction.

## Vessels and data quality

- **A** — dolomite + NaCl, equilibrated (Ca, Mg released by the solid); the zero-charge reference.
- **B** — NaCl blank, pH pre-adjusted (runs 1–4 HCl, run 5 none, runs 6–9 NaOH, in cumulative
  0.001 mol/L steps).
- **C** — the A + B mixture after re-equilibration.

Two conservation budgets flag which runs are real titrations:

- **Acid budget.** Dissolution by acid consumes ~2 protons per divalent cation pair. Runs 1–4 release
  20–40× more divalent cation than the 0.001–0.004 mol/L of HCl can explain (run 4: 2 meq/L acid vs
  83 meq/L cations), so those points are **not** acid-driven titrations and are dropped from the curve.
- **Carbon budget.** The carbon released in runs 1–4 is likewise unaccounted for. Runs 5–9 close both
  budgets to within about a factor of two, so the titration curve is built from runs 5–9, plotted
  against the measured pH and compared to the 0.0166 mmol/m² (~10 sites/nm²) site-density ceiling.

## Notes

1. **No surface complexation model** — the surface charge is the aqueous charge-sum difference of
   Eqn 1, from the measured cations, alkalinity and pH only.
2. **Activity coefficients (Davies)** and Plummer–Busenberg / Nordstrom formation constants, as in
   the paper; ionic strength is iterated to convergence.
3. **Open vs closed carbonate** is a one-line switch; closed (measured alkalinity) is the default
   because open atmospheric pCO₂ is inconsistent with the measured DIC above pH 9.
4. The Mathematica `.nb`/`.wl` mirror the Python method but have not been executed in-engine here;
   the Python notebook is the verified reference.
