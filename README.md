# Dolomite surface titration — speciation and surface charge

A Wolfram (Mathematica) notebook that reproduces the batch surface titration of
dolomite in NaCl solutions following Pokrovsky et al. (1999). It solves the
aqueous speciation, calculates pH from the charge balance at a fixed CO₂
partial pressure, and derives the surface charge directly from the charge
balance — no surface complexation model is fitted or assumed.

## Files

| File | Purpose |
|------|---------|
| `DolomiteSurfaceTitration.ipynb` | Google Colab / Jupyter notebook (Python). Open in Colab and run all. |
| `DolomiteSurfaceTitration.nb` | The notebook to open in Mathematica. |
| `DolomiteSurfaceTitration.wl` | The same Mathematica code as a plain package/script. |

### Running in Google Colab

Open [colab.research.google.com](https://colab.research.google.com), choose **File ▸ Upload
notebook**, pick `DolomiteSurfaceTitration.ipynb`, then **Runtime ▸ Run all**. Every
library it uses (numpy, scipy, matplotlib, pandas) ships with Colab, so there is nothing
to install. The Python notebook carries the current method (corrected charge coefficients and
the proton mass balance surface charge); the Mathematica files still hold the earlier
charge-balance version and have not yet been updated to match.

## What is calculated

**Unknowns (free ion concentrations):** H⁺, Ca²⁺, Mg²⁺, Na⁺, Cl⁻. Carbonate
(CO₂(aq), HCO₃⁻, CO₃²⁻) is set by the fixed pCO₂ together with H⁺, so the
system is open with respect to CO₂.

**Law of mass action.** Every complex is written from its dissociation constant
(25 °C, taken from the supplied TOUGHREACT/EQ3-6 database):

- Carbonate: CO₂(aq)+H₂O = H⁺+HCO₃⁻ (log K = −6.345); HCO₃⁻ = H⁺+CO₃²⁻ (−10.329);
  H₂O = H⁺+OH⁻ (−13.995); [CO₂(aq)] = 10^(−1.468)·pCO₂.
- Calcium: CaCl⁺, CaCO₃(aq), CaHCO₃⁺, CaOH⁺.
- Magnesium: MgCl⁺, MgCO₃(aq), MgHCO₃⁺, MgOH⁺.
- Sodium: NaCO₃⁻, NaHCO₃(aq), NaCl(aq), NaOH(aq).

**Mass balances** are enforced for Ca, Mg, Na and Cl, using the measured totals
(ppm converted to mol/L).

**pH from the charge balance.** For each sample the notebook finds the pH at
which the reactive charge balance closes:

```
pos = H⁺ + 2·Ca²⁺ + 2·Mg²⁺ + CaHCO₃⁺ + CaOH⁺ + MgHCO₃⁺ + MgOH⁺
neg = OH⁻ + HCO₃⁻ + 2·CO₃²⁻ + NaCO₃⁻
```

The background electrolyte (free Na⁺, Cl⁻) is deliberately left out of `pos`/`neg`;
its imbalance is what the mineral surface carries. NaCO₃⁻ carries a single negative
charge and NaHCO₃(aq) is neutral, so it drops out of the sum. Every pH is calculated;
the electrode reading is kept only for comparison.

**Surface charge (proton mass balance).** The calculated pH is fixed by the charge
balance above, so that same sum is zero for every solution and cannot itself be the
surface charge. Instead the surface charge is the net proton exchange between the
conservative A+B mixture and the reacted suspension C:

```
σ_H = ( TOTH_mix − TOTH_C ) / St ,     TOTH_mix = ( TOTH_A + TOTH_B ) / 2
```

`TOTH` is the total proton excess of a solution relative to the reference components
(H₂O, CO₃²⁻, Ca²⁺, Mg²⁺, Na⁺, Cl⁻). The fixed-pCO₂ term is identical in every solution
and cancels in the difference, so the open CO₂ reservoir contributes nothing to the
surface charge. `St` is the surface area per litre of vessel C (specific area × loading,
0.84 m²/g × 30 g/L = 25.2 m²/L). The notebook reports σ_H per m² and plots it against
the calculated pH of vessel C: positive as pH drops, negative at high pH.

## Vessels

- **A** — dolomite + NaCl, equilibrated (Ca, Mg released by the solid).
- **B** — NaCl blank, pH pre-adjusted (no Ca, no Mg; the calculated pH reflects
  only the water/carbonate system because the strong acid/base titrant is not a
  species in the balance).
- **C** — the A + B mixture after re-equilibration; this vessel gives the clearest
  surface-charge trend, positive at low pH and negative at high pH.

## Assumptions and knobs

1. Concentrations are used directly (activity coefficients = 1). At the high
   ionic strength of these NaCl solutions a Davies/SIT/Pitzer correction would
   shift the numbers; add γ factors inside `allSpecies` if required.
2. `pCO2` is a single parameter near the top of the notebook (default 10^(−3.5)
   atm). To use the measured alkalinity instead, replace the carbonate lines in
   `allSpecies` with a carbonate mass balance.
3. No surface complexation constants are used — the surface charge is a proton
   mass balance derived from the measured cations and the aqueous speciation only.
4. Forcing an open system together with a calculated pH detaches the model pH from
   the electrode (large dissolved Ca/Mg in equilibrium with atmospheric CO₂ must be
   alkaline), so the calculated pH is an internal reference rather than a prediction
   of the measured value. Runs 1–4 show strong dissolution and their σ_H is dominated
   by carbonate released into solution, not by simple surface protonation.
