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
to install. The Python notebook and the Mathematica notebook implement the same method and
give the same results.

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
neg = OH⁻ + HCO₃⁻ + 2·CO₃²⁻ + 2·NaCO₃⁻ + NaHCO₃
```

The background electrolyte (free Na⁺, Cl⁻) is deliberately left out of `pos`/`neg`;
its imbalance is what the mineral surface carries.

**Surface charge.** `σ = pos − neg`. It is zero at the calculated pH (no surface
interaction) and becomes negative or positive when the real, measured pH lies
above or below that value. The notebook plots σ(pH) for vessels A and C and the
per-sample surface charge at the measured pH, normalised per gram of dolomite
(6 g of solid; 100 mL in A, 200 mL in C).

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
3. No surface complexation constants are used — the surface charge is purely a
   charge-balance quantity.
