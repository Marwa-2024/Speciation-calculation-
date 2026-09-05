# Dolomite in NaCl — Visual MINTEQ-style speciation and surface charge

Reproduces the Visual MINTEQ speciation for the dolomite surface-titration solutions, then computes
the net surface charge (Pokrovsky et al., 1999) from it. No surface complexation model is fitted.

## Files

| File | Purpose |
|------|---------|
| `DolomiteSurfaceTitration.ipynb` | Google Colab / Jupyter notebook (Python). The executed, verified reference. |
| `DolomiteSurfaceTitration.nb` | Mathematica notebook, mirrors the Python; every equation written out. |
| `DolomiteSurfaceTitration.wl` | The same Mathematica code as a plain package/script. |

## Method (matches MINTEQ option 1)

- **Open carbonate, fixed pCO₂.** CO₂(aq) is held at a fixed activity, {CO₂(aq)} = 10⁻⁴·⁹ (MINTEQ's
  reported value), so HCO₃⁻ and CO₃²⁻ follow from pCO₂ and pH.
- **pH calculated from the mass balance**, not the electrode. The free H⁺ starts from a neutral guess
  (10⁻⁷) and is solved from the proton condition `P(pH) = −2·C_CO3`, where `C_CO3` is the entered total
  CO₃²⁻ (the measured carbonate). Because Ca²⁺, Mg²⁺, Na⁺ and Cl⁻ carry no proton, dissolved Ca/Mg do
  not drag the pH the wrong way — vessel C now falls with acid and rises with base, as MINTEQ gives.
- **Davies activity coefficients** (with the 0.1·I Setchenow term for neutrals) and the
  TOUGHREACT / EQ3-6 (Plummer–Busenberg) constants, exactly as in MINTEQ.
- **Every equation is listed** in both notebooks: the reactions with their log K, the law of mass
  action for each species, the four element mass balances (Ca, Mg, Na, Cl), and the proton condition.

The calculated pH reproduces the MINTEQ output to within about 0.1 pH unit across all 27 solutions
(A ≈ 8; B acid runs ≈ 5.5, base runs ≈ 8.3; C acid runs ≈ 5.5, base runs ≈ 8.6). The small residual is
the activity-model detail — MINTEQ's ion pairing gives a slightly lower ionic strength.

## Species

H⁺, OH⁻, CO₂(aq), HCO₃⁻, CO₃²⁻; Ca²⁺, CaCl⁺, CaCO₃(aq), CaHCO₃⁺, CaOH⁺; Mg²⁺, MgCl⁺, MgCO₃(aq),
MgHCO₃⁺, MgOH⁺; Na⁺, NaCl(aq), NaCO₃⁻, NaHCO₃(aq), NaOH(aq); Cl⁻.

## Surface charge (Pokrovsky Eqn 1)

At the calculated pH, `σ_T = ( ½(q_A + q_B) − q_C ) / S`, with `q = Σ z_k [k]` over the reactive
species (Na⁺, Cl⁻ cancel in the difference) and `S` the dolomite area per litre of reactor C
(0.84 m²/g × 30 g/L = 25.2 m²/L).

## Notes

1. The Python notebook is executed and verified against the MINTEQ output. The Mathematica `.nb`/`.wl`
   mirror the same equations but have not been run in-engine here.
2. Vessel B run 6 differs from the MINTEQ sheet because that MINTEQ run carried no carbonate (its
   RUN5/RUN6 columns are identical); the model value is what a B run with 48 mg/L CO₃ gives.
