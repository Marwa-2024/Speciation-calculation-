# Dolomite in NaCl — aqueous speciation and net surface charge

Computes the aqueous speciation of the dolomite surface-titration solutions from **published
thermodynamic constants only**, then derives the net surface charge of dolomite
(Pokrovsky et al., 1999) from that speciation. No surface complexation model is fitted, and no value
is taken from any software output — a Visual MINTEQ run of the same solutions is used for comparison
only.

## Files

| File | Purpose |
|------|---------|
| `DolomiteSurfaceTitration.ipynb` | Google Colab / Jupyter notebook (Python). The executed, verified reference. |
| `DolomiteSurfaceTitration.nb` | Mathematica notebook, mirrors the Python; every equation written out. |
| `DolomiteSurfaceTitration.wl` | The same Mathematica code as a plain package/script. |

## Method

- **Published constants.** All log K values are from the TOUGHREACT / EQ3-6 thermodynamic database
  (Plummer & Busenberg 1982; Nordstrom et al. 1990; Wolery 1992). Nothing is taken from the software.
- **Open carbonate, fixed pCO₂.** CO₂(aq) is held at a fixed activity set by atmospheric pCO₂
  (10⁻³·⁵ atm) and Henry's law (log K_H = −1.469), so HCO₃⁻ and CO₃²⁻ follow from pCO₂ and pH. The
  dissolved inorganic carbon is an output, not a fixed input.
- **pH calculated from the mass balance**, not the electrode. The free H⁺ starts from a neutral guess
  (10⁻⁷) and is solved from the proton condition `P(pH) = −2·C_CO3`, where `C_CO3` is the entered total
  CO₃²⁻ (the measured carbonate). Because Ca²⁺, Mg²⁺, Na⁺ and Cl⁻ carry no proton, dissolved Ca/Mg do
  not drag the pH the wrong way — vessel C falls with acid and rises with base.
- **Davies activity coefficients** with A = 0.509 and salting term b = 0.5; neutral species take
  γ₀ = 1.
- **Every equation is listed** in both notebooks: the reactions with their log K, the law of mass
  action for each species, the four element mass balances (Ca, Mg, Na, Cl), and the proton condition.

The calculated pH and speciation are then **compared** with a Visual MINTEQ run of the same solutions.
On the base runs the two agree closely; on the acid runs they differ by roughly 20%, the expected
offset between two independent constant sets where the charge sum is dominated by the large Ca/Mg
release.

## Species

H⁺, OH⁻, CO₂(aq), HCO₃⁻, CO₃²⁻; Ca²⁺, CaCl⁺, CaCO₃(aq), CaHCO₃⁺, CaOH⁺; Mg²⁺, MgCl⁺, MgCO₃(aq),
MgHCO₃⁺, MgOH⁺; Na⁺, NaCl(aq), NaCO₃⁻, NaHCO₃(aq), NaOH(aq); Cl⁻.

## Surface charge (Pokrovsky Eqn 1)

At the calculated pH, `σ_T = ( ½(q_A + q_B) − q_C ) / S`, with `q = Σ z_k [k]` over the reactive
species (free Na⁺, Cl⁻ cancel in the difference) and `S` the dolomite area per litre of reactor C
(0.76 m²/g × 60 g/L = 45.6 m²/L).

## Notes

1. The Python notebook is executed and verified. The Mathematica `.nb`/`.wl` mirror the same equations
   but have not been run in-engine here.
2. The acid runs (pH ≈ 5.5) are dominated by dolomite dissolution, so their large negative σ_T is a
   dissolution signal more than a surface-protonation signal, consistent with the pH ≈ 6.5–11.5 working
   range of Pokrovsky et al. (1999).
3. Vessel B run 6 differs from the MINTEQ comparison because that MINTEQ run carried no carbonate; the
   model value is what a B run with 48 mg/L CO₃ gives.

## References

- Charlet, L., Wersin, P., Stumm, W. (1990). *Geochim. Cosmochim. Acta* 54, 2329–2336.
- Davies, C. W. (1962). *Ion Association*. Butterworths, London.
- Nordstrom, D. K. et al. (1990). Revised chemical equilibrium data for major water–mineral reactions.
- Plummer, L. N., Busenberg, E. (1982). *Geochim. Cosmochim. Acta* 46, 1011–1040.
- Pokrovsky, O. S., Schott, J., Thomas, F. (1999). *Geochim. Cosmochim. Acta* 63, 3133–3143.
- Wolery, T. J. (1992). EQ3/6 thermodynamic database (distributed with TOUGHREACT).
