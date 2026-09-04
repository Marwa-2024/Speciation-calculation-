(* ::Package:: *)

(* ==================================================================== *)
(*  Dolomite surface titration in NaCl  (after Pokrovsky et al., 1999)  *)
(*                                                                      *)
(*  What this notebook does                                             *)
(*   - Writes the law of mass action for every aqueous complex.         *)
(*   - Uses the free ion concentrations as the unknowns                 *)
(*        ( H+, Ca2+, Mg2+, Na+, Cl- ; carbonate set by a fixed pCO2 ). *)
(*   - Solves the mass balances for Ca, Mg, Na, Cl.                     *)
(*   - Calculates pH from the charge balance (it is NOT read from the   *)
(*        measured pH); pCO2 is held fixed.                             *)
(*   - Lists the full speciation for vessels A, B and C.                *)
(*   - Computes the surface charge directly from the charge balance     *)
(*        ( no surface complexation model ) and plots it against pH.    *)
(* ==================================================================== *)

ClearAll["Global`*"];

(* -------------------------------------------------------------------- *)
(*  1.  Molar masses (g/mol) and a ppm (mg/L) -> mol/L converter        *)
(* -------------------------------------------------------------------- *)
mmNa = 22.990;  mmCl = 35.453;  mmCa = 40.078;  mmMg = 24.305;

ppmToM[ppm_, mm_] := ppm/1000./mm;      (* mg/L divided by g/mol -> mol/L *)

(* -------------------------------------------------------------------- *)
(*  2.  Equilibrium constants at 25 C  (log10 K)                        *)
(*      Reactions are written exactly as in the TOUGHREACT / EQ3-6      *)
(*      database supplied (complex on the left, basis species on the    *)
(*      right).  The number is log K of that dissociation.              *)
(* -------------------------------------------------------------------- *)

(* Carbonate system -- open system, pCO2 fixed                          *)
logKH = -1.468;    (* CO2(g)  = CO2(aq)      : [CO2aq] = 10^logKH * pCO2 *)
logK1 = -6.345;    (* CO2(aq) + H2O = H+ + HCO3-                          *)
logK2 = -10.329;   (* HCO3-        = H+ + CO3^2-                          *)
logKw = -13.995;   (* H2O          = H+ + OH-                             *)

(* Metal - ligand complexes ( log K of the dissociation as written )    *)
logKCaCl   =  0.696;   (* CaCl+    = Ca2+ + Cl-        *)
logKCaCO3  =  7.002;   (* CaCO3(aq)= Ca2+ + HCO3- - H+ *)
logKCaHCO3 = -1.047;   (* CaHCO3+  = Ca2+ + HCO3-      *)
logKCaOH   = 12.850;   (* CaOH+    = Ca2+ + H2O - H+   *)

logKMgCl   =  0.135;   (* MgCl+    = Mg2+ + Cl-        *)
logKMgCO3  =  7.350;   (* MgCO3(aq)= Mg2+ + HCO3- - H+ *)
logKMgHCO3 = -1.036;   (* MgHCO3+  = Mg2+ + HCO3-      *)
logKMgOH   = 11.785;   (* MgOH+    = Mg2+ + H2O - H+   *)

logKNaCO3  =  9.815;   (* NaCO3-    = Na+ + HCO3- - H+ *)
logKNaHCO3 = -0.154;   (* NaHCO3(aq)= Na+ + HCO3-      *)
logKNaCl   =  0.777;   (* NaCl(aq)  = Na+ + Cl-        *)
logKNaOH   = 14.180;   (* NaOH(aq)  = Na+ + H2O - H+   *)

(* -------------------------------------------------------------------- *)
(*  3.  Law of mass action for every species                            *)
(*      Given the free ion concentrations h, cCa, cMg, cNa, cCl and a   *)
(*      fixed pCO2, each concentration below follows from its K.        *)
(*      Water activity is taken as 1 (concentration basis, no activity  *)
(*      coefficients -- see the note at the end).                       *)
(* -------------------------------------------------------------------- *)

(* the pCO2 argument is called pc_ so the definition never clashes with the
   global pCO2 value set further down (a symbol that already holds a value
   cannot be used as a pattern name). *)
allSpecies[h_?NumericQ, cCa_?NumericQ, cMg_?NumericQ, cNa_?NumericQ,
           cCl_?NumericQ, pc_?NumericQ] :=
 Module[{co2, hco3, co3, oh},
  (* carbonate fixed by pCO2 and h *)
  co2  = 10^logKH * pc;
  hco3 = 10^logK1 * co2 / h;      (* CO2 + H2O = H+ + HCO3-  *)
  co3  = 10^logK2 * hco3 / h;     (* HCO3-     = H+ + CO3^2- *)
  oh   = 10^logKw / h;            (* H2O       = H+ + OH-    *)
  <|
   "H"      -> h,
   "OH"     -> oh,
   "CO2"    -> co2,
   "HCO3"   -> hco3,
   "CO3"    -> co3,
   (* calcium *)
   "Ca"     -> cCa,
   "CaCl"   -> cCa cCl / 10^logKCaCl,
   "CaCO3"  -> cCa hco3 / (h 10^logKCaCO3),
   "CaHCO3" -> cCa hco3 / 10^logKCaHCO3,
   "CaOH"   -> cCa / (h 10^logKCaOH),
   (* magnesium *)
   "Mg"     -> cMg,
   "MgCl"   -> cMg cCl / 10^logKMgCl,
   "MgCO3"  -> cMg hco3 / (h 10^logKMgCO3),
   "MgHCO3" -> cMg hco3 / 10^logKMgHCO3,
   "MgOH"   -> cMg / (h 10^logKMgOH),
   (* sodium *)
   "Na"     -> cNa,
   "NaCO3"  -> cNa hco3 / (h 10^logKNaCO3),
   "NaHCO3" -> cNa hco3 / 10^logKNaHCO3,
   "NaCl"   -> cNa cCl / 10^logKNaCl,
   "NaOH"   -> cNa / (h 10^logKNaOH),
   (* chloride *)
   "Cl"     -> cCl
  |>
 ];

(* Analytical totals implied by a set of free concentrations *)
totCa[h_, cCa_, cMg_, cNa_, cCl_, p_] := With[{s = allSpecies[h, cCa, cMg, cNa, cCl, p]},
   s["Ca"] + s["CaCl"] + s["CaCO3"] + s["CaHCO3"] + s["CaOH"]];
totMg[h_, cCa_, cMg_, cNa_, cCl_, p_] := With[{s = allSpecies[h, cCa, cMg, cNa, cCl, p]},
   s["Mg"] + s["MgCl"] + s["MgCO3"] + s["MgHCO3"] + s["MgOH"]];
totNa[h_, cCa_, cMg_, cNa_, cCl_, p_] := With[{s = allSpecies[h, cCa, cMg, cNa, cCl, p]},
   s["Na"] + s["NaCO3"] + s["NaHCO3"] + s["NaCl"] + s["NaOH"]];
totCl[h_, cCa_, cMg_, cNa_, cCl_, p_] := With[{s = allSpecies[h, cCa, cMg, cNa, cCl, p]},
   s["Cl"] + s["CaCl"] + s["MgCl"] + s["NaCl"]];

(* -------------------------------------------------------------------- *)
(*  4.  Speciation solver.                                              *)
(*      For a given pH and pCO2, solve the four mass balances for the   *)
(*      free ion concentrations, then return every species.            *)
(* -------------------------------------------------------------------- *)

solveSpeciation[caT_, mgT_, naT_, clT_, pH_, pc_] :=
 Module[{h = 10.^(-pH), cCa, cMg, cNa, cCl, sol},
  sol = FindRoot[
    {
     totCa[h, cCa, cMg, cNa, cCl, pc] == caT,
     totMg[h, cCa, cMg, cNa, cCl, pc] == mgT,
     totNa[h, cCa, cMg, cNa, cCl, pc] == naT,
     totCl[h, cCa, cMg, cNa, cCl, pc] == clT
    },
    {
     {cCa, Max[caT, 1.*^-12]},
     {cMg, Max[mgT, 1.*^-12]},
     {cNa, Max[naT, 1.*^-12]},
     {cCl, Max[clT, 1.*^-12]}
    },
    MaxIterations -> 500];
  allSpecies[h, cCa /. sol, cMg /. sol, cNa /. sol, cCl /. sol, pc]
 ];

(* -------------------------------------------------------------------- *)
(*  5.  Charge balance and surface charge.                             *)
(*                                                                      *)
(*   pos and neg are the reactive-species charge sums requested         *)
(*   (the indifferent electrolyte Na+ / Cl- is left out on purpose:     *)
(*    their imbalance is what the mineral surface carries).             *)
(* -------------------------------------------------------------------- *)

posCharge[s_] := s["H"] + 2 s["Ca"] + 2 s["Mg"] +
   s["CaHCO3"] + s["CaOH"] + s["MgHCO3"] + s["MgOH"];

negCharge[s_] := s["OH"] + s["HCO3"] + 2 s["CO3"] +
   2 s["NaCO3"] + s["NaHCO3"];

sigma[s_] := posCharge[s] - negCharge[s];   (* surface charge, eq/L *)

(* pH from the charge balance: the pH where sigma = 0, at fixed pCO2.       *)
(* A coarse scan locates the first sign change, then bisection refines it.  *)
pHfromCharge[caT_, mgT_, naT_, clT_, pc_] :=
 Module[{f, ps, vals, k, lo, hi, mid},
  f[p_?NumericQ] := sigma[solveSpeciation[caT, mgT, naT, clT, p, pc]];
  ps   = Range[2., 12., 0.25];
  vals = f /@ ps;
  k = SelectFirst[Range[Length[ps] - 1],
        Sign[vals[[#]]] =!= Sign[vals[[# + 1]]] &, $Failed];
  If[k === $Failed, Return[$Failed]];
  lo = ps[[k]]; hi = ps[[k + 1]];
  Do[
    mid = (lo + hi)/2.;
    If[Sign[f[mid]] === Sign[f[lo]], lo = mid, hi = mid],
    {50}];
  (lo + hi)/2.
 ];

(* -------------------------------------------------------------------- *)
(*  6.  Fixed conditions                                                *)
(* -------------------------------------------------------------------- *)

pCO2  = 10^-3.5;   (* atmospheric CO2 partial pressure (atm); change freely *)

(* solid / solution ratios for the surface-charge normalisation          *)
massSolid = 6.0;      (* g of dolomite                          *)
volA = 0.100;         (* L in vessel A                          *)
volC = 0.200;         (* L in vessel C (100 mL A + 100 mL B)    *)

(* -------------------------------------------------------------------- *)
(*  7.  Measured data                                                   *)
(*      columns: pH(meas), Cl_ppm, Na_ppm, Ca_ppm, Mg_ppm               *)
(*      (measured carbonate/alkalinity are kept only for reference;      *)
(*       the carbonate speciation here is set by pCO2, not by them)      *)
(* -------------------------------------------------------------------- *)

(* Vessel A : dolomite + NaCl, equilibrated *)
dataA = {
  {8.8, 39648, 23604.3, 30.7, 42.7},
  {8.8, 32751, 24657.7, 31.0, 46.6},
  {8.8, 36206, 26775.3, 35.0, 49.5},
  {8.9, 32807, 23506.8, 35.7, 45.0},
  {8.9, 48942, 27749.8, 29.4, 43.3},
  {8.9, 42883, 34199.3, 28.3, 41.2},
  {9.0, 38612, 29672.0, 28.2, 41.6},
  {9.0, 33835, 24947.5, 30.1, 43.8},
  {9.0, 30515, 21705.0, 30.6, 45.1}
};

(* Vessel B : NaCl blank, pH pre-adjusted (no Ca, no Mg) *)
(* columns: pH(meas), Cl_ppm, Na_ppm  ( Ca = Mg = 0 )    *)
dataB = {
  {1.8, 16777, 24113},
  {2.1, 17017, 32927},
  {2.1, 16827, 23180},
  {2.4, 16937, 53900},
  {5.5, 16838, 23718},
  {10.5, 16879, 21836},
  {10.6, 16854, 21302},
  {10.6, 16829, 21073},
  {10.8, 17022, 23484}
};

(* Vessel C : mixture A + B after re-equilibration *)
dataC = {
  {7.0, 26800, 16953.4, 1325.3, 233.9},
  {7.2, 36300, 25729.3,  755.1, 180.4},
  {7.5, 32258, 23296.4,  560.6, 145.9},
  {7.5, 24901, 19593.1,  398.9, 114.2},
  {8.8, 24244, 16466.5,   37.2,  24.4},
  {9.8, 23699, 16183.2,   15.6,  13.6},
  {10.2, 21815, 14759.0,  16.3,   4.4},
  {10.2, 24118, 16450.8,  17.0,   3.6},
  {10.6, 22678, 15393.8,   9.8,   1.0}
};

(* -------------------------------------------------------------------- *)
(*  8.  Solve every sample                                              *)
(*      returns: measured pH, calculated pH, full speciation, and the   *)
(*      surface charge (evaluated at the calculated pH and, for         *)
(*      reference, at the measured pH).                                 *)
(* -------------------------------------------------------------------- *)

(* NOTE: the parameter is named pc_ (not pCO2_).  By this point pCO2 already
   holds a numeric value, so writing pCO2_ here would be read as
   Pattern[<number>, _] and Mathematica would reject the whole definition. *)
analyseRow[{pHmeas_, clPpm_, naPpm_, caPpm_, mgPpm_}, pc_] :=
 Module[{caT, mgT, naT, clT, pHc, sCalc, sMeas},
  caT = ppmToM[caPpm, mmCa];
  mgT = ppmToM[mgPpm, mmMg];
  naT = ppmToM[naPpm, mmNa];
  clT = ppmToM[clPpm, mmCl];
  pHc = pHfromCharge[caT, mgT, naT, clT, pc];
  sCalc = solveSpeciation[caT, mgT, naT, clT, pHc, pc];
  sMeas = solveSpeciation[caT, mgT, naT, clT, pHmeas, pc];
  <|
   "pHmeas"  -> pHmeas,
   "pHcalc"  -> pHc,
   "spec"    -> sCalc,           (* speciation at the calculated pH *)
   "Qcalc"   -> sigma[sCalc],    (* = 0 by construction (reference) *)
   "Qmeas"   -> sigma[sMeas]     (* surface charge at the measured pH, eq/L *)
  |>
 ];

(* form for vessel B (no Ca/Mg) reuses the same routine *)
analyseRowB[{pHmeas_, clPpm_, naPpm_}, pc_] :=
  analyseRow[{pHmeas, clPpm, naPpm, 0., 0.}, pc];

resultsA = analyseRow[#, pCO2] & /@ dataA;
resultsB = analyseRowB[#, pCO2] & /@ dataB;
resultsC = analyseRow[#, pCO2] & /@ dataC;

(* -------------------------------------------------------------------- *)
(*  9.  Speciation tables for A, B and C                                *)
(* -------------------------------------------------------------------- *)

speciesOrder = {"H", "OH", "CO2", "HCO3", "CO3",
   "Ca", "CaCl", "CaCO3", "CaHCO3", "CaOH",
   "Mg", "MgCl", "MgCO3", "MgHCO3", "MgOH",
   "Na", "NaCO3", "NaHCO3", "NaCl", "NaOH", "Cl"};

speciationTable[results_, label_] :=
 Module[{header, rows},
  header = Prepend[Table["S" <> ToString[i], {i, Length[results]}], "species (mol/L)"];
  rows = Table[
    Prepend[
     Table[ScientificForm[results[[j, "spec", sp]], 3], {j, Length[results]}],
     sp],
    {sp, speciesOrder}];
  Labeled[
   Grid[Prepend[rows, header], Frame -> All, Alignment -> Left,
    Background -> {None, {{Lighter[Yellow, 0.8], None}}}],
   Style["Vessel " <> label <> " : speciation at the calculated pH", Bold, 14],
   Top]
 ];

pHtable[results_, label_] :=
 Module[{header, rows},
  header = {"sample", "pH measured", "pH calculated", "Q at measured pH (eq/L)"};
  rows = Table[
    {i, results[[i, "pHmeas"]], NumberForm[results[[i, "pHcalc"]], {4, 2}],
     ScientificForm[results[[i, "Qmeas"]], 3]},
    {i, Length[results]}];
  Labeled[
   Grid[Prepend[rows, header], Frame -> All, Alignment -> Left],
   Style["Vessel " <> label <> " : measured vs calculated pH", Bold, 14], Top]
 ];

(* Print the tables *)
Print[pHtable[resultsA, "A"]];
Print[pHtable[resultsB, "B"]];
Print[pHtable[resultsC, "C"]];
Print[speciationTable[resultsA, "A"]];
Print[speciationTable[resultsB, "B"]];
Print[speciationTable[resultsC, "C"]];

(* -------------------------------------------------------------------- *)
(* 10.  Surface charge plots                                            *)
(*                                                                      *)
(*   sigma(pH) is the reactive charge imbalance the surface must carry. *)
(*   It crosses zero at the calculated pH (no surface interaction) and  *)
(*   is negative/positive when the real (measured) pH sits above/below. *)
(* -------------------------------------------------------------------- *)

meanTotals[data_] := Mean[
   (# /. {pHm_, cl_, na_, ca_, mg_} :>
        {ppmToM[ca, mmCa], ppmToM[mg, mmMg], ppmToM[na, mmNa], ppmToM[cl, mmCl]}) & /@
    (If[Length[First[data]] == 3,
       (Join[#, {0., 0.}] &) /@ data,   (* pad B (pH,Cl,Na) with Ca=Mg=0 *)
       data])];

(* p_?NumericQ keeps this unevaluated during Plot's symbolic pass *)
sigmaCurve[{caT_, mgT_, naT_, clT_}, p_?NumericQ] :=
  sigma[solveSpeciation[caT, mgT, naT, clT, p, pCO2]];

(* per-gram conversion (mol charge per g of solid); B has no solid *)
perGram[Q_, vol_] := Q vol / massSolid;

(* --- (a) surface charge titration curve, mean composition of each vessel --- *)
tA = meanTotals[dataA];
tC = meanTotals[dataC];

plotAC = Plot[
   {perGram[sigmaCurve[tA, p], volA]*1000,
    perGram[sigmaCurve[tC, p], volC]*1000},
   {p, 6, 11},
   PlotRange -> All,
   Frame -> True,
   FrameLabel -> {"pH (calculated)", "surface charge  (mmol / g)"},
   PlotLegends -> {"Vessel A (mean)", "Vessel C (mean)"},
   PlotLabel -> "Surface charge from the charge balance  ( pCO2 fixed )",
   GridLines -> Automatic, ImageSize -> 520];

(* --- (b) surface charge of each sample at its measured pH --- *)
ptsA = {#["pHmeas"], perGram[#["Qmeas"], volA]*1000} & /@ resultsA;
ptsC = {#["pHmeas"], perGram[#["Qmeas"], volC]*1000} & /@ resultsC;

plotPts = ListPlot[{ptsA, ptsC},
   Frame -> True,
   FrameLabel -> {"pH (measured)", "surface charge  (mmol / g)"},
   PlotLegends -> {"Vessel A", "Vessel C"},
   PlotStyle -> {Red, Blue},
   PlotMarkers -> {Automatic, 12},
   PlotLabel -> "Surface charge at the measured pH",
   GridLines -> Automatic, ImageSize -> 520];

Print[plotAC];
Print[plotPts];

(* -------------------------------------------------------------------- *)
(*  Notes                                                               *)
(*   1. Concentrations are used directly (activity coefficients = 1).   *)
(*      At the very high ionic strength of these NaCl solutions a       *)
(*      Davies / SIT / Pitzer correction would shift the numbers; add   *)
(*      gamma factors inside allSpecies if needed.                      *)
(*   2. Carbonate is fixed by pCO2. To use the measured alkalinity      *)
(*      instead, replace hco3/co3 in allSpecies by a carbonate mass     *)
(*      balance.                                                        *)
(*   3. The surface charge here is purely a charge-balance quantity --  *)
(*      no surface complexation constants are fitted or assumed.        *)
(* -------------------------------------------------------------------- *)
