(* ::Package:: *)

(* ==================================================================== *)
(*  Dolomite surface titration in NaCl  (after Pokrovsky et al., 1999)  *)
(*                                                                      *)
(*  Net surface charge from the charge-sum difference of Charlet et al. *)
(*  (1990), Pokrovsky Eqn 1:                                            *)
(*     sigmaT = (1/S) ( Sum_k z_k [k]_0  -  Sum_k z_k [k]_f )           *)
(*  [k]_0 : conservative 1:1 mix of reactors A and B (zero-charge ref)  *)
(*  [k]_f : reactor C, computed from the MEASURED pH, alkalinity, Ca,Mg *)
(*  S     : dolomite surface area per litre of reactor C                *)
(*                                                                      *)
(*  No surface complexation model is fitted. Davies activity coeffs and *)
(*  Plummer-Busenberg / Nordstrom formation constants, as in the paper. *)
(*  Solved step by step: gammas -> mass action (gamma substitution) ->  *)
(*  free ions from the mass balances -> charge balance -> surface charge*)
(*                                                                      *)
(*  NOTE: this file mirrors the Python (Colab) notebook, which is the   *)
(*  executed and verified reference.                                    *)
(* ==================================================================== *)

ClearAll["Global`*"];

(* -------------------------------------------------------------------- *)
(*  1.  Constants and thermodynamic data                                *)
(* -------------------------------------------------------------------- *)
ADavies = 0.509;            (* Davies equation constant                 *)
KW      = 1.0*^-14;         (* water dissociation, [H+][OH-]            *)
logKH   = -1.47;            (* CO2(g) = CO2(aq): [CO2aq] = 10^logKH pCO2 *)
pCO2atm = 10^-3.5;          (* atmospheric CO2 partial pressure (atm)   *)

(* formation constants (log10 K), complex on the RIGHT of the reaction  *)
logKa2    = -10.329;   (* HCO3-         = H+ + CO3^2-   *)
logKco2   =  -6.35;    (* CO2(aq) + H2O = H+ + HCO3-    *)
logKCaCO3 =   3.224;   (* Ca2+ + CO3^2- = CaCO3(aq)     *)
logKCaHCO3=   1.106;   (* Ca2+ + HCO3-  = CaHCO3+       *)
logKCaOH  = -12.78;    (* Ca2+ + H2O    = CaOH+ + H+    *)
logKMgCO3 =   2.98;    (* Mg2+ + CO3^2- = MgCO3(aq)     *)
logKMgHCO3=   1.07;    (* Mg2+ + HCO3-  = MgHCO3+       *)
logKMgOH  = -11.44;    (* Mg2+ + H2O    = MgOH+ + H+    *)
logKNaCO3 =   1.27;    (* Na+  + CO3^2- = NaCO3-        *)
logKNaHCO3=  -0.25;    (* Na+  + HCO3-  = NaHCO3(aq)    *)

(* molar masses (g/mol) *)
mmCa = 40.078; mmMg = 24.305; mmCO3 = 60.008; mmHCO3 = 61.016; mmNa = 22.99; mmCl = 35.45;

(* surface area of dolomite per litre of reactor C  (S in Eqn 1) *)
ssa     = 0.84;            (* specific surface area, m2/g              *)
loading = 30.0;            (* 6 g of solid in 0.2 L of reactor C       *)
St      = ssa*loading;     (* = 25.2 m2 / L                            *)

carbonate    = "closed";   (* "closed" = measured alkalinity (Pokrovsky) ; "open" = fixed pCO2 *)
siteCeiling  = 0.0166;     (* site-density ceiling, mmol/m2 (~10 sites/nm2) *)

(* -------------------------------------------------------------------- *)
(*  2.  Step 1 -- Davies activity coefficients                          *)
(*      log gamma_z = -A z^2 ( sqrt I /(1+ sqrt I) - 0.3 I )            *)
(* -------------------------------------------------------------------- *)
gammas[ii_?NumericQ] := If[ii <= 0, {1.0, 1.0},
  Module[{f = -ADavies (Sqrt[ii]/(1 + Sqrt[ii]) - 0.3 ii)},
   {10^f, 10^(4 f)}]];       (* {gamma1 (z=1), gamma2 (z=2)} *)

(* -------------------------------------------------------------------- *)
(*  3.  Steps 2-3 -- mass action (gamma substituted) and free ions      *)
(*      For a MEASURED pH, build every species from the free carbonate  *)
(*      cCO3 and the mass-balance denominators, iterating on I.         *)
(*      "open"  : CO2(aq) fixed by pCO2, cCO3 follows from pH.           *)
(*      "closed": cCO3 scaled until total inorganic carbon = DIC.        *)
(* -------------------------------------------------------------------- *)
buildFromCO3[cCO3_, aH_, g1_, g2_, CaT_, MgT_, NaT_, ClT_] :=
 Module[{aCO3, aHCO3, cHCO3, cCO2, dCa, dMg, dNa, cCa, cMg, cNa,
         CaCO3, CaHCO3, MgCO3, MgHCO3, NaCO3, NaHCO3, CaOH, MgOH, DICc, cH},
  cH    = aH/g1;
  aCO3  = g2 cCO3;
  aHCO3 = aCO3 aH/10^logKa2;
  cHCO3 = aHCO3/g1;
  cCO2  = aHCO3 aH/10^logKco2;
  (* free metals from the mass-balance denominators *)
  dCa = 1 + 10^logKCaCO3 g2 g2 cCO3 + 10^logKCaHCO3 g2 cHCO3 + 10^logKCaOH g2/(aH g1);
  dMg = 1 + 10^logKMgCO3 g2 g2 cCO3 + 10^logKMgHCO3 g2 cHCO3 + 10^logKMgOH g2/(aH g1);
  dNa = 1 + 10^logKNaCO3 g2 cCO3 + 10^logKNaHCO3 g1 g1 cHCO3;
  cCa = CaT/dCa; cMg = MgT/dMg; cNa = NaT/dNa;
  (* rebuild every complex from the free ions *)
  CaCO3  = 10^logKCaCO3  g2 g2 cCa cCO3;
  CaHCO3 = 10^logKCaHCO3 g2 cCa cHCO3;
  MgCO3  = 10^logKMgCO3  g2 g2 cMg cCO3;
  MgHCO3 = 10^logKMgHCO3 g2 cMg cHCO3;
  NaCO3  = 10^logKNaCO3  g2 cNa cCO3;
  NaHCO3 = 10^logKNaHCO3 g1 g1 cNa cHCO3;
  CaOH   = 10^logKCaOH g2 cCa/(aH g1);
  MgOH   = 10^logKMgOH g2 cMg/(aH g1);
  DICc   = cCO3 + cHCO3 + cCO2 + CaCO3 + CaHCO3 + MgCO3 + MgHCO3 + NaCO3 + NaHCO3;
  <|"H" -> cH, "HCO3" -> cHCO3, "CO2" -> cCO2, "Ca" -> cCa, "Mg" -> cMg, "Na" -> cNa,
    "CaCO3" -> CaCO3, "CaHCO3" -> CaHCO3, "MgCO3" -> MgCO3, "MgHCO3" -> MgHCO3,
    "NaCO3" -> NaCO3, "NaHCO3" -> NaHCO3, "CaOH" -> CaOH, "MgOH" -> MgOH, "DICc" -> DICc|>];

speciate[pH_, CaT_, MgT_, NaT_, ClT_, DIC_: Automatic] :=
 Module[{aH = 10^-pH, ii, g1, g2, cCO3, s, cOH, iiNew, k, aHCO3},
  ii = 0.5 (NaT + ClT + 4 CaT + 4 MgT);
  cCO3 = Max[If[DIC === Automatic, 1.*^-3, DIC] 0.1, 1.*^-12];
  Do[
    {g1, g2} = gammas[ii];
    If[carbonate === "open",
      Module[{cCO2 = 10^logKH pCO2atm},
        aHCO3 = 10^logKco2 cCO2/aH;
        cCO3 = (10^logKa2 aHCO3/aH)/g2];
      s = buildFromCO3[cCO3, aH, g1, g2, CaT, MgT, NaT, ClT],
      (* closed: inner loop scales cCO3 to the measured DIC *)
      Do[ s = buildFromCO3[cCO3, aH, g1, g2, CaT, MgT, NaT, ClT];
          Module[{new = cCO3 (DIC/s["DICc"])},
            If[Abs[new - cCO3] < 1.*^-13 Max[cCO3, 1.*^-15], cCO3 = new; Break[]];
            cCO3 = new], {kk, 300}];
      s = buildFromCO3[cCO3, aH, g1, g2, CaT, MgT, NaT, ClT]];
    cOH = KW/aH/g1;
    iiNew = 0.5 (s["Na"] + ClT + s["H"] + cOH + 4 s["Ca"] + 4 s["Mg"] + 4 cCO3
                 + s["HCO3"] + s["CaHCO3"] + s["CaOH"] + s["MgHCO3"] + s["MgOH"] + s["NaCO3"]);
    If[Abs[iiNew - ii] < 1.*^-11, ii = iiNew; Break[]];
    ii = 0.5 ii + 0.5 iiNew,
    {k, 300}];
  <|"pH" -> pH, "I" -> ii, "H" -> s["H"], "OH" -> cOH,
    "Ca" -> s["Ca"], "Mg" -> s["Mg"], "Na" -> s["Na"], "Cl" -> ClT,
    "CO3" -> cCO3, "HCO3" -> s["HCO3"], "CO2" -> s["CO2"],
    "CaHCO3" -> s["CaHCO3"], "CaOH" -> s["CaOH"], "MgHCO3" -> s["MgHCO3"],
    "MgOH" -> s["MgOH"], "NaCO3" -> s["NaCO3"], "DIC" -> s["DICc"]|>];

(* -------------------------------------------------------------------- *)
(*  4.  Step 4 -- charge balance  Sum_k z_k [k]  (reactive)             *)
(*      Na+ and Cl- are omitted: conserved on mixing, they cancel in    *)
(*      the [k]_0 - [k]_f difference.                                   *)
(* -------------------------------------------------------------------- *)
qCharge[s_] :=
  (+1 s["H"] - 1 s["OH"])                              (* q_H   *)
  + (-1 s["HCO3"] - 2 s["CO3"] - 1 s["NaCO3"])         (* q_DIC *)
  + (+2 s["Ca"] + 1 s["CaHCO3"] + 1 s["CaOH"])         (* q_Ca  *)
  + (+2 s["Mg"] + 1 s["MgHCO3"] + 1 s["MgOH"]);        (* q_Mg  *)

(* -------------------------------------------------------------------- *)
(*  5.  Measured data (run-aligned): {pH, Cl, Na, Ca, Mg, CO3, HCO3} mg/L*)
(* -------------------------------------------------------------------- *)
dataA = {{8.9,48942,27749.8,29.4,43.3,27.9,119.8},{8.8,39648,23604.3,30.7,42.7,23.4,127.8},
  {8.9,42883,34199.3,28.3,41.2,32.1,107.8},{9.0,38612,29672.0,28.2,41.6,31.7,116.1},
  {9.0,33835,24947.5,30.1,43.8,30.4,109.7},{8.8,36206,26775.3,35.0,49.5,23.5,125.5},
  {8.8,32751,24657.7,31.0,46.6,20.8,131.5},{9.0,30515,21705.0,30.6,45.1,35.0,121.0},
  {8.9,32807,23506.8,35.7,45.0,22.9,128.3}};
dataB = {{2.4,16937,53900,0,0,0,0},{2.1,17017,32927,0,0,0,0},{2.1,16827,23180,0,0,0,0},
  {1.8,16777,24113,0,0,0,0},{5.5,16838,23718,0,0,0,0},{10.5,16879,21836,0,0,48.0,0},
  {10.6,16854,21302,0,0,46.6,0},{10.6,16829,21073,0,0,55.0,0},{10.8,17022,23484,0,0,80.1,0}};
dataC = {{7.5,24901,19593.1,398.9,114.2,0,407.6},{7.5,32258,23296.4,560.6,145.9,0,413.3},
  {7.2,36300,25729.3,755.1,180.4,0,250.1},{7.0,26800,16953.4,1325.3,233.9,0,527.5},
  {8.8,24244,16466.5,37.2,24.4,18.2,106.1},{9.8,23699,16183.2,15.6,13.6,92.8,7.2},
  {10.2,21815,14759.0,16.3,4.4,120.2,0},{10.2,24118,16450.8,17.0,3.6,134.0,0},
  {10.6,22678,15393.8,9.8,1.0,136.4,0}};

(* strong titrant added to each B run (mol/L): runs 1-4 HCl, runs 6-9 NaOH *)
HClB  = {0.001,0.002,0.003,0.004,0,0,0,0,0};
NaOHB = {0,0,0,0,0,0.001,0.002,0.003,0.004};

totals[row_] := Module[{},
  {row[[4]]/mmCa/1000., row[[5]]/mmMg/1000., row[[3]]/mmNa/1000., row[[2]]/mmCl/1000.,
   (row[[6]]/mmCO3 + row[[7]]/mmHCO3)/1000.}];

(* -------------------------------------------------------------------- *)
(*  6.  Step 5 -- surface charge (Pokrovsky Eqn 1)                       *)
(* -------------------------------------------------------------------- *)
speciateRow[row_] := Module[{t = totals[row]},
  speciate[row[[1]], t[[1]], t[[2]], t[[3]], t[[4]], t[[5]]]];

surfaceCharge := Table[
  Module[{sA = speciateRow[dataA[[i]]], sB = speciateRow[dataB[[i]]],
          sC = speciateRow[dataC[[i]]], q0, qf},
   q0 = 0.5 (qCharge[sA] + qCharge[sB]);   (* Sum z_k [k]_0 *)
   qf = qCharge[sC];                       (* Sum z_k [k]_f *)
   <|"run" -> i, "pHC" -> dataC[[i, 1]], "sigmaT" -> (q0 - qf)/St*1000.|>], {i, 9}];

(* -------------------------------------------------------------------- *)
(*  7.  Data-quality budgets  (runs 1-4 fail the acid balance)          *)
(* -------------------------------------------------------------------- *)
budgets := Table[
  Module[{tA = totals[dataA[[i]]], tC = totals[dataC[[i]]], acid, dCat},
   dCat = 2 ((tC[[1]] + tC[[2]]) - 0.5 (tA[[1]] + tA[[2]]));
   acid = (HClB[[i]] - NaOHB[[i]])/2;
   <|"run" -> i, "pHC" -> dataC[[i, 1]], "acidMeq" -> acid*1000.,
     "dCatMeq" -> dCat*1000.,
     "ratio" -> If[acid != 0, Abs[dCat/acid], Missing[]],
     "coherent" -> If[i >= 5, "yes", "no (acid budget fails)"]|>], {i, 9}];

(* -------------------------------------------------------------------- *)
(*  8.  Titration curve -- coherent runs 5-9 vs the MEASURED pH          *)
(* -------------------------------------------------------------------- *)
plotCurve := Module[{sc = surfaceCharge, coh, drp},
  coh = Select[sc, #["run"] >= 5 &];
  drp = Select[sc, #["run"] <= 4 &];
  Show[
   ListLinePlot[{#["pHC"], #["sigmaT"]} & /@ coh, PlotMarkers -> Automatic,
     PlotStyle -> RGBColor[0.17, 0.37, 0.54], PlotLegends -> {"runs 5-9 (coherent)"}],
   ListPlot[{#["pHC"], #["sigmaT"]} & /@ drp, PlotStyle -> Gray],
   Plot[{siteCeiling, -siteCeiling}, {x, 6.8, 10.8},
     PlotStyle -> Directive[Red, Dotted]],
   Frame -> True, FrameLabel -> {"measured pH of reactor C", "sigmaT (mmol/m2)"},
   PlotLabel -> "Dolomite surface charge -- Pokrovsky Eqn 1 (" <> carbonate <> ")"]];

(* Evaluate:  Dataset[surfaceCharge] , Dataset[budgets] , plotCurve *)
