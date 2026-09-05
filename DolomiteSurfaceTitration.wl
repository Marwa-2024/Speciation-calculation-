(* ::Package:: *)

(* ==================================================================== *)
(*  Dolomite in NaCl -- Visual MINTEQ-style speciation                  *)
(*  (fixed pCO2, pH calculated from the mass balance)                   *)
(*                                                                      *)
(*  Reproduces the Visual MINTEQ speciation for the dolomite titration. *)
(*   - Carbonate is OPEN: {CO2(aq)} is fixed (fixed pCO2).              *)
(*   - pH is CALCULATED from the mass balance (MINTEQ option 1), with   *)
(*       the free H+ started from a neutral guess 1E-7.                 *)
(*   - Davies activity coefficients, TOUGHREACT / EQ3-6 constants.      *)
(*   - No surface complexation model.                                   *)
(*                                                                      *)
(*  ALL equations are written out below: the reactions, the law of      *)
(*  mass action for every species, the element mass balances, and the   *)
(*  proton condition that fixes the pH.  This mirrors the Python (Colab) *)
(*  notebook, which is the executed and verified reference.             *)
(* ==================================================================== *)

ClearAll["Global`*"];

(* -------------------------------------------------------------------- *)
(*  1.  REACTIONS  (dissociation; complex on the left)  and  log10 K    *)
(* -------------------------------------------------------------------- *)
(*  log K are Visual MINTEQ's own values, recovered from the MINTEQ output. *)
(*  Water / carbonate                                                    *)
(*    H2O            = H+ + OH-                    logK = -14.011         *)
(*    CO2(aq) + H2O  = H+ + HCO3-                  logK =  -6.346         *)
(*    HCO3-          = H+ + CO3^2-                 logK = -10.335         *)
(*    CO2(g)         = CO2(aq)  -> {CO2(aq)} fixed by pCO2                *)
(*  Complexes                                                            *)
(*    CaCl+          = Ca2+ + Cl-                  logK =  -0.400         *)
(*    CaCO3(aq) + H+ = Ca2+ + HCO3-                logK =   6.971         *)
(*    CaHCO3+        = Ca2+ + HCO3-                logK =  -1.091         *)
(*    CaOH+  + H+    = Ca2+ + H2O                  logK =  12.711         *)
(*    MgCl+          = Mg2+ + Cl-                  logK =  -0.600         *)
(*    MgCO3(aq) + H+ = Mg2+ + HCO3-                logK =   7.355         *)
(*    MgHCO3+        = Mg2+ + HCO3-                logK =  -1.060         *)
(*    MgOH+  + H+    = Mg2+ + H2O                  logK =  11.798         *)
(*    NaCl(aq)       = Na+ + Cl-                   logK =   0.300         *)
(*    NaCO3- + H+    = Na+ + HCO3-                 logK =   9.065         *)
(*    NaHCO3(aq)     = Na+ + HCO3-                 logK =   0.306         *)
(*    NaOH(aq) + H+  = Na+ + H2O                   logK =  13.911         *)
(* -------------------------------------------------------------------- *)
(* log K recovered from the Visual MINTEQ output itself (its own database) *)
lkco2=-6.346; lkco3=10.335; lkoh=14.011;
lkcacl=-0.400; lkcaco3=6.971; lkcahco3=-1.091; lkcaoh=12.711;
lkmgcl=-0.600; lkmgco3=7.355; lkmghco3=-1.060; lkmgoh=11.798;
lknacl=0.300; lknaco3=9.065; lknahco3=0.306; lknaoh=13.911;

aCO2 = 10^-4.90;          (* fixed {CO2(aq)} activity == fixed pCO2 (MINTEQ log activity = -4.9) *)
ADavies = 0.509;
daviesB = 0.5;            (* Davies b-term recovered from MINTEQ (its default here, not 0.3) *)
mmCa=40.078; mmMg=24.305; mmNa=22.99; mmCl=35.45; mmCO3=60.008; mmHCO3=61.016;
St = 0.84*30.0;          (* dolomite area per litre of reactor C = 25.2 m2/L *)

(* -------------------------------------------------------------------- *)
(*  2.  Activity coefficients (Davies)                                  *)
(*      log g_z = -A z^2 ( sqrt I/(1+sqrt I) - 0.3 I ) ;  g0 = 10^(0.1 I)*)
(* -------------------------------------------------------------------- *)
gammas[ii_?NumericQ] := If[ii <= 0, {1.,1.,1.},
  Module[{f = -ADavies (Sqrt[ii]/(1+Sqrt[ii]) - daviesB ii)}, {10^f, 10^(4 f), 10^(0.1 ii)}]];

(* -------------------------------------------------------------------- *)
(*  3.  LAW OF MASS ACTION  and  ELEMENT MASS BALANCES                  *)
(*                                                                      *)
(*  Carbonate activities (aH = {H+}, aCO2 fixed):                       *)
(*     {HCO3-} = 10^lkco2 aCO2 / aH                                     *)
(*     {CO3-2} = {HCO3-} / (10^lkco3 aH)                                *)
(*     {OH-}   = 10^-lkoh / aH                                          *)
(*                                                                      *)
(*  Complex activities, e.g.:                                           *)
(*     {CaHCO3+} = 10^-lkcahco3 {Ca2+}{HCO3-}                           *)
(*     {CaCl+}   = {Ca2+}{Cl-} / 10^lkcacl                              *)
(*     {CaOH+}   = {Ca2+} / (10^lkcaoh aH)                              *)
(*                                                                      *)
(*  Mass balances (free ion = total - bound; no surface terms):         *)
(*     CaT = [Ca2+]+[CaCl+]+[CaCO3]+[CaHCO3+]+[CaOH+]                    *)
(*     MgT = [Mg2+]+[MgCl+]+[MgCO3]+[MgHCO3+]+[MgOH+]                    *)
(*     NaT = [Na+] +[NaCl] +[NaCO3-]+[NaHCO3]+[NaOH]                     *)
(*     ClT = [Cl-] +[CaCl+]+[MgCl+]+[NaCl]                              *)
(*  Each metal total gives  aMe = MeT/dMe ; aCl follows from ClT.       *)
(*  aCa,aMg,aNa,aCl and the ionic strength are iterated to convergence. *)
(* -------------------------------------------------------------------- *)
speciate[pH_, CaT_, MgT_, NaT_, ClT_] :=
 Module[{aH = 10^-pH, ii, g1, g2, g0, aHCO3, aCO3, aOH, aCa, aMg, aNa, aCl,
         dCa, dMg, dNa, dCl, aClNew, k, kk,
         cH, cOH, cHCO3, cCO3, cCO2, cCa, cMg, cNa, cCl,
         cCaCl, cCaCO3, cCaHCO3, cCaOH, cMgCl, cMgCO3, cMgHCO3, cMgOH,
         cNaCl, cNaCO3, cNaHCO3, cNaOH, iiNew},
  ii = 0.5 (NaT + ClT + 4 CaT + 4 MgT);
  aCa = aMg = aNa = aCl = 0.;
  Do[
    {g1, g2, g0} = gammas[ii];
    aHCO3 = 10^lkco2 aCO2/aH;
    aCO3  = aHCO3/(10^lkco3 aH);
    aOH   = 10^-lkoh/aH;
    cH = aH/g1; cOH = aOH/g1; cHCO3 = aHCO3/g1; cCO3 = aCO3/g2; cCO2 = aCO2/g0;
    If[aCl == 0., aCa = g2 CaT; aMg = g2 MgT; aNa = g1 NaT; aCl = g1 ClT];
    Do[
      dCa = 1/g2 + aCl/(g1 10^lkcacl) + aHCO3/(g0 aH 10^lkcaco3)
            + aHCO3 10^(-lkcahco3)/g1 + 1/(g1 aH 10^lkcaoh);
      dMg = 1/g2 + aCl/(g1 10^lkmgcl) + aHCO3/(g0 aH 10^lkmgco3)
            + aHCO3 10^(-lkmghco3)/g1 + 1/(g1 aH 10^lkmgoh);
      dNa = 1/g1 + aHCO3/(g1 aH 10^lknaco3) + aHCO3 10^(-lknahco3)/g0
            + aCl/(g0 10^lknacl) + 1/(g0 aH 10^lknaoh);
      aCa = CaT/dCa; aMg = MgT/dMg; aNa = NaT/dNa;
      dCl = 1/g1 + aCa/(g1 10^lkcacl) + aMg/(g1 10^lkmgcl) + aNa/(g0 10^lknacl);
      aClNew = ClT/dCl;
      If[Abs[aClNew - aCl] < 1.*^-12, aCl = aClNew; Break[]];
      aCl = aClNew, {kk, 100}];
    cCaCl = aCa aCl/(g1 10^lkcacl); cCaCO3 = aCa aHCO3/(g0 aH 10^lkcaco3);
    cCaHCO3 = aCa aHCO3 10^(-lkcahco3)/g1; cCaOH = aCa/(g1 aH 10^lkcaoh);
    cMgCl = aMg aCl/(g1 10^lkmgcl); cMgCO3 = aMg aHCO3/(g0 aH 10^lkmgco3);
    cMgHCO3 = aMg aHCO3 10^(-lkmghco3)/g1; cMgOH = aMg/(g1 aH 10^lkmgoh);
    cNaCl = aNa aCl/(g0 10^lknacl); cNaCO3 = aNa aHCO3/(g1 aH 10^lknaco3);
    cNaHCO3 = aNa aHCO3 10^(-lknahco3)/g0; cNaOH = aNa/(g0 aH 10^lknaoh);
    cCa = aCa/g2; cMg = aMg/g2; cNa = aNa/g1; cCl = aCl/g1;
    iiNew = 0.5 (cNa + cCl + cH + cOH + 4 cCa + 4 cMg + 4 cCO3 + cHCO3
                 + cCaCl + cCaHCO3 + cCaOH + cMgCl + cMgHCO3 + cMgOH + cNaCO3);
    If[Abs[iiNew - ii] < 1.*^-10, ii = iiNew; Break[]];
    ii = 0.5 ii + 0.5 iiNew, {k, 200}];
  <|"pH" -> pH, "I" -> ii, "H" -> cH, "OH" -> cOH, "CO2" -> cCO2, "HCO3" -> cHCO3, "CO3" -> cCO3,
    "Ca" -> cCa, "CaCl" -> cCaCl, "CaCO3" -> cCaCO3, "CaHCO3" -> cCaHCO3, "CaOH" -> cCaOH,
    "Mg" -> cMg, "MgCl" -> cMgCl, "MgCO3" -> cMgCO3, "MgHCO3" -> cMgHCO3, "MgOH" -> cMgOH,
    "Na" -> cNa, "NaCl" -> cNaCl, "NaCO3" -> cNaCO3, "NaHCO3" -> cNaHCO3, "NaOH" -> cNaOH, "Cl" -> cCl|>];

(* -------------------------------------------------------------------- *)
(*  4.  pH FROM THE MASS BALANCE  (proton condition, fixed pCO2)        *)
(*                                                                      *)
(*  P(pH) = [H+] - [OH-] - [HCO3-] - 2[CO3-2]                           *)
(*          - [NaHCO3] - 2[NaCO3-] - [CaHCO3+] - 2[CaCO3] - [CaOH+]     *)
(*          - [MgHCO3+] - 2[MgCO3] - [MgOH+] - [NaOH]                   *)
(*  MINTEQ sets total H+ from the entered carbonate, so pH solves       *)
(*        P(pH) = -2 Ccarb    (Ccarb = measured total CO3-2).           *)
(*  Ca2+, Mg2+, Na+, Cl- carry no proton, so dissolved Ca/Mg do not     *)
(*  drag the pH the wrong way.  H+ is started from the neutral 1E-7.     *)
(* -------------------------------------------------------------------- *)
protonCondition[s_] :=
  s["H"] - s["OH"] - s["HCO3"] - 2 s["CO3"]
  - s["NaHCO3"] - 2 s["NaCO3"] - s["CaHCO3"] - 2 s["CaCO3"] - s["CaOH"]
  - s["MgHCO3"] - 2 s["MgCO3"] - s["MgOH"] - s["NaOH"];

solvePH[CaT_, MgT_, NaT_, ClT_, Ccarb_] :=
 Module[{f}, f[pH_?NumericQ] := protonCondition[speciate[pH, CaT, MgT, NaT, ClT]] + 2 Ccarb;
  pH /. FindRoot[f[pH] == 0, {pH, 7., 3., 11.5}]];

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

totals[row_] := {row[[4]]/mmCa/1000., row[[5]]/mmMg/1000., row[[3]]/mmNa/1000.,
                 row[[2]]/mmCl/1000., row[[6]]/mmCO3/1000.};   (* Ca,Mg,Na,Cl, Ccarb *)

speciateRow[row_] := Module[{t = totals[row], pH},
  pH = solvePH[t[[1]], t[[2]], t[[3]], t[[4]], t[[5]]];
  Append[speciate[pH, t[[1]], t[[2]], t[[3]], t[[4]]], "pHmeas" -> row[[1]]]];

(* -------------------------------------------------------------------- *)
(*  6.  Net surface charge -- Pokrovsky Eqn 1 (at the calculated pH)     *)
(*      sigmaT = ( (qA+qB)/2 - qC ) / S ,  q = Sum z_k [k] (Na,Cl cancel)*)
(* -------------------------------------------------------------------- *)
qCharge[s_] :=
  (s["H"] - s["OH"]) + (-s["HCO3"] - 2 s["CO3"] - s["NaCO3"])
  + (2 s["Ca"] + s["CaCl"] + s["CaHCO3"] + s["CaOH"])
  + (2 s["Mg"] + s["MgCl"] + s["MgHCO3"] + s["MgOH"]);

surfaceCharge := Table[
  Module[{sA = speciateRow[dataA[[i]]], sB = speciateRow[dataB[[i]]], sC = speciateRow[dataC[[i]]]},
   <|"run" -> i, "pHC" -> sC["pH"],
     "sigmaT" -> (0.5 (qCharge[sA] + qCharge[sB]) - qCharge[sC])/St*1000.|>], {i, 9}];

(* Evaluate:  Dataset[Table[speciateRow[dataC[[i]]], {i,9}]]  ;  Dataset[surfaceCharge] *)
