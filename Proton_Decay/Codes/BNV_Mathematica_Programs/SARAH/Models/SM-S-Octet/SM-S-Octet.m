Off[General::spell]

Model`Name = "SM_S_Oc";
Model`NameLaTeX ="Standard Model extended by real singlet and color octet";
Model`Authors = "F.Staub, M. Goodsell";
Model`Date = "2016-02-11";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False};
Gauge[[2]]={WB, SU[2], left,        g2,True};
Gauge[[3]]={G,  SU[3], color,       g3,False};


(* Chiral Superfields *)

FermionFields[[1]] = {q, 3, {uL, dL},   1/6, 2, 3};  
FermionFields[[2]] = {l, 3, {vL, eL},  -1/2, 2, 1};
FermionFields[[3]] = {d, 3, conj[dR],   1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],  -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],     1, 1,  1};

ScalarFields[[1]] =  {H,  1, {Hp, H0},    1/2, 2, 1};
ScalarFields[[2]] =  {s,  1, Sing,          0, 1, 1};
ScalarFields[[3]] =  {oc, 1, {OcP,Oc0},    1/2, 2, 8};

RealScalars = {s};


(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][Additional]= {
	{LagHC, {Overwrite->False, AddHC->True}},
	{LagNoHC,{Overwrite->False, AddHC->False}}
};


(* The O^3 and O^4 terms play hardly a role *)
(* to speed up the code, let's put them to zero *)
Lambda4=0;
Lambda5=0;
Lambda6=0;
Lambda7=0;
Lambda8=0;
Lambda9=0;
Lambda10=0;
Lambda11=0;


LagNoHC = -(mu2 conj[H].H + MS2/2 s.s + LambdaS s.s.s.s + LambdaH conj[H].H.conj[H].H + Kappa1  conj[H].H.s.s + \
 MO2 conj[oc].oc + Kappa2 s.s.conj[oc].oc + Lambda1 conj[H].H.conj[oc].oc + Lambda2 Delta[col3,col4] Delta[lef1,lef4] Delta[lef2,lef3] conj[H].H.conj[oc].oc + \
 Lambda6/16 Delta[lef1,lef2] Delta[lef3,lef4] sum[IndC1,1,3] sum[IndC2,1,3] sum[IndC3,1,3] sum[IndC4,1,3] conj[Lam[col1,IndC1,IndC2]] Lam[col2,IndC2,IndC3] conj[Lam[col3,IndC3,IndC4]] Lam[col4,IndC4,IndC1] conj[oc].oc.conj[oc].oc + \
 Lambda7/16 Delta[lef1,lef4] Delta[lef3,lef2] sum[IndC1,1,3] sum[IndC2,1,3] sum[IndC3,1,3] sum[IndC4,1,3] conj[Lam[col1,IndC1,IndC2]] Lam[col2,IndC2,IndC3] conj[Lam[col3,IndC3,IndC4]] Lam[col4,IndC4,IndC1] conj[oc].oc.conj[oc].oc + \
 Lambda8 conj[oc].oc.conj[oc].oc + Lambda9 Delta[lef1,lef4] Delta[lef2,lef3] Delta[col1,col2] Delta[col3,col4]  conj[oc].oc.conj[oc].oc + \
 Lambda9 Delta[lef1,lef3] Delta[lef2,lef4] Delta[col1,col2] Delta[col3,col4]  oc.oc.conj[oc].conj[oc] + \
 Lambda11/16 Delta[lef1,lef3] Delta[lef2,lef4] sum[IndC1,1,3] sum[IndC2,1,3] sum[IndC3,1,3] sum[IndC4,1,3]   Lam[col1,IndC1,IndC2] Lam[col2,IndC2,IndC3] conj[Lam[col3,IndC3,IndC4]] conj[Lam[col4,IndC4,IndC1]]  oc.oc.conj[oc].conj[oc]);


LagHC = - ( Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q + 
 Lambda3 Delta[lef1,lef3] Delta[lef2,lef4] Delta[col3,col4] conj[H].conj[H].oc.oc + 
 Lambda4/8 Delta[left1,lef4] Delta[lef2,lef3] sum[IndC1,1,3] sum[IndC2,1,3] sum[IndC3,1,3] conj[Lam[col2,IndC1,IndC2]] Lam[col2,IndC2,IndC3] Lam[col3,IndC3,IndC1] conj[H].conj[oc].oc.oc + \
 Lambda5/8 Delta[left1,lef3] Delta[lef2,lef4] sum[IndC1,1,3] sum[IndC2,1,3] sum[IndC3,1,3] conj[Lam[col2,IndC1,IndC2]] Lam[col2,IndC2,IndC3] Lam[col3,IndC3,IndC1] conj[H].conj[oc].oc.oc  );


(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};     
        
        
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
{     {H0, {v, 1/Sqrt[2]}, {Ah, \[ImaginaryI]/Sqrt[2]},{phiH, 1/Sqrt[2]}},
      {Sing, {vS, 1}, {0, 0},{phiS, 1}},
      {Oc0, {0,0}, {OcI, \[ImaginaryI]/Sqrt[2]},{OcR, 1/Sqrt[2]} }};
 

DEFINITION[EWSB][MatterSector]=   
    {{{phiH,phiS},{hh,ZH}},
     {{{dL}, {conj[dR]}}, {{FDL,Vd}, {FDR,Ud}}},
     {{{uL}, {conj[uR]}}, {{FUL,Vu}, {FUR,Uu}}},
     {{{eL}, {conj[eR]}}, {{FEL,Ve}, {FER,Ue}}}};  


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  FDL, conj[FDR]},
 Fe ->{  FEL, conj[FER]},
 Fu ->{  FUL, conj[FUR]},
 Fv ->{  vL, 0}};

DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};


(* Model not supported by CalcHep *)
SetOptions[MakeAll, IncludeCalcHep->False];

