(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEFusaoka`",{"REAP`RGESymbol`"}];

ClearAll[RGEFusaokaYukawa];

RGEFusaokaYukawa::usage="RGEFusaokaYukawa[Model,Particle] returns the Yukawa matrix for the particle type <Particle> in the <Model> at the scale <Scale>. The possible options for <Scale> are 'GUT' and 'MZ'";



 Begin["`Private`"];

 Map[Needs,{"REAP`RGESymbol`"}];
(* data of paper arxiv:hep-ph/9712201 v2 from Fusaoka and Koide *)
(* they used
MassHiggs=246.2;
tan\[Beta]=10;
(* GUT: 2*10^16 *)
The matrices are given at the scale MX
*)

Begin["`SM`"];
Begin["`GUT`"];
ClearAll[Dagger];
Dagger[px_]:=Transpose[Conjugate[px]];
ClearAll[UnitaryMatrix];
UnitaryMatrix[pMatrix_]:=Dagger[pMatrix]+pMatrix;

(* SM *)
MassTop=84.2;
MassBottom=1.071;

(* quark mass matrices *)
ClearAll[Mu];
Mu[pMassTop_]:=pMassTop*DiagonalMatrix[{1.11*10^-5,3.23*10^-3,1}];
ClearAll[Md];
Md[pMassBottom_]:=pMassBottom*UnitaryMatrix[{{0.0035/2,0.0074*Exp[-1.2\[ImaginaryI]],0.0035*Exp[-95.3\[ImaginaryI]]},{0,0.0363/2,0.0418*Exp[0.03\[ImaginaryI]]},{0,0,0.9982/2}}];
(* Mass["s"]<0 *)
ClearAll[Md2];
Md2[pMassBottom_]:=pMassBottom*UnitaryMatrix[{{-1.9*10^(-5)/2,-0.0082*Exp[1.1\[ImaginaryI]],0.0035*Exp[-84.1\[ImaginaryI]]},{0,-0.0324/2,0.0447*Exp[-0.04\[ImaginaryI]]},{0,0,0.9980/2}}];

(* CKM matrix *)
VCKM={{0.9754,0.2206,-0.0035\[ImaginaryI]},{-0.2203,0.9745,0.0433},{0.0101*Exp[-19\[ImaginaryI]],-0.0422*Exp[1.0\[ImaginaryI]],0.9991}};


(* Yukawa coupling matrices *)
Yu=Sqrt[2]/RGEvEW*Mu[MassTop];
Yd=Sqrt[2]/RGEvEW*Md[MassBottom];
Yd2=Sqrt[2]/RGEvEW*Md2[MassBottom];
End[];


Begin["`MZ`"];
ClearAll[Dagger];
Dagger[px_]:=Transpose[Conjugate[px]];
ClearAll[UnitaryMatrix];
UnitaryMatrix[pMatrix_]:=Dagger[pMatrix]+pMatrix;

(* SM *)
MassTop=180;
MassBottom=3.00;

(* quark mass matrices *)
ClearAll[Mu];
Mu[pMassTop_]:=pMassTop*DiagonalMatrix[{1.29*10^-5,3.75*10^-3,1}];
ClearAll[Md];
(* Mass["s"]>0 *)
Md[pMassBottom_,\[Delta]_]:=pMassBottom*UnitaryMatrix[{{3.01*10^-3/2,(6.36+0.11*Exp[-\[ImaginaryI] \[Delta]])*10^-3,(-0.24+2.97*Exp[-\[ImaginaryI] \[Delta]])*10^-3},{0,0.0310/2,0.0362},{0,0,0.9986/2}}];


(* Yukawa coupling matrices *)
Yu=Sqrt[2]/RGEvEW*Mu[MassTop];
Yd=Sqrt[2]/RGEvEW*Md[MassBottom,0];
End[];
End[];


Begin["`MSSM`"];
Begin["`GUT`"];
ClearAll[Dagger];
Dagger[px_]:=Transpose[Conjugate[px]];
ClearAll[UnitaryMatrix];
UnitaryMatrix[pMatrix_]:=Dagger[pMatrix]+pMatrix;

(* MSSM *)
MassTop=129.3;
MassBottom=0.997;

(* quark mass matrices *)
ClearAll[Mu];
Mu[pMassTop_]:=pMassTop*DiagonalMatrix[{8.0*10^-6,2.33*10^-3,1}];
ClearAll[Md];
Md[pMassBottom_]:=pMassBottom*UnitaryMatrix[{{0.0026/2,0.0054*Exp[-0.9\[ImaginaryI]],0.0025*Exp[-93.9\[ImaginaryI]]},{0,0.0263/2,0.0310*Exp[0.03\[ImaginaryI]]},{0,0,0.9990/2}}];
(* Mass["s"]<0 *)
ClearAll[Md2];
Md2[pMassBottom_]:=pMassBottom*UnitaryMatrix[{{-1.6*10^(-5)/2,-0.0060*Exp[0.8\[ImaginaryI]],0.0026*Exp[-85.8\[ImaginaryI]]},{0,-0.0241/2,0.0326*Exp[-0.03\[ImaginaryI]]},{0,0,0.9990/2}}];

(* CKM matrix *)
VCKM={{0.9754,0.2205,-0.0026\[ImaginaryI]},{-0.2203*Exp[0.03\[ImaginaryI]],0.9749,0.0318},{0.0075*Exp[-19\[ImaginaryI]],-0.0311*Exp[1.0\[ImaginaryI]],0.9995}};

(* Yukawa coupling matrices *)
Yu=Sqrt[2]/RGEvEW/Sin[N[ArcTan[RGEtan\[Beta]]]]*Mu[MassTop];
Yd=Sqrt[2]/RGEvEW/Cos[N[ArcTan[RGEtan\[Beta]]]]*Md[MassBottom];
Yd2=Sqrt[2]/RGEvEW/Cos[N[ArcTan[RGEtan\[Beta]]]]*Md2[MassBottom];

End[];
End[];


RGEFusaokaYukawa[pModel_,pParticle_,pScale_:"GUT"]:=Block[{},
    Switch[pScale,
	"GUT",
		Switch[pModel,
		  "SM", Switch[pParticle,
		     "u",Return[`SM`GUT`Yu],
		     "d",Return[`SM`GUT`Yd]],
		  "MSSM", Switch[pParticle,
		     "u",Return[`MSSM`GUT`Yu],
		     "d",Return[`MSSM`GUT`Yd]]
		 ],
	"MZ",
		Switch[pModel,
		  "SM", Switch[pParticle,
		     "u",Return[`SM`MZ`Yu],
		     "d",Return[`SM`MZ`Yd]],
		  "MSSM", Switch[pParticle,
		     "u",Return[`SM`MZ`Yu],
		     "d",Return[`SM`MZ`Yd]]
		 ]
	];
];


 End[];

 Protect[RGEFusaokaYukawa];

 EndPackage[];

