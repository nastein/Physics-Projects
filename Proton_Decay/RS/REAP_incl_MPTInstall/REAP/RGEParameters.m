(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEParameters`"];

ClearAll[RGEMass];

RGEMass::usage="RGEMass[particle name] returns the mass of the particle at the Z
boson mass. The possible particles are u,d,s,c,b,t,e,\[Mu],\[Tau] and Z.";


ClearAll[gMZ];

RGEgMZ::usage="RGEgMZ[coupling] returns the value of the coupling constant <coupling> (coupling=1,2,3) at the mass of the Z boson. Be careful RGEgMZ[1] is the coupling constant of U(1) Hypercharge with GUT charge normalization, i.e. RGEgMZ[1]=Sqrt[5/3]*g1";



 Begin["`Private`"];

 Map[Needs,{}];
ThetaW = ArcSin[Sqrt[0.231]];(*Weinberg angle*)
ThetaC = ThetaQuark12; (*Cabibbo Angle*)
ThetaQuark12 = ArcSin[0.223];
ThetaQuark13 = ArcSin[0.0036];
ThetaQuark23 = ArcSin[0.041];
DeltaQuark = 1.02;


RGEMass[s_] := 
	Which[
		s == "u", 4*10^(-3), 
		s == "d", 7*10^(-3),
		s == "s", 125*10^(-3), 
		s == "c", 1.25,
		s == "b", 4.2, 
		s == "t", 174, 
		s == "e", 511*10^(-6), 
		s == "\<\[Mu]\>", 105.7*10^(-3), 
		s == "\<\[Tau]\>", 1.777,
		s == "Z", 91.19
		]; 


RGEgMZ[i_] := 
	Which[	  i == 1, Sqrt[5/3]*Sqrt[4*Pi/(128*Cos[ThetaW]^2)], 
		  i == 2, Sqrt[4*Pi/(128*Sin[ThetaW]^2)], 
		  i == 3, Sqrt[4*Pi*0.118]];


 End[];

 Protect[RGEMass,gMZ];

 EndPackage[];

