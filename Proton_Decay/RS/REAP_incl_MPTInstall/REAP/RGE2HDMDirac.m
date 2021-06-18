(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGE2HDMDirac`",{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGESMDirac`","REAP`RGEMSSMDirac`"}];


(* register 2HDMDirac *)
RGERegisterModel["2HDMDirac","REAP`RGE2HDMDirac`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGEMu->`Private`GetMu,RGE\[Alpha]->`Private`Get\[Alpha],RGEYu->`Private`GetYu,RGERaw->`Private`GetRawSolution,RGEM\[Nu]->`Private`GetM\[Nu],RGEMd->`Private`GetMd,RGEYd->`Private`GetYd,RGEYe->`Private`GetYe,RGEPoleMTop->`Private`GetPoleMTop,RGEMe->`Private`GetMe,RGE\[Lambda]->`Private`Get\[Lambda],RGEAll->`Private`GetRawSolution,RGEY\[Nu]->`Private`GetY\[Nu],RGECoupling->`Private`GetCoupling},
{{"2HDMDirac",`Private`Trans2HDM},{"MSSMDirac",`Private`TransMSSM}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGESMDirac`","REAP`RGEMSSMDirac`"}];
ModelName="2HDMDirac";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of 2HDMDirac *)
(* GetRawY\[Nu] is not a function of 2HDMDirac *)
ClearAll[GetCoupling];
GetCoupling[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the coupling constants *)
   Return[({g1[pScale],g2[pScale],g3[pScale]}/.pSolution)[[1]]];
];

ClearAll[Get\[Alpha]];
Get\[Alpha][pScale_,pSolution_,pOpts___]:=Block[{lg},
(* returns the fine structure constants *)
    lg=({g1[pScale],g2[pScale],g3[pScale]}/.pSolution)[[1]];
    Return[lg^2/(4*Pi)];
];

ClearAll[Get\[Lambda]];
Get\[Lambda][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the higgs couplings *)
   Return[({\[Lambda]1[pScale],\[Lambda]2[pScale],\[Lambda]3[pScale],\[Lambda]4[pScale],\[Lambda]5[pScale]}/.pSolution)[[1]]];
];

(* GetM\[CapitalDelta]2 is not a function of 2HDMDirac *)
(* GetM\[CapitalDelta] is not a function of 2HDMDirac *)
ClearAll[GetYe];
GetYe[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa coupling matrix of the charged leptons *)
    Return[(Ye[pScale]/.pSolution)[[1]]];
];

ClearAll[GetYu];
GetYu[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa coupling matrix of the up-type quarks *)
    Return[(Yu[pScale]/.pSolution)[[1]]];
];

ClearAll[GetYd];
GetYd[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa coupling matrix of the down-type quarks *)
    Return[(Yd[pScale]/.pSolution)[[1]]];
];

(* GetRawY\[CapitalDelta] is not a function of 2HDMDirac *)
ClearAll[GetY\[Nu]];
GetY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa couplings of the neutrinos *)
	Return[(Y\[Nu][pScale]/.pSolution)[[1]]];
];

(* Get\[Kappa] is not a function of 2HDMDirac *)
(* Get\[Kappa]1 is not a function of 2HDMDirac *)
(* Get\[Kappa]2 is not a function of 2HDMDirac *)
(* GetM\[Nu]r is not a function of 2HDMDirac *)
ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu],lM},
(* returns the mass matrix of the neutrinos *)
	lY\[Nu]=(Y\[Nu][pScale]/.pSolution)[[1]];

        lOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];

(*	lvu=RGEvEW/.Options[lOpts,RGEvEW];
	l\[Beta]=N[ArcTan[RGEtan\[Beta]]]/.Options[lOpts,RGEtan\[Beta]];
	\[Nu]=(RGEz\[Nu]/.Options[lOpts,RGEz\[Nu]]).{Cos[l\[Beta]],Sin[l\[Beta]]};
	Return[(l\[Nu]/Sqrt[2]*lY\[Nu])*10^9];*)
	lM=(10^9*RGEv\[Nu]/Sqrt[2]*lY\[Nu])/.Options[lOpts];
	Return[lM];
];

ClearAll[GetMu];
GetMu[pScale_,pSolution_,pOpts___]:=Block[{lMu,lvu,l\[Beta]},
(* returns the mass matrix of the up-type quarks *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   l\[Beta]=ArcTan[RGEtan\[Beta]]/.Options[lOpts,RGEtan\[Beta]];
   lvu=N[(RGEvEW/.Options[lOpts,RGEvEW])*(RGEzu/.Options[lOpts,RGEzu]).{Cos[l\[Beta]],Sin[l\[Beta]]}];
   lMu=lvu/Sqrt[2]*(Yu[pScale]/.pSolution)[[1]];
   Return[lMu];
];

ClearAll[GetPoleMTop];
GetPoleMTop[pScale_,pSolution_,pOpts___]:=Block[{lMtop,lg3,lScale,lPrecision,lCount},
(* returns the mass matrix of the up-type quarks *)
   lPrecision=RGEPrecision/.pOpts/.RGETransition->6;
   lCountMax=RGEMaxNumberIterations/.pOpts/.RGEMaxNumberIterations->20;
   lScale=Exp[pScale];
   lg3=RGEGetSolution[lScale,RGECoupling][[3]];
   lMu=RGEGetSolution[lScale,RGEMu];
   lMtop=Re[Sqrt[Max[Eigenvalues[Dagger[lMu].lMu]]*(1+lg3^2/3/Pi)]];
   lCount=0;
   While[RGEFloor[Abs[lMtop-lScale],RGEPrecision->lPrecision]>0,
	lScale=lMtop;
	lg3=RGEGetSolution[lScale,RGECoupling][[3]];
	lMu=RGEGetSolution[lScale,RGEMu];
	lMtop=Re[Sqrt[Max[Eigenvalues[Dagger[lMu].lMu]]*(1+lg3^2/3/Pi)]];
	lCount++;
	If[lCount>lCountMax,
		Print["RGEGetsolution[pScale,RGEPoleMTop]: algorithm to search transitions does not converge. There have been ",lCount," iterations so far. Returning: ",N[Sort[lTransitions,Greater],lPrecision]];
		Return[lMtop];
            ];
	];
   Return[lMtop];
];

ClearAll[GetMd];
GetMd[pScale_,pSolution_,pOpts___]:=Block[{lMd,lvd,l\[Beta]},
(* returns the mass matrix of the down-type quarks *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   l\[Beta]=ArcTan[RGEtan\[Beta]]/.Options[lOpts,RGEtan\[Beta]];
   lvd=N[(RGEvEW/.Options[lOpts,RGEvEW])*(RGEzd/.Options[lOpts,RGEzd]).{Cos[l\[Beta]],Sin[l\[Beta]]}];
   lMd=lvd/Sqrt[2]*(Yd[pScale]/.pSolution)[[1]];
   Return[lMd];
];

ClearAll[GetMe];
GetMe[pScale_,pSolution_,pOpts___]:=Block[{lMe,lvd,l\[Beta]},
(* returns the mass matrix of the charged leptons *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   l\[Beta]=ArcTan[RGEtan\[Beta]]/.Options[lOpts,RGEtan\[Beta]];
   lvd=N[(RGEvEW/.Options[lOpts,RGEvEW])*Cos[l\[Beta]]];
   lMe=lvd/Sqrt[2]*(Ye[pScale]/.pSolution)[[1]];
   Return[lMe];
];

(* GetSolution is not a function of 2HDMDirac *)
(* GetMixingParameters is not a function of 2HDMDirac *)
(* GetTwistingParameters is not a function of 2HDMDirac *)
(* GetM1Tilde is not a function of 2HDMDirac *)
(* Get\[Epsilon]1Max is not a function of 2HDMDirac *)
(* Get\[Epsilon]1 is not a function of 2HDMDirac *)
(* GetGWCond is not a function of 2HDMDirac *)
(* GetGWConditions is not a function of 2HDMDirac *)
(* GetVEVratio is not a function of 2HDMDirac *)
(* GetVEVratios is not a function of 2HDMDirac *)
ClearAll[Dagger];
RGEvu:=N[RGEvEW*RGEzu.{1/Sqrt[1+RGEtan\[Beta]^2],RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]}];
RGEvd:=N[RGEvEW*RGEzd.{1/Sqrt[1+RGEtan\[Beta]^2],RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]}];
RGEv\[Nu]:=N[RGEvEW*RGEz\[Nu].{1/Sqrt[1+RGEtan\[Beta]^2],RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]}];
RGEve:=N[RGEvEW/Sqrt[1+RGEtan\[Beta]^2]];

ClearAll[Dagger];
(*shortcuts*)
Dagger[x_] := Transpose[Conjugate[x]];

ClearAll[GetParameters];
GetParameters[]:= Block[{},
(* returns the parameters of the model *)
   Return[ParameterSymbols];
];

ClearAll[ModelSetOptions];
ModelSetOptions[pOpts_]:= Block[{},
(* sets the options of the model *)
    SetOptions[RGEOptions,RGEFilterOptions[RGEOptions,pOpts]];
];

ClearAll[ModelGetOptions];
ModelGetOptions[]:= Block[{},
(* returns the options *)
   Return[Options[RGEOptions]];
];

ClearAll[GetInitial];
GetInitial[pOpts___]:=Block[{lSuggestion,lIndexSuggestion,lInitial,lParameters,lParameterRepl},
(* returns the suggested initial values *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lSuggestion=(RGESuggestion/.pOpts)/.{RGESuggestion->"*"};
   lIndexSuggestion=1;
   While[lIndexSuggestion<=Length[Initial] && !StringMatchQ[Initial[[lIndexSuggestion,1]],lSuggestion],	lIndexSuggestion++ ];
   lInitial=Initial[[ lIndexSuggestion,2 ]];
   lParameters=ParameterSymbols/.pOpts/.lInitial/.pOpts/.lInitial/.pOpts/.lInitial/.Options[lOpts];
   Return[Table[ParameterSymbols[[li]]->lParameters[[li]], {li,Length[ParameterSymbols]}]];
   ];

ClearAll[SolveModel];
SolveModel[{pUp_?NumericQ,pUpModel_,pUpOptions_},{pDown_?NumericQ,pDownModel_,pDownOptions_},pDirection_?NumericQ,pBoundary_?NumericQ,pInitial_,pNDSolveOpts_,pOpts___]:=Block[{lSolution,lInitial,lIndexModel,lODE},
(* solves the model and returns the solution, pDown and 0, because the function does not add any models*)
	lNDSolveOpts;
	Options[lNDSolveOpts]=Options[NDSolve];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,Options[RGEOptions]]];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,pOpts]];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,Sequence[pNDSolveOpts]]];

	lOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];

	lInitial=SetInitial[pBoundary,pInitial];
	lIndexModel=Flatten[
	Position[ModelVariants,(RGEModelVariant/.Options[lOpts])] ][[ 1 ]];
	lODE=RGE[[lIndexModel]]/.Options[lOpts];

	lSolution=NDSolve[lODE ~Join~ lInitial, Parameters,{t,pDown,pUp}, Sequence[Options[lNDSolveOpts]]];
	Return[{lSolution,pDown,0}];
];


(* definitions for the Two Higgs Doublet Model (2HDM) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={   RGEModelVariant->"1Loop", (* different variations of a model can be set here *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
			RGEzu->{0,1},
			RGEzd->{1,0},
			RGEz\[Nu]->{0,1}, (* options used to distinguish between different 2HD models *)
                        RGEvEW->246, (* vev of the SM Higgs *)
                        RGEtan\[Beta]->50, (* tan \[Beta]=v2/v1 *)  
			RGE\[Lambda]1->0.75, (* initial value for \[Lambda]1 *)
			RGE\[Lambda]2->0.75, (* initial value for \[Lambda]2 *)
			RGE\[Lambda]3->0.2, (* initial value for \[Lambda]3 *)
			RGE\[Lambda]4->0.2, (* initial value for \[Lambda]4 *)
			RGE\[Lambda]5->0.25, (* initial value for \[Lambda]5 *)
			Method->StiffnessSwitching (* option of NDSolve *)
			}; (* options of the model *)

                        
Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[Nu],\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5};
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[Nu],RGE\[Lambda]1,RGE\[Lambda]2,RGE\[Lambda]3,RGE\[Lambda]4,RGE\[Lambda]5};

(* initial values of 2HDM *)
ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.5828902259929809,
	RGEg2->0.5264896882619359,
	RGEg3->0.5269038670286043,
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGE\[Lambda]1->0.75,
	RGE\[Lambda]2->0.75,
	RGE\[Lambda]3->0.2,
	RGE\[Lambda]4->0.2,
	RGE\[Lambda]5->0.25,
	RGEYd->RGEGetYd[RGEyd,RGEys,RGEyb,RGEq\[Theta]12,RGEq\[Theta]13,RGEq\[Theta]23,RGEq\[Delta],RGEq\[Delta]e,RGEq\[Delta]\[Mu],RGEq\[Delta]\[Tau],RGEq\[CurlyPhi]1,RGEq\[CurlyPhi]2],
	RGEYu->DiagonalMatrix[{RGEyu,RGEyc,RGEyt}],
	RGEq\[Theta]12 -> 12.5216 Degree,
	RGEq\[Theta]13 -> 0.219376 Degree, 
	RGEq\[Theta]23 -> 2.48522 Degree,
	RGEq\[Delta] -> 353.681 Degree,
	RGEq\[CurlyPhi]1 -> 0 Degree,
	RGEq\[CurlyPhi]2 -> 0 Degree, 
	RGEq\[Delta]e -> 0 Degree,
	RGEq\[Delta]\[Mu] -> 0 Degree,
	RGEq\[Delta]\[Tau] -> 0 Degree,
	RGEyu -> 0.94*10^-3*Sqrt[2]/RGEvu,
	RGEyc -> 0.272*Sqrt[2]/RGEvu,
	RGEyt -> 84*Sqrt[2]/RGEvu,
	RGEyd -> 1.94*10^-3*Sqrt[2]/RGEvd,
	RGEys -> 38.7*10^-3*Sqrt[2]/RGEvd,
	RGEyb -> 1.07*Sqrt[2]/RGEvd,
	RGEye -> 0.49348567*10^-3*Sqrt[2]/RGEve,
	RGEy\[Mu] -> 104.15246*10^-3*Sqrt[2]/RGEve,
	RGEy\[Tau] -> 1.7706*Sqrt[2]/RGEve,
	RGEMassHierarchy -> "n",
	RGE\[Theta]12 -> 20 Degree,
	RGE\[Theta]13 -> 0 Degree, 
	RGE\[Theta]23 -> 45 Degree,
	RGE\[Delta] -> 0 Degree,
	RGE\[Delta]e -> 0 Degree,
	RGE\[Delta]\[Mu] -> 0 Degree,
	RGE\[Delta]\[Tau] -> 0 Degree,
	RGE\[CurlyPhi]1 -> 0 Degree,
	RGE\[CurlyPhi]2 -> 0 Degree, 
	RGEMlightest -> 0.05,
	RGE\[CapitalDelta]m2atm -> 3 10^-3, 
	RGE\[CapitalDelta]m2sol -> 1.5 10^-4,
	RGEY\[Nu]->RGEGetDiracY\[Nu][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEv\[Nu]]
	}
},


{"MZ",{
	RGEg1->RGEgMZ[1],
	RGEg2->RGEgMZ[2],
	RGEg3->RGEgMZ[3],
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGEY\[Nu]->RGEGetDiracY\[Nu][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEv\[Nu]],
	RGE\[Lambda]1->0.75,
	RGE\[Lambda]2->0.75,
	RGE\[Lambda]3->0.2,
	RGE\[Lambda]4->0.2,
	RGE\[Lambda]5->0.25,
	RGEYd->RGEGetYd[RGEyd,RGEys,RGEyb,RGEq\[Theta]12,RGEq\[Theta]13,RGEq\[Theta]23,RGEq\[Delta],RGEq\[Delta]e,RGEq\[Delta]\[Mu],RGEq\[Delta]\[Tau],RGEq\[CurlyPhi]1,RGEq\[CurlyPhi]2],
	RGEYu->DiagonalMatrix[{RGEyu,RGEyc,RGEyt}],
	RGEq\[Theta]12 -> 12.7652 Degree, 
	RGEq\[Theta]13 -> 0.170675 Degree,
	RGEq\[Theta]23 -> 2.14069 Degree,
	RGEq\[Delta] -> 0 Degree,
	RGEq\[CurlyPhi]1 -> 0 Degree,
	RGEq\[CurlyPhi]2 -> 0 Degree, 
	RGEq\[Delta]e -> 0 Degree,
	RGEq\[Delta]\[Mu] -> 0 Degree,
	RGEq\[Delta]\[Tau] -> 0 Degree,
	RGEyu -> 2.33*10^-3*Sqrt[2]/RGEvu,
	RGEyc -> 0.677*Sqrt[2]/RGEvu,
	RGEyt -> 181*Sqrt[2]/RGEvu,
	RGEyd -> 4.69*10^-3*Sqrt[2]/RGEvd,
	RGEys -> 93.4*10^-3*Sqrt[2]/RGEvd,
	RGEyb -> 3.00*Sqrt[2]/RGEvd,
	RGEye -> 0.48684727*10^-3*Sqrt[2]/RGEve,
	RGEy\[Mu] -> 0.10275138*Sqrt[2]/RGEve,
	RGEy\[Tau] -> 1.7467*Sqrt[2]/RGEve,
	RGEMassHierarchy -> "n",
	RGE\[Theta]12 -> 33 Degree,
	RGE\[Theta]13 -> 0 Degree, 
	RGE\[Theta]23 -> 45 Degree,
	RGE\[Delta] -> 0 Degree,
	RGE\[Delta]e -> 0 Degree,
	RGE\[Delta]\[Mu] -> 0 Degree,
	RGE\[Delta]\[Tau] -> 0 Degree,
	RGE\[CurlyPhi]1 -> 0 Degree,
	RGE\[CurlyPhi]2 -> 0 Degree, 
	RGEMlightest -> 0.05,
	RGE\[CapitalDelta]m2atm -> 2.5 10^-3, 
	RGE\[CapitalDelta]m2sol -> 8 10^-5
	}
}

}; (* a list containing suggestions for initial values *)


ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[Y\[Nu][t],t]==BetaY\[Nu][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[\[Lambda]1[t],t]==Beta\[Lambda]1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[\[Lambda]2[t],t]==Beta\[Lambda]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[\[Lambda]3[t],t]==Beta\[Lambda]3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[\[Lambda]4[t],t]==Beta\[Lambda]4[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]],
		D[\[Lambda]5[t],t]==Beta\[Lambda]5[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda]1[t],\[Lambda]2[t],\[Lambda]3[t],\[Lambda]4[t],\[Lambda]5[t]]
};

(* Beta functions of the 2HDM *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[Nu], Beta\[Lambda]1,Beta\[Lambda]2,Beta\[Lambda]3,Beta\[Lambda]4,Beta\[Lambda]5];

(* 1 loop contributions *)

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] :=
	(21/5) * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] :=
	(-3) * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] :=
	(-7) * 1/(16*Pi^2) * g3^3;


BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := Block[{li},
	Return[1/(16*Pi^2) * (
          Yd.(
          + 3/2*Dagger[Yd].Yd
	  + (1/2-2*RGEzu.RGEzd)*Dagger[Yu].Yu
          )
          + (
          - (1/4)*g1^2
	  - 9/4*g2^2
	  - 8*g3^2
	  + RGEzd.{Tr[Dagger[Ye].Ye]+3 Tr[Dagger[Yd].Yd],3 Tr[Dagger[Yd].Yd]}
	  + RGEzd.RGEz\[Nu]*Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*RGEzd.RGEzu*Tr[Dagger[Yu].Yu]
          )*Yd
          )]
	  ];

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := Block[{li},
	Return[1/(16*Pi^2) * (
          Yu.(
          + 3/2*Dagger[Yu].Yu
	  + (1/2-2*RGEzu.RGEzd)*Dagger[Yd].Yd
          )
          + (
          - (17/20)*g1^2
	  - 9/4*g2^2
	  - 8*g3^2
	  + RGEzu.{Tr[Dagger[Ye].Ye]+3 Tr[Dagger[Yu].Yu],3 Tr[Dagger[Yu].Yu]}
	  + RGEzu.RGEz\[Nu]*Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*RGEzu.RGEzd*Tr[Dagger[Yd].Yd]
          )*Yu
          )]
	  ];

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := Block[{li},
	  Return[1/(16*Pi^2) * (
          Ye.(
          + 3/2*Dagger[Ye].Ye
	  + (1/2-2*RGEz\[Nu].{1,0})*Dagger[Y\[Nu]].Y\[Nu]
          )
          + (
          - (9/4)*g1^2
	  - 9/4*g2^2
          + Tr[Dagger[Ye].Ye]
	  + (RGEz\[Nu].{1,0})* Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*(RGEzd.{1,0})* Tr[Dagger[Yd].Yd]
	  + 3*(RGEzu.{1,0})* Tr[Dagger[Yu].Yu]
          )*Ye
          )]
	  ];

BetaY\[Nu][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := Block[{li},
	 Return[1/(16*Pi^2) * (
          Y\[Nu].(
          + 3/2*Dagger[Y\[Nu]].Y\[Nu]
	  + (1/2-2*(RGEz\[Nu].{1,0}))*Dagger[Ye].Ye
          )
          + (
          - (9/20)*g1^2
	  - 9/4*g2^2
	  +
	 RGEz\[Nu].{Tr[Dagger[Ye].Ye]+Tr[Dagger[Y\[Nu]].Y\[Nu]],Tr[Dagger[Y\[Nu]].Y\[Nu]]}
	 + 3 RGEz\[Nu].RGEzd *Tr[Dagger[Yd].Yd]
	 + 3 RGEz\[Nu].RGEzu *Tr[Dagger[Yu].Yu]
          )*Y\[Nu]
          )]
	  ];

Beta\[Lambda]1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := 1/(16*Pi^2) * (
		 + 6 * \[Lambda]1^2
		 + 8 * \[Lambda]3^2
		 + 6 * \[Lambda]3 * \[Lambda]4
		 + \[Lambda]5^2
		 - 3 *\[Lambda]1 * (3*g2^2 +3/5*g1^2)
		 + 3 * g2^4
		 + 3/2 * (3/5 *g1^2 + g2^2)^2
		 + 4 \[Lambda]1 * (
		 + Tr[Dagger[Ye].Ye]
		 + (RGEz\[Nu].{1,0}) Tr[Dagger[Y\[Nu]].Y\[Nu]]
		 + 3*(RGEzd.{1,0}) Tr[Dagger[Yd].Yd]
		 + 3*(RGEzu.{1,0}) Tr[Dagger[Yu].Yu]
		 )
		 - 8 * (
		 + Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
		 + (RGEz\[Nu].{1,0})* Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
		 + 3 * (RGEzd.{1,0})* Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
		 + 3 * (RGEzu.{1,0})* Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
		 )
	);


Beta\[Lambda]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := 1/(16*Pi^2) * (
		 + 6 * \[Lambda]2^2
		 + 8 * \[Lambda]3^2
		 + 6 * \[Lambda]3 * \[Lambda]4
		 + \[Lambda]5^2
		 - 3 *\[Lambda]2 * (3*g2^2 +3/5*g1^2)
		 + 3 * g2^4
		 + 3/2 * (3/5 *g1^2 + g2^2)^2
		 + 4 \[Lambda]2 * (
		 + (RGEz\[Nu].{0,1}) Tr[Dagger[Y\[Nu]].Y\[Nu]]
		 + 3*(RGEzd.{0,1}) Tr[Dagger[Yd].Yd]
		 + 3*(RGEzu.{0,1}) Tr[Dagger[Yu].Yu]
		 )
		 - 8 * (
		 + (RGEz\[Nu].{0,1})* Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
		 + 3 * (RGEzd.{0,1})* Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
		 + 3 * (RGEzu.{0,1})* Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
		 )
	);


Beta\[Lambda]3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := 1/(16*Pi^2) * (
		 + (\[Lambda]1 + \[Lambda]2) * (3*\[Lambda]3+\[Lambda]4)
		 + 4 * \[Lambda]3^2
		 + 2 * \[Lambda]4^2
		 + 1/2 * \[Lambda]5^2
		 - 3*\[Lambda]3 * (3*g2^2 +3/5* g1^2)
		 +9/4 * g2^4
		 +27/100 *g1^4
		 -9/10*g1^2*g2^2
		 +4 \[Lambda]3*(
		 + Tr[Dagger[Ye].Ye]
		 +Tr[Dagger[Y\[Nu]].Y\[Nu]]
		 +3*Tr[Dagger[Yd].Yd]
		 +3*Tr[Dagger[Yu].Yu]
		 )
		 -4* (
		 + (RGEz\[Nu].{0,1})*Tr[Dagger[Ye].Ye.Dagger[Y\[Nu]].Y\[Nu]]
		 + 3 *((RGEzd.{1,0})*(RGEzu.{0,1})+(RGEzd.{0,1})*(RGEzu.{1,0})) *Tr[Dagger[Yd].Yd.Dagger[Yu].Yu]
		 )
	);

Beta\[Lambda]4[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := 1/(16*Pi^2) * (
		 + 2*(\[Lambda]1 + \[Lambda]2) * \[Lambda]4
		 + 4* (2*\[Lambda]3+\[Lambda]4)*\[Lambda]4
		 + 8* \[Lambda]5^2
		 -3*\[Lambda]4 * (3g2^2+3/5*g1^2)
		 +9/5*g1^2*g2^2
		 +4 \[Lambda]4*(
		 + Tr[Dagger[Ye].Ye]
		 +Tr[Dagger[Y\[Nu]].Y\[Nu]]
		 +3*Tr[Dagger[Yd].Yd]
		 +3*Tr[Dagger[Yu].Yu]
		 )
		 -4* (
		 + (RGEz\[Nu].{0,1})*Tr[Dagger[Ye].Ye.Dagger[Y\[Nu]].Y\[Nu]]
		 + 3 *((RGEzd.{1,0})*(RGEzu.{0,1})+(RGEzd.{0,1})*(RGEzu.{1,0})) *Tr[Dagger[Yd].Yd.Dagger[Yu].Yu]
		 )
	);


Beta\[Lambda]5[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]1_,\[Lambda]2_,\[Lambda]3_,\[Lambda]4_,\[Lambda]5_] := 1/(16*Pi^2) * (
             + \[Lambda]5 * (
	     + \[Lambda]1
	     + \[Lambda]2
	     + 8*\[Lambda]3
	     + 12*\[Lambda]4
	     -6*(3/5*g1^2+3g2^2)
	     +2*(
	     +Tr[Dagger[Ye].Ye]
	     +Tr[Dagger[Y\[Nu]].Y\[Nu]]
	     +3*Tr[Dagger[Yd].Yd]
	     +3*Tr[Dagger[Yu].Yu]
	     )
	     )
	
	);
	



(* transition functions *)


ClearAll[Trans2HDM];
Trans2HDM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]1,l\[Lambda]2,l\[Lambda]3,l\[Lambda]4,l\[Lambda]5},
(* make a transition from the 2HDM to the 2HDM *)

(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEOptions];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];

        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];


(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]1,l\[Lambda]2,l\[Lambda]3,l\[Lambda]4,l\[Lambda]5}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[Nu]->lY\[Nu],RGE\[Lambda]1->l\[Lambda]1,RGE\[Lambda]2->l\[Lambda]2,RGE\[Lambda]3->l\[Lambda]3,RGE\[Lambda]4->l\[Lambda]4,RGE\[Lambda]5->l\[Lambda]5}];
];


ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]1,l\[Lambda]2,l\[Lambda]3,l\[Lambda]4,l\[Lambda]5,lcb,lsb,lTocb,lTosb,le,lu,ld,l\[Nu]},
(* make a transition from the 2HDM to the 2HDM *)

(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["MSSMDirac"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];

        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];


(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]1,l\[Lambda]2,l\[Lambda]3,l\[Lambda]4,l\[Lambda]5}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=(RGEzu/.Options[lFromOpts,RGEzu]).{lcb,lsb}/lTosb;
	ld=(RGEzd/.Options[lFromOpts,RGEzd]).{lcb,lsb}/lTocb;
	l\[Nu]=(RGEz\[Nu]/.Options[lFromOpts,RGEz\[Nu]]).{lcb,lsb}/lTosb;
	le=lcb/lTocb;

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,RGEY\[Nu]->lY\[Nu]*l\[Nu]}];
];


(* internal functions *)


ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[Nu][pScale],\[Lambda]1[pScale],\[Lambda]2[pScale],\[Lambda]3[pScale],\[Lambda]4[pScale],\[Lambda]5[pScale]};

ClearAll[SetInitial];
SetInitial[pBoundary_?NumericQ,pInitial_]:=Block[{},
(* sets the initial values *)
   Return[		{g1[pBoundary]==RGEg1,
			g2[pBoundary]==RGEg2,
			g3[pBoundary]==RGEg3,
			Yu[pBoundary]==RGEYu,
			Yd[pBoundary]==RGEYd,
			Ye[pBoundary]==RGEYe,
			Y\[Nu][pBoundary]==RGEY\[Nu],
			\[Lambda]1[pBoundary]==RGE\[Lambda]1,
			\[Lambda]2[pBoundary]==RGE\[Lambda]2,
			\[Lambda]3[pBoundary]==RGE\[Lambda]3,
			\[Lambda]4[pBoundary]==RGE\[Lambda]4,
			\[Lambda]5[pBoundary]==RGE\[Lambda]5
			}//.pInitial
			];
];

End[]; (* end of `Private`*)


EndPackage[];
