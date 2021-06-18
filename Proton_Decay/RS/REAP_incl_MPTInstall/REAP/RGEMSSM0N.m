(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEMSSM0N`",{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGESM0N`","REAP`RGE2HDM0N`" (* transtions to these models are possible *)
}];


(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! package depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
(* dependencies are marked by !!! in the file *)

(* register MSSM0N *)
RGERegisterModel["MSSM0N","REAP`RGEMSSM0N`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGE\[Kappa]->`Private`Get\[Kappa],RGEMixingParameters->`Private`GetMixingParameters,RGEMd->`Private`GetMd,RGECoupling->`Private`GetCoupling,RGETwistingParameters->`Private`GetTwistingParameters,RGEM\[Nu]->`Private`GetM\[Nu],RGEM\[Nu]r->`Private`GetM\[Nu]r,RGEYd->`Private`GetYd,RGERaw->`Private`GetRawSolution,RGEPoleMTop->`Private`GetPoleMTop,RGEYe->`Private`GetYe,RGEMe->`Private`GetMe,RGEY\[Nu]->`Private`GetY\[Nu],RGEAll->`Private`GetSolution,RGEMu->`Private`GetMu,RGE\[Alpha]->`Private`Get\[Alpha],RGEYu->`Private`GetYu},
{{"2HDM",`Private`Trans2HDM},{"MSSM0N",`Private`TransMSSM0N},{"2HDM0N",`Private`Trans2HDM0N},{"SM",`Private`TransSM},{"MSSM",`Private`TransMSSM},{"SM0N",`Private`TransSM0N}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGESM0N`","REAP`RGE2HDM0N`"}];

ModelName="MSSM0N";
ModelVariants={"1Loop","2Loop"};
RGE={RGE1Loop,RGE2Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of MSSM0N *)
(* GetRawY\[Nu] is not a function of MSSM0N *)
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

(* Get\[Lambda] is not a function of MSSM0N *)
(* GetM\[CapitalDelta]2 is not a function of MSSM0N *)
(* GetM\[CapitalDelta] is not a function of MSSM0N *)
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

(* GetRawY\[CapitalDelta] is not a function of MSSM0N *)
ClearAll[GetY\[Nu]];
GetY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu]high,lCutoff},
(* returns the mass matrix of the heavy neutrinos *)
            lCutoff=RGEGetCutoff[Exp[pScale],1];
            lY\[Nu]high=RGEGetSolution[lCutoff,RGEY\[Nu],1];
        Return[lY\[Nu]high];
];

ClearAll[Get\[Kappa]];
Get\[Kappa][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa]},
(* returns \[Kappa] *)
	l\[Kappa]=(\[Kappa][pScale]/.pSolution)[[1]];
        Return[l\[Kappa]];
];

(* Get\[Kappa]1 is not a function of MSSM0N *)
(* Get\[Kappa]2 is not a function of MSSM0N *)
ClearAll[GetM\[Nu]r];
GetM\[Nu]r[pScale_,pSolution_,pOpts___]:=Block[{lMhigh,lCutoff},
(* returns the mass matrix of the heavy neutrinos *)
            lCutoff=RGEGetCutoff[Exp[pScale],1];
            lMhigh=RGEGetSolution[lCutoff,RGEM\[Nu]r,1];
        Return[lMhigh];
];

ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa],lvu},
(* returns the mass matrix of the neutrinos *)
	l\[Kappa]=(\[Kappa][pScale]/.pSolution)[[1]];
	lOpts;
	Options[lOpts]=Options[RGEOptions];
	SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
	lvu=RGEv\[Nu]/.Options[lOpts];

	Return[-lvu^2*(1/2)^2*10^9*l\[Kappa]];
];

ClearAll[GetMu];
GetMu[pScale_,pSolution_,pOpts___]:=Block[{lMu,lvu},
(* returns the mass matrix of the up-type quarks *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lvu=RGEvu/.Options[lOpts];
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
GetMd[pScale_,pSolution_,pOpts___]:=Block[{lMd,lvd},
(* returns the mass matrix of the down-type quarks *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lvd=RGEvd/.Options[lOpts];
   lMd=lvd/Sqrt[2]*(Yd[pScale]/.pSolution)[[1]];
   Return[lMd];
];

ClearAll[GetMe];
GetMe[pScale_,pSolution_,pOpts___]:=Block[{lMe,lvd},
(* returns the mass matrix of the charged leptons *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lvd=RGEve/.Options[lOpts];
   lMe=lvd/Sqrt[2]*(Ye[pScale]/.pSolution)[[1]];
   Return[lMe];
];

ClearAll[GetSolution];
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lM,lY\[Nu],l\[Kappa]},
(* returns all parameters of the SM *)
        {lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[pScale]/.pSolution)[[1]];
	lM=GetM\[Nu]r[pScale,pSolution,pOpts];
	lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
        Return[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM}];
];

ClearAll[GetMixingParameters];
GetMixingParameters[pScale_,pSolution_,pOpts___]:=Block[{lM,lYe},
(* returns the leptonic mixing parameters *)
   lM=GetM\[Nu][pScale,pSolution,pOpts];
   lYe=GetYe[pScale,pSolution,pOpts];
   Return[MNSParameters[lM,lYe]];
];

ClearAll[GetTwistingParameters];
GetTwistingParameters[pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu],lYe},
(* returns the mixing parameters of the twisting matrix between Ye^\dagger Ye and Y\nu^\dagger Y\nu *)
   lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
   lYe=GetYe[pScale,pSolution,pOpts];
    Return[MNSParameters[Dagger[lY\[Nu]].lY\[Nu],Dagger[lYe].lYe]];
];

(* GetM1Tilde is not a function of MSSM0N *)
(* Get\[Epsilon]1Max is not a function of MSSM0N *)
(* Get\[Epsilon]1 is not a function of MSSM0N *)
(* GetGWCond is not a function of MSSM0N *)
(* GetGWConditions is not a function of MSSM0N *)
(* GetVEVratio is not a function of MSSM0N *)
(* GetVEVratios is not a function of MSSM0N *)
ClearAll[Dagger];
(*RGEvu:=N[RGEvEW*Sin[ArcTan[RGEtan\[Beta]]]];*)
RGEvu:=N[RGEvEW*RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]];
(*RGEvd:=N[RGEvEW*Cos[ArcTan[RGEtan\[Beta]]]];*)
RGEvd:=N[RGEvEW/Sqrt[1+RGEtan\[Beta]^2]];
RGEv\[Nu]:=RGEvu;
RGEve:=RGEvd;

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
SolveModel[{pUp_?NumericQ,pUpModel_,pUpOptions_},{pDown_?NumericQ,pDownModel_,pDownOptions_},pDirection_?NumericQ,pBoundary_?NumericQ,pInitial_,pNDSolveOpts_,pOpts___]:=Block[{lSolution,lInitial,lODE,lIndexModel},
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
	lIndexModel=Flatten[ Position[ModelVariants,(RGEModelVariant/.Options[lOpts])] ][[ 1 ]];
	lODE=RGE[[lIndexModel]]/.Options[lOpts];

	lSolution=NDSolve[lODE ~Join~ lInitial, Parameters,{t,pDown,pUp}, Sequence[Options[lNDSolveOpts]]];
	Return[{lSolution,pDown,0}];
];




(* definitions for the Minimal Supersymmetric Standard Model (MSSM) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={   RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
		        RGEvEW->246, (* vev of the SM Higgs *)
                        RGEtan\[Beta]->50, (* tan \[Beta]=vu/vd *)  
                        Method->StiffnessSwitching (* option of NDSolve *)
                        }; (* options of the MSSM w/o heavy neutrinos *)


Parameters={g1,g2,g3,Yu,Yd,Ye,\[Kappa]};
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGE\[Kappa]};


ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.7044110331165641,
	RGEg2->0.6965468498179075,
	RGEg3->0.6983661130877465,
	RGEYd->RGEGetYd[RGEyd,RGEys,RGEyb,RGEq\[Theta]12,RGEq\[Theta]13,RGEq\[Theta]23,RGEq\[Delta],RGEq\[Delta]e,RGEq\[Delta]\[Mu],RGEq\[Delta]\[Tau],RGEq\[CurlyPhi]1,RGEq\[CurlyPhi]2],
	RGEYu->DiagonalMatrix[{RGEyu,RGEyc,RGEyt}],
	RGEq\[Theta]12 -> 12.651 Degree,
	RGEq\[Theta]13 -> 0.147249 Degree, 
	RGEq\[Theta]23 -> 1.82387 Degree,
	RGEq\[Delta] -> 293.06 Degree,
	RGEq\[CurlyPhi]1 -> 0 Degree,
	RGEq\[CurlyPhi]2 -> 0 Degree, 
	RGEq\[Delta]e -> 0 Degree,
	RGEq\[Delta]\[Mu] -> 0 Degree,
	RGEq\[Delta]\[Tau] -> 0 Degree,
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
	RGEMlightest -> 0.04,
	RGE\[CapitalDelta]m2atm -> 3*10^-3, 
	RGE\[CapitalDelta]m2sol -> 2 10^-4,
	RGEyu -> 0.00104/RGEvu*Sqrt[2],
	RGEyc -> 0.302/RGEvu*Sqrt[2],
	RGEyt -> 129/RGEvu*Sqrt[2],
	RGEyd -> 0.00133/RGEvd*Sqrt[2],
	RGEys -> 0.0265/RGEvd*Sqrt[2],
	RGEyb ->  1.00/RGEvd*Sqrt[2],
	RGEye -> 0.32502032*10^-3*Sqrt[2]/RGEve,
	RGEy\[Mu] -> 68.59813*10^-3*Sqrt[2]/RGEve,
	RGEy\[Tau] -> 1171.4*10^-3*Sqrt[2]/RGEve,
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGE\[Kappa]->RGEGet\[Kappa][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu]
}
},
{"MZ",{
	RGE\[Kappa] -> RGEGet\[Kappa][RGE\[Theta]12, RGE\[Theta]13, RGE\[Theta]23, RGE\[Delta], RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2, RGEMlightest, RGE\[CapitalDelta]m2atm, RGE\[CapitalDelta]m2sol, RGEMassHierarchy, RGEvu], 
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGEg1 -> RGEgMZ[1],
	RGEg2 -> RGEgMZ[2],
	RGEg3 -> RGEgMZ[3],
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
	RGEyt -> 172.7*Sqrt[2]/RGEvu,
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
}; (* a list containing the suggestions of initial values *)



ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]
}; (* renormalization group equations of the MSSM ( 1 Loop ) *)


(* Beta functions of the MSSM *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, Beta\[Kappa], Beta\[Lambda]];

(* 1 loop contributions *)

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=
	(33/5) * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=
	(1) * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=
	(-3) * 1/(16*Pi^2) * g3^3;


BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] := 1/(16*Pi^2) * (
          Yd.(
          + 3*Dagger[Yd].Yd
	  + Dagger[Yu].Yu
          )
          + (
          - (7/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + 3*Tr[Dagger[Yd].Yd]
	  + Tr[Dagger[Ye].Ye]
          )*Yd
          );

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] := 1/(16*Pi^2) * (
          Yu.(
          + Dagger[Yd].Yd
	  + 3*Dagger[Yu].Yu
          )
          + (
          - (13/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + 3*Tr[Dagger[Yu].Yu]
          )*Yu
          );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] := 1/(16*Pi^2) * (
          Ye.(
          + 3*Dagger[Ye].Ye
          )
          + (
          - (9/5)*g1^2
	  - 3*g2^2
	  + 3*Tr[Dagger[Yd].Yd]
	  + Tr[Dagger[Ye].Ye]
          )*Ye
          );


Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=1/(16*Pi^2) * (
	\[Kappa].(
        + Dagger[Ye].Ye
        )
	+ (
        + Transpose[Ye].Conjugate[Ye]
        ).\[Kappa]
        + (
	- (6/5)*g1^2
	- 6*g2^2
	+ 6*Tr[Dagger[Yu].Yu]
        )*\[Kappa]
        );

ClearAll[RGE2Loop];
RGE2Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]+BetaYu2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]+BetaYd2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]+BetaYe2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]+Beta\[Kappa]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t]]
}; (* renormalization group equations of the MSSM ( 2 Loop ) *)


ClearAll[BetaYu2, BetaYd2, BetaYe2, Beta\[Kappa]2];


(* 2 loop contributions *)

BetaYu2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=1/(4*Pi)^4 * (
	 + Yu.(
         - 2*Dagger[Yd].Yd.Dagger[Yd].Yd
	 - 2*Dagger[Yd].Yd.Dagger[Yu].Yu
	 - 4*Dagger[Yu].Yu.Dagger[Yu].Yu
	 - 3*Dagger[Yd].Yd*Tr[Yd.Dagger[Yd]]
	 - Dagger[Yd].Yd*Tr[Ye.Dagger[Ye]]
	 - 9*Dagger[Yu].Yu*Tr[Yu.Dagger[Yu]]
         + (2/5)*g1^2*Dagger[Yd].Yd
	 + (2/5)*g1^2*Dagger[Yu].Yu
	 + 6*g2^2*Dagger[Yu].Yu
         )
         + Yu*(
	 - 3*Tr[Dagger[Yu].Yd.Dagger[Yd].Yu]
	 - 9*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
         + (4/5)*g1^2*Tr[Dagger[Yu].Yu]
	 + 16*g3^2*Tr[Dagger[Yu].Yu]
         + (2743/450)*g1^4
         + g1^2*g2^2
	 + (15/2)*g2^4
	 + (136/45)*g1^2*g3^2
	 + 8*g2^2*g3^2
	 - (16/9)*g3^4
         )
         );


BetaYd2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=1/(4*Pi)^4 * (
	+Yd.(
        - 4*Dagger[Yd].Yd.Dagger[Yd].Yd
	- 2*Dagger[Yu].Yu.Dagger[Yd].Yd
	- 2*Dagger[Yu].Yu.Dagger[Yu].Yu
	- 9*Dagger[Yd].Yd*Tr[Yd.Dagger[Yd]]
	- 3*Dagger[Yd].Yd*Tr[Ye.Dagger[Ye]]
	- 3*Dagger[Yu].Yu*Tr[Yu.Dagger[Yu]]
	+ 6*Dagger[Yd].Yd*g2^2
	+ (4/5)*g1^2*Dagger[Yd].Yd
	+ (4/5)*g1^2*Dagger[Yu].Yu
        )
	+Yd*(
	- 9*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
	- 3*Tr[Dagger[Yd].Yu.Dagger[Yu].Yd]
	- 3*Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
        - (2/5)*g1^2*Tr[Dagger[Yd].Yd]
	+ (6/5)*g1^2*Tr[Dagger[Ye].Ye]
	+ 16*g3^2*Tr[Dagger[Yd].Yd]
        + (287/90)*g1^4
        + g1^2*g2^2
	+ (15/2)*g2^4
	+ (8/9)*g1^2*g3^2
	+ 8*g2^2*g3^2
	- (16/9)*g3^4
        )
        );


BetaYe2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=1/(4*Pi)^4 * (
	+ Ye.(
        - 4*Dagger[Ye].Ye.Dagger[Ye].Ye
	- 9*Dagger[Ye].Ye*Tr[Yd.Dagger[Yd]]
	- 3*Dagger[Ye].Ye*Tr[Ye.Dagger[Ye]]
	+ 6*g2^2*Dagger[Ye].Ye
        )
        +Ye*(
	- 9*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
	- 3*Tr[Dagger[Yd].Yu.Dagger[Yu].Yd]
	- 3*Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
	+ (6/5)*g1^2*Tr[Dagger[Ye].Ye]
        - (2/5)*g1^2*Tr[Dagger[Yd].Yd]
	+ 16*g3^2*Tr[Dagger[Yd].Yd]
	+ (27/2)*g1^4
        + (9/5)*g1^2*g2^2
	+ (15/2)*g2^4
        )
        );


Beta\[Kappa]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_] :=1/(4*Pi)^4 * (
	+ \[Kappa].(
        - 2*Dagger[Ye].Ye.Dagger[Ye].Ye
	+ ((6/5)*g1^2 - Tr[Ye.Dagger[Ye]]- 3*Tr[Yd.Dagger[Yd]])*Transpose[Ye].Conjugate[Ye]
        )
        + (
	- 2*Transpose[Ye].Conjugate[Ye].Transpose[Ye].Conjugate[Ye]
	+ ((6/5)*g1^2 - Tr[Ye.Dagger[Ye]]- 3*Tr[Yd.Dagger[Yd]])*Transpose[Ye].Conjugate[Ye]
        ).\[Kappa]
	+ \[Kappa]*(
	- 6*Tr[Dagger[Yu].Yd.Dagger[Yd].Yu]
	- 18*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
        + (8/5)*g1^2*Tr[Dagger[Yu].Yu]
	+ 32*g3^2*Tr[Dagger[Yu].Yu]
        + (207/25)*g1^4
	+ (18/5)*g1^2*g2^2
	+ 15*g2^4
        )
        );





	(* transition functions *)

ClearAll[TransMSSM0N];
TransMSSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa]}];
];

ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];

ClearAll[Trans2HDM0N];
Trans2HDM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lsb,lcb,lTosb,lTocb,lu,le,ld,l\[Nu]},
(* make a transition from the MSSM to the 2HDM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=lsb/((RGEzu/.Options[lToOpts,RGEzu]).{lTocb,lTosb});
	ld=lcb/((RGEzd/.Options[lToOpts,RGEzd]).{lTocb,lTosb});
	l\[Nu]=lsb/((RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]]).{lTocb,lTosb});
	le=lcb/lTocb;

	
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,
	RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,
	RGE\[Kappa]1->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[1]]*l\[Kappa]*l\[Nu]^2,
	RGE\[Kappa]2->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[2]]*l\[Kappa]*l\[Nu]^2,RGE\[Lambda]1->(lg2^2+lg1^2)/2.,RGE\[Lambda]2->(lg2^2+lg1^2)/2.,RGE\[Lambda]3->(lg2^2-lg1^2)/4,RGE\[Lambda]4->-(lg2^2)/2,RGE\[Lambda]5->0.}];
];

ClearAll[Trans2HDM];
Trans2HDM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lTosb,lTocb,lu,le,ld,l\[Nu]},

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
	lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=lsb/((RGEzu/.Options[lToOpts,RGEzu]).{lTocb,lTosb});
	ld=lcb/((RGEzd/.Options[lToOpts,RGEzd]).{lTocb,lTosb});
	l\[Nu]=lsb/((RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]]).{lTocb,lTosb});
	le=lcb/lTocb;

	
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,
	RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,
	RGE\[Kappa]1->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[1]]*l\[Kappa]*l\[Nu]^2,
	RGE\[Kappa]2->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[2]]*l\[Kappa]*l\[Nu]^2,
	RGEM\[Nu]r->{},RGE\[Lambda]1->(lg2^2+lg1^2)/2.,RGE\[Lambda]2->(lg2^2+lg1^2)/2.,RGE\[Lambda]3->(lg2^2-lg1^2)/4,RGE\[Lambda]4->-(lg2^2)/2,RGE\[Lambda]5->0.}];
];


ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb},
(* make a transition from the MSSM to the SM w/o \[Nu] *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
	lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}]; (* !!! Parameters !!!*)
];


ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb},
(* make a transition from the MSSM to the SM w/o \[Nu] *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
	lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->{},RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}]; (* !!! Parameters !!!*)
];


(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],\[Kappa][pScale]};
ClearAll[SetIntial];
SetInitial[pBoundary_?NumericQ,pInitial_]:=Block[{},
(* sets the initial values *)
   Return[		{g1[pBoundary]==RGEg1,
			g2[pBoundary]==RGEg2,
			g3[pBoundary]==RGEg3,
			Yu[pBoundary]==RGEYu,
			Yd[pBoundary]==RGEYd,
			Ye[pBoundary]==RGEYe,
			\[Kappa][pBoundary]==RGE\[Kappa]
			}//.pInitial
			];
];


End[]; (* end of `Private`*)


EndPackage[];
