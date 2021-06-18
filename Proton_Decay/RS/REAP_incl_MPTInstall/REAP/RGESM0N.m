(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGESM0N`",{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEParameters`","REAP`RGEUtilities`","MixingParameterTools`MPT3x3`","REAP`RGEInitial`","REAP`RGEMSSM0N`"}];


(* register SM0N *)
RGERegisterModel["SM0N","REAP`RGESM0N`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGEAll->`Private`GetSolution,RGEMixingParameters->`Private`GetMixingParameters,RGERaw->`Private`GetRawSolution,RGECoupling->`Private`GetCoupling,RGEM\[Nu]r->`Private`GetM\[Nu]r,RGEMu->`Private`GetMu,RGEYd->`Private`GetYd,RGEMe->`Private`GetMe,RGEYu->`Private`GetYu,RGEMd->`Private`GetMd,RGE\[Kappa]->`Private`Get\[Kappa],RGEYe->`Private`GetYe,RGEY\[Nu]->`Private`GetY\[Nu],RGEPoleMTop->`Private`GetPoleMTop,RGE\[Lambda]->`Private`Get\[Lambda],RGEM\[Nu]->`Private`GetM\[Nu],RGETwistingParameters->`Private`GetTwistingParameters,RGE\[Alpha]->`Private`Get\[Alpha]},
{{"SM0N",`Private`TransSM0N},{"MSSM",`Private`TransMSSM},{"SM",`Private`TransSM},{"MSSM0N",`Private`TransMSSM0N}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEParameters`","REAP`RGEUtilities`","MixingParameterTools`MPT3x3`","REAP`RGEInitial`","REAP`RGEMSSM0N`"}];

ModelName="SM0N";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of SM0N *)
(* GetRawY\[Nu] is not a function of SM0N *)
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
(* returns the coupling constants *)
   Return[(\[Lambda][pScale]/.pSolution)[[1]]];
];

(* GetM\[CapitalDelta]2 is not a function of SM0N *)
(* GetM\[CapitalDelta] is not a function of SM0N *)
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

(* GetRawY\[CapitalDelta] is not a function of SM0N *)
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

(* Get\[Kappa]1 is not a function of SM0N *)
(* Get\[Kappa]2 is not a function of SM0N *)
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
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lM,lY\[Nu],l\[Kappa],l\[Lambda]},
(* returns all parameters of the SM *)
        {lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]}=(ParametersFunc[pScale]/.pSolution)[[1]];
	lM=GetM\[Nu]r[pScale,pSolution,pOpts];
	lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
        Return[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM,l\[Lambda]}];
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

(* GetM1Tilde is not a function of SM0N *)
(* Get\[Epsilon]1Max is not a function of SM0N *)
(* Get\[Epsilon]1 is not a function of SM0N *)
(* GetGWCond is not a function of SM0N *)
(* GetGWConditions is not a function of SM0N *)
(* GetVEVratio is not a function of SM0N *)
(* GetVEVratios is not a function of SM0N *)
ClearAll[Dagger];
RGEvu:=N[RGEvEW];
RGEvd:=RGEvu;
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


(* definitions for the Standard Model (SM) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={		  RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
				  RGEvEW->246, (* vev for the electroweak transition *)
				  RGE\[Lambda]->0.5, (* initial value for \[Lambda] *)
				  Method->StiffnessSwitching  (* option of NDSolve *)
};

Parameters={g1,g2,g3,Yu,Yd,Ye,\[Kappa],\[Lambda]}; 
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGE\[Kappa],RGE\[Lambda]};


ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.5787925294736758,
	RGEg2->0.5214759925514961,
	RGEg3->0.5269038649895842,
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
	RGE\[CapitalDelta]m2atm -> 5*10^-3, 
	RGE\[CapitalDelta]m2sol -> 2 10^-4,
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGE\[Kappa] -> RGEGet\[Kappa][RGE\[Theta]12, RGE\[Theta]13, RGE\[Theta]23, RGE\[Delta], RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2, RGEMlightest, RGE\[CapitalDelta]m2atm, RGE\[CapitalDelta]m2sol, RGEMassHierarchy, RGEvu], 
	RGE\[Lambda]->0.5
	}
},
{"MZ",{
	RGE\[Kappa] -> RGEGet\[Kappa][RGE\[Theta]12, RGE\[Theta]13, RGE\[Theta]23, RGE\[Delta], RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2, RGEMlightest, RGE\[CapitalDelta]m2atm, RGE\[CapitalDelta]m2sol, RGEMassHierarchy, RGEvu], 
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGEg1 -> RGEgMZ[1],
	RGEg2 -> RGEgMZ[2],
	RGEg3 -> RGEgMZ[3],
	RGE\[Lambda] -> 0.5,
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
};


ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]],
		D[\[Lambda][t],t]==Beta\[Lambda][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],\[Kappa][t],\[Lambda][t]]};


(* Beta Functions of the Standardmodel *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, Beta\[Kappa], Beta\[Lambda]];


Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=
	41/10 * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=
	-19/6 * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=
	-7 * 1/(16*Pi^2) * g3^3;

BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] := 1/(16*Pi^2) * (
	(3/2)*Yd.Dagger[ Yd ].Yd
	- (3/2)*Yd.Dagger[ Yu ].Yu
	+ (
	- (1/4)*g1^2
	- (9/4)*g2^2
	- 8*g3^2
	+ 3*Tr[Dagger[ Yd ].Yd ]
	+ 3*Tr[Dagger[ Yu ].Yu ]
	+ Tr[Dagger[ Ye ].Ye ])*Yd
	);

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=1/(16*Pi^2) * (
       (-3/2)*Yu.Dagger[Yd].Yd
       + (3/2)*Yu.Dagger[Yu].Yu
       + (
       - (17/20)*g1^2
       - (9/4)*g2^2
       - 8*g3^2
       + Tr[Dagger[Ye].Ye]
       + 3*Tr[Dagger[Yu].Yu]
       + 3*Tr[Dagger[Yd].Yd])*Yu
       );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=1/(16*Pi^2) * (
	(3/2)*Ye.Dagger[Ye].Ye
	+ (
	- (9/4)*g1^2
	- (9/4)*g2^2
	+ 3*Tr[Dagger[ Yd ].Yd ]
	+ 3*Tr[Dagger[ Yu ].Yu ]
	+ Tr[Dagger[ Ye ].Ye ])*Ye
	);

      
Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] :=1/(16*Pi^2) * (
	(-3/2)*\[Kappa].Dagger[Ye].Ye
	- (3/2)*Transpose[Ye].Conjugate[Ye].\[Kappa]
	+ (
	- 3*g2^2
	+ 6*Tr[Dagger[Yu].Yu]
	+ 6*Tr[Dagger[Yd].Yd]
	+ 2*Tr[Dagger[Ye].Ye]
	+ \[Lambda])*\[Kappa]
	);


Beta\[Lambda][g1_,g2_,g3_,Yu_,Yd_,Ye_,\[Kappa]_,\[Lambda]_] := 1/(16*Pi^2) * (
    6*\[Lambda]^2
    - (9/5)*g1^2*\[Lambda]
    - 9*g2^2*\[Lambda]
    + (27/50)*g1^4
    + (18/10)*g1^2*g2^2
    + (9/2)*g2^4
    + 4*\[Lambda]*(
    + 3*Tr[Dagger[Yu].Yu]
    + 3*Tr[Dagger[Yd].Yd]
    + Tr[Dagger[Ye].Ye]
    )
    - 8*(
    + 3*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
    + 3*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
    + Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
    )
    );
    
ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda]}];
];

ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda],RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];

ClearAll[TransMSSM0N];
TransMSSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["MSSM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	l\[Beta]=N[ArcTan[RGEtan\[Beta]]]/.Options[lToOpts,RGEtan\[Beta]];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu/Sin[l\[Beta]],RGEYd->lYd/Cos[l\[Beta]],RGEYe->lYe/Cos[l\[Beta]],RGE\[Kappa]->l\[Kappa]/Sin[l\[Beta]]^2}];
];


ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["MSSM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	l\[Beta]=N[ArcTan[RGEtan\[Beta]]]/.Options[lToOpts,RGEtan\[Beta]];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu/Sin[l\[Beta]],RGEYd->lYd/Cos[l\[Beta]],RGEYe->lYe/Cos[l\[Beta]],RGE\[Kappa]->l\[Kappa]/Sin[l\[Beta]]^2,RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];


(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],\[Kappa][pScale],\[Lambda][pScale]};



ClearAll[SetInitial];
SetInitial[pBoundary_?NumericQ,pInitial_]:=Block[{},
(* sets the initial values *)
   Return[		{g1[pBoundary]==RGEg1,
			g2[pBoundary]==RGEg2,
			g3[pBoundary]==RGEg3,
			Yu[pBoundary]==RGEYu,
			Yd[pBoundary]==RGEYd,
			Ye[pBoundary]==RGEYe,
			\[Kappa][pBoundary]==RGE\[Kappa],
			\[Lambda][pBoundary]==RGE\[Lambda]
			}//.pInitial
			];
];


End[]; (* end of `Private` *)


EndPackage[];
