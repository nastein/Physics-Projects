(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGESMDirac`",{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEUtilities`","REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGEMSSMDirac`","REAP`RGE2HDMDirac`" (* transitions to these models are possible *)
}];

(* register SMDirac *)
RGERegisterModel["SMDirac","REAP`RGESMDirac`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGE\[Alpha]->`Private`Get\[Alpha],RGEMu->`Private`GetMu,RGEMd->`Private`GetMd,RGEPoleMTop->`Private`GetPoleMTop,RGECoupling->`Private`GetCoupling,RGERaw->`Private`GetRawSolution,RGEYd->`Private`GetYd,RGEY\[Nu]->`Private`GetY\[Nu],RGEYu->`Private`GetYu,RGE\[Lambda]->`Private`Get\[Lambda],RGEYe->`Private`GetYe,RGEMe->`Private`GetMe,RGEM\[Nu]->`Private`GetM\[Nu],RGEAll->`Private`GetRawSolution},
{{"MSSMDirac",`Private`TransMSSM},{"SMDirac",`Private`TransSM}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEUtilities`","REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGEMSSMDirac`","REAP`RGE2HDMDirac`"}];

ModelName="SMDirac";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of SMDirac *)
(* GetRawY\[Nu] is not a function of SMDirac *)
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

(* GetM\[CapitalDelta]2 is not a function of SMDirac *)
(* GetM\[CapitalDelta] is not a function of SMDirac *)
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

(* GetRawY\[CapitalDelta] is not a function of SMDirac *)
ClearAll[GetY\[Nu]];
GetY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa couplings of the neutrinos *)
	Return[(Y\[Nu][pScale]/.pSolution)[[1]]];
];

(* Get\[Kappa] is not a function of SMDirac *)
(* Get\[Kappa]1 is not a function of SMDirac *)
(* Get\[Kappa]2 is not a function of SMDirac *)
(* GetM\[Nu]r is not a function of SMDirac *)
ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu],lvu},
(* returns the mass matrix of the neutrinos *)
	lY\[Nu]=(Y\[Nu][pScale]/.pSolution)[[1]];

        lOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lvu=RGEv\[Nu]/.Options[lOpts];
	Return[(lvu/Sqrt[2]*lY\[Nu])*10^9];
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

(* GetSolution is not a function of SMDirac *)
(* GetMixingParameters is not a function of SMDirac *)
(* GetTwistingParameters is not a function of SMDirac *)
(* GetM1Tilde is not a function of SMDirac *)
(* Get\[Epsilon]1Max is not a function of SMDirac *)
(* Get\[Epsilon]1 is not a function of SMDirac *)
(* GetGWCond is not a function of SMDirac *)
(* GetGWConditions is not a function of SMDirac *)
(* GetVEVratio is not a function of SMDirac *)
(* GetVEVratios is not a function of SMDirac *)
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
Options[RGEOptions]={     RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
			  RGEvEW->246, (* vev for the electroweak transition *)
			  RGE\[Lambda]->0.5, (* initial value for \[Lambda] *)
			  Method->StiffnessSwitching (* option of NDSolve *)
};


Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[Nu],\[Lambda]}; (* These are the parameters of the model *)
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[Nu],RGE\[Lambda]};


ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.5787925294736758,
	RGEg2->0.5214759925514961,
	RGEg3->0.5269038649895842,
	RGE\[Lambda]->0.5,
	RGEY\[Nu]->RGEGetDiracY\[Nu][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu],
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
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
	RGE\[CapitalDelta]m2atm -> 2*10^-3, 
	RGE\[CapitalDelta]m2sol -> 10^-4
}
},
{"MZ",{
	RGEY\[Nu]->RGEGetDiracY\[Nu][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu],
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
}; (* This list contains the suggestion for the initial values *)


ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[Y\[Nu][t],t]==BetaY\[Nu][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]],
		D[\[Lambda][t],t]==Beta\[Lambda][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Lambda][t]]}; (* These are the Renormalization group equations *)

              
(* Beta Functions of the Standardmodel *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[Nu], Beta\[Lambda]];

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=
	41/10 * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=
	-19/6 * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=
	-7 * 1/(16*Pi^2) * g3^3;

BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] := 1/(16*Pi^2) * (
	(3/2)*Yd.(
        + Dagger[ Yd ].Yd
	- Dagger[ Yu ].Yu
        )
	+ (
	- (1/4)*g1^2
	- (9/4)*g2^2
	- 8*g3^2
	+ 3*Tr[Dagger[ Yd ].Yd ]
	+ 3*Tr[Dagger[ Yu ].Yu ]
	+ Tr[Dagger[ Ye ].Ye ]
	+ Tr[Dagger[ Y\[Nu] ].Y\[Nu] ]
        )*Yd
	);

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=1/(16*Pi^2) * (
       (3/2)*Yu.(
       - Dagger[Yd].Yd
       + Dagger[Yu].Yu
       )
       + (
       - (17/20)*g1^2
       - (9/4)*g2^2
       - 8*g3^2
       + Tr[Dagger[Y\[Nu]].Y\[Nu]]
       + Tr[Dagger[Ye].Ye]
       + 3*Tr[Dagger[Yu].Yu]
       + 3*Tr[Dagger[Yd].Yd]
       )*Yu
       );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=1/(16*Pi^2) * (
	(3/2)*Ye.(
        Dagger[Ye].Ye
	- Dagger[Y\[Nu]].Y\[Nu]
        )
	+ (
	- (9/4)*g1^2
	- (9/4)*g2^2
	+ 3*Tr[Dagger[ Yd ].Yd ]
	+ 3*Tr[Dagger[ Yu ].Yu ]
	+ Tr[Dagger[ Ye ].Ye ]
	+ Tr[Dagger[ Y\[Nu] ].Y\[Nu] ]
        )*Ye
	);

      
BetaY\[Nu][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] :=1/(16*Pi^2) * (
	(3/2)*Y\[Nu].(
        - Dagger[Ye].Ye
	+ Dagger[Y\[Nu]].Y\[Nu]
        )
	+ (
	- (9/20)*g1^2
	- (9/4)*g2^2
	+ Tr[Dagger[Y\[Nu]].Y\[Nu]]
	+ Tr[Dagger[Ye].Ye]
	+ 3*Tr[Dagger[Yu].Yu]
	+ 3*Tr[Dagger[Yd].Yd]
        )*Y\[Nu]
	);


Beta\[Lambda][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Lambda]_] := 1/(16*Pi^2) * (
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
    + Tr[Dagger[Y\[Nu]].Y\[Nu]]
    )
    - 8*(
    + 3*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
    + 3*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
    + Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
    + Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
    )
    );

ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]},
(* make a transition from the SM to the SM *)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[Nu]->lY\[Nu],RGE\[Lambda]->l\[Lambda]}];
];

ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda],l\[Beta]},
(* make a transition from the SM to the SM *)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["MSSMDirac"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Lambda]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	l\[Beta]=N[ArcTan[RGEtan\[Beta]]]/.Options[lToOpts,RGEtan\[Beta]];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu/Sin[l\[Beta]],RGEYd->lYd/Cos[l\[Beta]],RGEYe->lYe/Cos[l\[Beta]],RGEY\[Nu]->lY\[Nu]/Sin[l\[Beta]]}];
];


(* internal functions *)


ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[Nu][pScale],\[Lambda][pScale]};


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
			\[Lambda][pBoundary]==RGE\[Lambda]
			}//.pInitial
			];
];


End[]; (* end of `Private` *)


EndPackage[];
