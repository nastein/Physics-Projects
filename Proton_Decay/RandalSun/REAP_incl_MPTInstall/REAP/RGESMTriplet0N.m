(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGESMTriplet0N`",{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEUtilities`","REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGESMTriplet`", "REAP`RGESM`", "REAP`RGESM0N`" (* transtions to these models are possible *)
}];



(* register SMTriplet0N *)
RGERegisterModel["SMTriplet0N","REAP`RGESMTriplet0N`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGEAll->`Private`GetSolution,RGEPoleMTop->`Private`GetPoleMTop,RGEMe->`Private`GetMe,RGE\[Kappa]->`Private`Get\[Kappa],RGERaw->`Private`GetRawSolution,RGEY\[Nu]->`Private`GetY\[Nu],RGEMixingParameters->`Private`GetMixingParameters,RGETwistingParameters->`Private`GetTwistingParameters,RGERawM\[Nu]r->`Private`GetRawM\[Nu]r,RGEM\[CapitalDelta]2->`Private`GetM\[CapitalDelta]2,RGE\[Epsilon]->`Private`Get\[Epsilon],RGECoupling->`Private`GetCoupling,RGE\[Alpha]->`Private`Get\[Alpha],RGEMd->`Private`GetMd,RGEYe->`Private`GetYe,RGEMu->`Private`GetMu,RGEM\[Nu]r->`Private`GetM\[Nu]r,RGERawY\[Nu]->`Private`GetRawY\[Nu],RGEYd->`Private`GetYd,RGEM1Tilde->`Private`GetM1Tilde,RGEM\[Nu]->`Private`GetM\[Nu],RGE\[Lambda]->`Private`Get\[Lambda],RGE\[Epsilon]Max->`Private`Get\[Epsilon]Max,RGE\[Epsilon]1Max->`Private`Get\[Epsilon]1Max,RGEYu->`Private`GetYu,RGERawY\[CapitalDelta]->`Private`GetRawY\[CapitalDelta],RGE\[Epsilon]1->`Private`Get\[Epsilon]1},
{{"SMTriplet",`Private`TransSMTriplet},{"SM0N",`Private`TransSM0N},{"SM",`Private`TransSM},{"SMTriplet0N",`Private`TransSMTriplet0N}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEUtilities`","REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGESMTriplet`", "REAP`RGESM`", "REAP`RGESM0N`"}];

ModelName="SMTriplet0N";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

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

ClearAll[Get\[Lambda]];
Get\[Lambda][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the coupling constants *)
   Return[({\[Lambda][pScale],\[CapitalLambda]1[pScale],\[CapitalLambda]2[pScale],\[CapitalLambda]4[pScale],\[CapitalLambda]5[pScale],\[CapitalLambda]6[pScale],M\[CapitalDelta]2[pScale]}/.pSolution)[[1]]];
];

ClearAll[GetM\[CapitalDelta]2];
GetM\[CapitalDelta]2[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the coupling constants *)
   Return[(M\[CapitalDelta]2[pScale]/.pSolution)[[1]]];
];

(* GetM\[CapitalDelta] is not a function of SMTriplet0N *)
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

ClearAll[GetRawY\[CapitalDelta]];
GetRawY\[CapitalDelta][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa coupling matrix of the down-type quarks *)
    Return[(Y\[CapitalDelta][pScale]/.pSolution)[[1]]];
];

ClearAll[GetY\[Nu]];
GetY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu]high,lCutoff},
(* returns the mass matrix of the heavy neutrinos *)
            lCutoff=RGEGetCutoff[Exp[pScale],1];
            lY\[Nu]high=RGEGetSolution[lCutoff,RGEY\[Nu],1];
        Return[lY\[Nu]high];
];

ClearAll[GetY\[Nu]];
GetY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{lY\[Nu],lY\[Nu]Rotated,lM,lIntegratedOut,lLenM,lMhigh,lf,lg,lMEval,lCutoff,lUp},
(* returns the Yukawa couplings of the neutrinos *)
        lOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lIntegratedOut=(RGEIntegratedOut/.Options[lOpts]);

	lY\[Nu]=(Y\[Nu][pScale]/.pSolution)[[1]];

        Catch[
          If[lIntegratedOut>0,
            lCutoff=RGEGetCutoff[Exp[pScale],1];
	    lUp=RGEGetRange[][[2]];
	    If[lCutoff>Exp[lUp],
		{lMInitial,lY\[Nu]Rotated}=({RGEM\[Nu]r,RGEY\[Nu]}/.RGEGetInitial[][[2]]);
		lLenM=Length[lMInitial];
		If[lLenM>Length[lM],
		            {lMInitial,lY\[Nu]Rotated}=RGERotateM[lMInitial,lY\[Nu]Rotated];
			    lY\[Nu]=Table[If[lf<=lLenM-lIntegratedOut,lY\[Nu][[lf,lg]],lY\[Nu]Rotated[[lf,lg]]],{lf,lLenM},{lg,lLenM}];
		];

		Throw[lCutoff,RGEScaleTooBig]];
            lY\[Nu]Rotated=RGEGetSolution[lCutoff,RGEY\[Nu],1];
            lM=RGEGetSolution[lCutoff,RGEM\[Nu]r,1];
            lLenM=Length[lM];
            {lM,lY\[Nu]Rotated}=RGERotateM[lM,lY\[Nu]Rotated];
            lY\[Nu]=Table[If[lf<=lLenM-lIntegratedOut,lY\[Nu][[lf,lg]],lY\[Nu]Rotated[[lf,lg]]],{lf,lLenM},{lg,lLenM}];
          ];
        ,RGEScaleTooBig];
	Return[lY\[Nu]];
];

ClearAll[Get\[Kappa]];
Get\[Kappa][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa]},
(* returns \[Kappa] *)
	l\[Kappa]=(\[Kappa][pScale]/.pSolution)[[1]];
        Return[l\[Kappa]];
];

(* Get\[Kappa]1 is not a function of SMTriplet0N *)
(* Get\[Kappa]2 is not a function of SMTriplet0N *)
ClearAll[GetM\[Nu]r];
GetM\[Nu]r[pScale_,pSolution_,pOpts___]:=Block[{lMhigh,lCutoff},
(* returns the mass matrix of the heavy neutrinos *)
            lCutoff=RGEGetCutoff[Exp[pScale],1];
            lMhigh=RGEGetSolution[lCutoff,RGEM\[Nu]r,1];
        Return[lMhigh];
];

ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]6,lM\[CapitalDelta]2,lvu},
(* returns the mass matrix of the neutrinos *)
	{l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]6,lM\[CapitalDelta]2}=({\[Kappa][pScale],Y\[CapitalDelta][pScale],\[CapitalLambda]6[pScale],M\[CapitalDelta]2[pScale]}/.pSolution)[[1]];
	lOpts;
	Options[lOpts]=Options[RGEOptions];
	SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
	lvu=RGEv\[Nu]/.Options[lOpts];

	Return[-lvu^2*(1/2)^2*(10^9*l\[Kappa]-2 *l\[CapitalLambda]6/(lM\[CapitalDelta]2/10^9)*lY\[CapitalDelta])];
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
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lM,lY\[Nu],l\[Kappa],l\[Lambda],lY\[CapitalDelta],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2},
(* returns all parameters of the SM *)
        {lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}=(ParametersFunc[pScale]/.pSolution)[[1]];
	lM=GetM\[Nu]r[pScale,pSolution,pOpts];
	lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
        Return[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM,l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}];
];

(* GetGWCond is not a function of SMTriplet0N *)
(* GetGWConditions is not a function of SMTriplet0N *)
(* GetVEVratio is not a function of SMTriplet0N *)
(* GetVEVratios is not a function of SMTriplet0N *)
ClearAll[ModelName0N];
ModelName0N=ModelName<>"0N";

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



(* definitions for the Standard Model (SM0N) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={     RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEPrecision->6, (* precision to find transitions *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
                        RGEMaxNumberIterations->20, (* maximum number of iterations in the loops to search transitions *)
                        RGEvEW->246, (* vev of the SM Higgs *)
			Method->StiffnessSwitching (* option of NDSolve *)
};

Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[CapitalDelta],\[Kappa],\[Lambda],\[CapitalLambda]1,\[CapitalLambda]2,\[CapitalLambda]4,\[CapitalLambda]5,\[CapitalLambda]6, M\[CapitalDelta]2}; (* These are the parameters of the model *)
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[CapitalDelta],RGE\[Kappa],RGE\[Lambda],RGE\[CapitalLambda]1,RGE\[CapitalLambda]2,RGE\[CapitalLambda]4,RGE\[CapitalLambda]5,RGE\[CapitalLambda]6, RGEM\[CapitalDelta]2};


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
        RGEY\[CapitalDelta]->-RGEM\[CapitalDelta]2/2/RGE\[CapitalLambda]6*RGEGet\[Kappa][RGE\[Theta]12, RGE\[Theta]13, RGE\[Theta]23, RGE\[Delta], RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2, RGEMlightest, RGE\[CapitalDelta]m2atm, RGE\[CapitalDelta]m2sol, RGEMassHierarchy, RGEvu],
	RGE\[Kappa] -> 0*IdentityMatrix[3], 
	RGE\[Lambda]->0.5,
	RGE\[CapitalLambda]1->0.5,
	RGE\[CapitalLambda]2->0.5,
	RGE\[CapitalLambda]4->0.5,
	RGE\[CapitalLambda]5->0.5,
	RGE\[CapitalLambda]6->10^9,
	RGEM\[CapitalDelta]2->10^18
}
},
{"MZ",{
	RGE\[Kappa] -> 0*IdentityMatrix[3], 
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
        RGEY\[CapitalDelta]->-RGEM\[CapitalDelta]2/2/RGE\[CapitalLambda]6*RGEGet\[Kappa][RGE\[Theta]12, RGE\[Theta]13, RGE\[Theta]23, RGE\[Delta], RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2, RGEMlightest, RGE\[CapitalDelta]m2atm, RGE\[CapitalDelta]m2sol, RGEMassHierarchy, RGEvu],
	RGEg1 -> RGEgMZ[1],
	RGEg2 -> RGEgMZ[2],
	RGEg3 -> RGEgMZ[3],
	RGE\[Lambda] -> 0.5,
	RGE\[CapitalLambda]1->0.5,
	RGE\[CapitalLambda]2->0.5,
	RGE\[CapitalLambda]4->0.5,
	RGE\[CapitalLambda]5->0.5,
	RGE\[CapitalLambda]6->10^9,
	RGEM\[CapitalDelta]2->10^18,
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
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[Y\[CapitalDelta][t],t]==BetaY\[CapitalDelta][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[Lambda][t],t]==Beta\[Lambda][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[CapitalLambda]1[t],t]==Beta\[CapitalLambda]1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[CapitalLambda]2[t],t]==Beta\[CapitalLambda]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[CapitalLambda]4[t],t]==Beta\[CapitalLambda]4[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[CapitalLambda]5[t],t]==Beta\[CapitalLambda]5[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[\[CapitalLambda]6[t],t]==Beta\[CapitalLambda]6[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]],
		D[M\[CapitalDelta]2[t],t]==BetaM\[CapitalDelta]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[Lambda][t],\[CapitalLambda]1[t],\[CapitalLambda]2[t],\[CapitalLambda]4[t],\[CapitalLambda]5[t],\[CapitalLambda]6[t],M\[CapitalDelta]2[t]]}; (* These are the Renormalization group equations *)

              
(* Beta Functions of the Standardmodel *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[CapitalDelta], Beta\[Kappa], Beta\[Lambda], Beta\[CapitalLambda]1, Beta\[CapitalLambda]2, Beta\[CapitalLambda]4, Beta\[CapitalLambda]5, Beta\[CapitalLambda]6, BetaM\[CapitalDelta]2];

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=
	47/10 * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=
	-5/2 * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=
	-7 * 1/(16*Pi^2) * g3^3;

BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
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
        )*Yd
	);

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=1/(16*Pi^2) * (
       (3/2)*Yu.(
       - Dagger[Yd].Yd
       + Dagger[Yu].Yu
       )
       + (
       - (17/20)*g1^2
       - (9/4)*g2^2
       - 8*g3^2
       + Tr[Dagger[Ye].Ye]
       + 3*Tr[Dagger[Yu].Yu]
       + 3*Tr[Dagger[Yd].Yd]
       )*Yu
       );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=1/(16*Pi^2) * (
	(3/2)*Ye.(
        Dagger[Ye].Ye
	+ Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
        )
	+ (
	- (9/4)*g1^2
	- (9/4)*g2^2
	+ 3*Tr[Dagger[ Yd ].Yd ]
	+ 3*Tr[Dagger[ Yu ].Yu ]
	+ Tr[Dagger[ Ye ].Ye ]
        )*Ye
	);


BetaY\[CapitalDelta][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=1/(16*Pi^2) * (
	Y\[CapitalDelta].(
	+ (1/2) * Dagger[Ye].Ye
	+ (3/2) * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
        )
	+ (
	+ (1/2) * Transpose[Ye].Conjugate[Ye]
	- (3/2) *Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
        ).Y\[CapitalDelta]
	+ (
	- (9/10)*g1^2
	- (9/2)*g2^2
	+ Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
        )*Y\[CapitalDelta]
	);


Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] :=1/(16*Pi^2) * (
	\[Kappa].(
        (-3/2)*Dagger[Ye].Ye
	+ (3/2)*Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
        )
        + (
	- (3/2)*Transpose[Ye].Conjugate[Ye]
	+ (3/2)*Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
        ).\[Kappa]
	+ (
	- 3*g2^2
	+ 6*Tr[Dagger[Yu].Yu]
	+ 6*Tr[Dagger[Yd].Yd]
	+ 2*Tr[Dagger[Ye].Ye]
	+ \[Lambda]
        )*\[Kappa]
	);

Beta\[Lambda][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
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
    + 12 * \[CapitalLambda]4^2
    + 8 * \[CapitalLambda]5^2
    - 8*(
    + 3*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
    + 3*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
    + Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
    )
    );

Beta\[CapitalLambda]1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    \[CapitalLambda]1 *(
    (-36/5)*g1^2
    - 24 * g2^2
    + 14 * \[CapitalLambda]1
    + 4 * \[CapitalLambda]2
    + 4 * Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
    )
    +(108/25) * g1^4
    + 18 * g2^4
    + (72/5) * g1^2 * g2^2
    + 2 * \[CapitalLambda]2^2
    + 4 * \[CapitalLambda]4^2
    + 4 * \[CapitalLambda]5^2
    - 8 * Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta].Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]    
    );

Beta\[CapitalLambda]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    \[CapitalLambda]2 *(
    (-36/5)*g1^2
    - 24 * g2^2
    + 12 *\[CapitalLambda]1
    + 3 * \[CapitalLambda]2
    + 4 * Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
    )
    + 12 * g2^4
    - (144/5) * g1^2 * g2^2
    - 8 * \[CapitalLambda]5^2
    + 8 * Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta].Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]    
    );

Beta\[CapitalLambda]4[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    \[CapitalLambda]4 *(
    (-9/2)*g1^2
    - (33/2) * g2^2
    + 8 *\[CapitalLambda]1
    + 2 * \[CapitalLambda]2
    + 3 * \[Lambda]
    + 4 * \[CapitalLambda]4
    +2 * (
	+ Tr[Dagger[Ye].Ye]
	+ 3*Tr[Dagger[Yu].Yu]
	+ 3*Tr[Dagger[Yd].Yd]
        + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
        )
    )
    + 6 * g2^4
    + (27/25) * g1^4
    + 8 * \[CapitalLambda]5^2
    );

Beta\[CapitalLambda]5[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    \[CapitalLambda]5 *(
    (-9/2)*g1^2
    - (33/2) * g2^2
    + 2 *\[CapitalLambda]1
    - 2 * \[CapitalLambda]2
    + \[Lambda]
    + 8 * \[CapitalLambda]4
    +2 * (
	+ Tr[Dagger[Ye].Ye]
	+ 3*Tr[Dagger[Yu].Yu]
	+ 3*Tr[Dagger[Yd].Yd]
        + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
        )
    )
    - (18/5) * g1^2 * g2^2
    );


Beta\[CapitalLambda]6[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    \[CapitalLambda]6 *(
    (-27/10)*g1^2
    - (21/2) * g2^2
    + \[Lambda]
    - 4 * \[CapitalLambda]4
    + 8 * \[CapitalLambda]5
    + Tr[Dagger[Ye].Ye]
    + 3*Tr[Dagger[Yu].Yu]
    + 3*Tr[Dagger[Yd].Yd]
    - Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
    )
    );


BetaM\[CapitalDelta]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[Lambda]_,\[CapitalLambda]1_,\[CapitalLambda]2_,\[CapitalLambda]4_,\[CapitalLambda]5_,\[CapitalLambda]6_,M\[CapitalDelta]2_] := 1/(16*Pi^2) * (
    M\[CapitalDelta]2 *(
    (-18/5)*g1^2
    - 12 * g2^2
    + 8 * \[CapitalLambda]1
    + 2 * \[CapitalLambda]2
    + 2 * Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
    )
    - 4 * \[CapitalLambda]4 * \[Lambda] * RGEvu^2 
    + \[CapitalLambda]6*Conjugate[\[CapitalLambda]6]
    );

(* 4 \[Lambda]4 m^2 using 4 m^2 = -\[Lambda] v^2 *)


ClearAll[TransSMTriplet0N];
TransSMTriplet0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda],RGE\[CapitalLambda]1->l\[CapitalLambda]1,RGE\[CapitalLambda]2->l\[CapitalLambda]2,RGE\[CapitalLambda]4->l\[CapitalLambda]4,RGE\[CapitalLambda]5->l\[CapitalLambda]5,RGE\[CapitalLambda]6->l\[CapitalLambda]6,RGEM\[CapitalDelta]2->lM\[CapitalDelta]2}];
];

ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

        l\[Lambda]+=2 l\[CapitalLambda]6*Conjugate[l\[CapitalLambda]6]/lM\[CapitalDelta]2;
        l\[Kappa]-=2 l\[CapitalLambda]6/lM\[CapitalDelta]2 * lY\[CapitalDelta];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda]}];
];


ClearAll[TransSMTriplet];
TransSMTriplet[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda],RGE\[CapitalLambda]1->l\[CapitalLambda]1,RGE\[CapitalLambda]2->l\[CapitalLambda]2,RGE\[CapitalLambda]4->l\[CapitalLambda]4,RGE\[CapitalLambda]5->l\[CapitalLambda]5,RGE\[CapitalLambda]6->l\[CapitalLambda]6,RGEM\[CapitalDelta]2->lM\[CapitalDelta]2,RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];



ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2},
(* make a transition from the SM to the SM *)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[Lambda],l\[CapitalLambda]1,l\[CapitalLambda]2,l\[CapitalLambda]4,l\[CapitalLambda]5,l\[CapitalLambda]6,lM\[CapitalDelta]2}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

        l\[Lambda]+=2 l\[CapitalLambda]6*Conjugate[l\[CapitalLambda]6]/lM\[CapitalDelta]2;
        l\[Kappa]-=2 l\[CapitalLambda]6/lM\[CapitalDelta]2 * lY\[CapitalDelta];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda],RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];

(* internal functions *)


ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[CapitalDelta][pScale],\[Kappa][pScale],\[Lambda][pScale],\[CapitalLambda]1[pScale],\[CapitalLambda]2[pScale],\[CapitalLambda]4[pScale],\[CapitalLambda]5[pScale],\[CapitalLambda]6[pScale],M\[CapitalDelta]2[pScale]};


ClearAll[SetInitial];
SetInitial[pBoundary_?NumericQ,pInitial_]:=Block[{},
(* sets the initial values *)
   Return[		{g1[pBoundary]==RGEg1,
			g2[pBoundary]==RGEg2,
			g3[pBoundary]==RGEg3,
			Yu[pBoundary]==RGEYu,
			Yd[pBoundary]==RGEYd,
			Ye[pBoundary]==RGEYe,
			Y\[CapitalDelta][pBoundary]==RGEY\[CapitalDelta],
			\[Kappa][pBoundary]==RGE\[Kappa],
			\[Lambda][pBoundary]==RGE\[Lambda],
			\[CapitalLambda]1[pBoundary]==RGE\[CapitalLambda]1,
			\[CapitalLambda]2[pBoundary]==RGE\[CapitalLambda]2,
			\[CapitalLambda]4[pBoundary]==RGE\[CapitalLambda]4,
			\[CapitalLambda]5[pBoundary]==RGE\[CapitalLambda]5,
			\[CapitalLambda]6[pBoundary]==RGE\[CapitalLambda]6,
			M\[CapitalDelta]2[pBoundary]==RGEM\[CapitalDelta]2
			}//.pInitial
			];

];
End[]; (* end of `Private` *)


EndPackage[];
