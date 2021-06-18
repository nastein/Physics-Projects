(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEMSSMTriplet0N`",{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGESMTriplet0N`", "REAP`RGESM0N`","REAP`RGE2HDM0N`" (* transtions to these models are possible *)
}];


(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! package depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
(* dependencies are marked by !!! in the file *)

(* register MSSMTriplet0N *)
RGERegisterModel["MSSMTriplet0N","REAP`RGEMSSMTriplet0N`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGEMd->`Private`GetMd,RGEPoleMTop->`Private`GetPoleMTop,RGE\[Kappa]->`Private`Get\[Kappa],RGERawY\[CapitalDelta]->`Private`GetRawY\[CapitalDelta],RGEYu->`Private`GetYu,RGETwistingParameters->`Private`GetTwistingParameters,RGECoupling->`Private`GetCoupling,RGEMixingParameters->`Private`GetMixingParameters,RGEM\[Nu]->`Private`GetM\[Nu],RGEYd->`Private`GetYd,RGEYe->`Private`GetYe,RGERaw->`Private`GetRawSolution,RGEAll->`Private`GetSolution,RGEMu->`Private`GetMu,RGEM\[CapitalDelta]->`Private`GetM\[CapitalDelta],RGEM\[Nu]r->`Private`GetM\[Nu]r,RGE\[Alpha]->`Private`Get\[Alpha],RGEY\[Nu]->`Private`GetY\[Nu],RGEMe->`Private`GetMe},
{{"SM",`Private`TransSM},{"SM0N",`Private`TransSM0N},{"MSSM0N",`Private`TransMSSM0N},{"2HDM",`Private`Trans2HDM},{"MSSMTriplet",`Private`TransMSSMTriplet},{"2HDM0N",`Private`Trans2HDM0N},{"MSSM",`Private`TransMSSM},{"MSSMTriplet0N",`Private`TransMSSMTriplet0N},{"SMTriplet",`Private`TransSMTriplet},{"SMTriplet0N",`Private`TransSMTriplet0N}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`", "REAP`RGESMTriplet0N`", "REAP`RGESM0N`","REAP`RGE2HDM0N`"}];

ModelName="MSSMTriplet0N";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of MSSMTriplet0N *)
(* GetRawY\[Nu] is not a function of MSSMTriplet0N *)
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

(* Get\[Lambda] is not a function of MSSMTriplet0N *)
(* GetM\[CapitalDelta]2 is not a function of MSSMTriplet0N *)
ClearAll[GetM\[CapitalDelta]];
GetM\[CapitalDelta][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the coupling constants *)
   Return[(M\[CapitalDelta][pScale]/.pSolution)[[1]]];
];

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

(* Get\[Kappa]1 is not a function of MSSMTriplet0N *)
(* Get\[Kappa]2 is not a function of MSSMTriplet0N *)
ClearAll[GetM\[Nu]r];
GetM\[Nu]r[pScale_,pSolution_,pOpts___]:=Block[{lMhigh,lCutoff},
(* returns the mass matrix of the heavy neutrinos *)
            lCutoff=RGEGetCutoff[Exp[pScale],1];
            lMhigh=RGEGetSolution[lCutoff,RGEM\[Nu]r,1];
        Return[lMhigh];
];

ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,lM\[CapitalDelta],lvu},
(* returns the mass matrix of the neutrinos *)
	{l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,lM\[CapitalDelta]}=({\[Kappa][pScale],Y\[CapitalDelta][pScale],\[CapitalLambda]u[pScale],M\[CapitalDelta][pScale]}/.pSolution)[[1]];
	lOpts;
	Options[lOpts]=Options[RGEOptions];
	SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
	lvu=RGEv\[Nu]/.Options[lOpts];

	Return[-lvu^2*(1/2)^2*(10^9*l\[Kappa]-2 *l\[CapitalLambda]u/(lM\[CapitalDelta]/10^9)*lY\[CapitalDelta])];
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
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lM,lY\[Nu],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2},
(* returns all parameters of the SM *)
        {lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2}=(ParametersFunc[pScale]/.pSolution)[[1]];
	lM=GetM\[Nu]r[pScale,pSolution,pOpts];
	lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
        Return[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2}];
];

(* GetM1Tilde is not a function of MSSMTriplet0N *)
(* Get\[Epsilon]1Max is not a function of MSSMTriplet0N *)
(* Get\[Epsilon]1 is not a function of MSSMTriplet0N *)
(* GetGWCond is not a function of MSSMTriplet0N *)
(* GetGWConditions is not a function of MSSMTriplet0N *)
(* GetVEVratio is not a function of MSSMTriplet0N *)
(* GetVEVratios is not a function of MSSMTriplet0N *)
ClearAll[ModelName0N];
ModelName0N=ModelName<>"0N";

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


Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[CapitalDelta],\[Kappa],\[CapitalLambda]u,\[CapitalLambda]d,M\[CapitalDelta]};
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[CapitalDelta],RGE\[Kappa],RGE\[CapitalLambda]u,RGE\[CapitalLambda]d,RGEM\[CapitalDelta]};


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
	RGE\[Kappa]->0 * IdentityMatrix[3],
        RGE\[CapitalLambda]u->0.5,
        RGE\[CapitalLambda]d->0.5,
        RGEM\[CapitalDelta]->10^9,
	RGEY\[CapitalDelta]->-RGEM\[CapitalDelta]/2/RGE\[CapitalLambda]u*RGEGet\[Kappa][RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu]
}
}
}; (* a list containing the suggestions of initial values *)


ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Y\[CapitalDelta][t],t]==BetaY\[CapitalDelta][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[CapitalLambda]u[t],t]==Beta\[CapitalLambda]u[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[CapitalLambda]d[t],t]==Beta\[CapitalLambda]d[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[M\[CapitalDelta][t],t]==BetaM\[CapitalDelta][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[CapitalDelta][t],\[Kappa][t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]]
}; (* renormalization group equations of the MSSM ( 1 Loop ) *)


(* Beta functions of the MSSM *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[CapitalDelta], Beta\[Kappa], Beta\[CapitalLambda]u, Beta\[CapitalLambda]d, BetaM\[CapitalDelta]];

(* 1 loop contributions *)


Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(51/5) * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(5) * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(-3) * 1/(16*Pi^2) * g3^3;


BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Yd.(
          + 3*Dagger[Yd].Yd
	  + Dagger[Yu].Yu
          )
          + (
          - (7/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + 3*Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + 3*Tr[Dagger[Yd].Yd]
	  + Tr[Dagger[Ye].Ye]
          )*Yd
          );

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Yu.(
          + Dagger[Yd].Yd
	  + 3*Dagger[Yu].Yu
          )
          + (
          - (13/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + 3*Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + 3*Tr[Dagger[Yu].Yu]
          )*Yu
          );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Ye.(
          + 3*Dagger[Ye].Ye
	  + 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
          )
          + (
          - (9/5)*g1^2
	  - 3*g2^2
	  + 3*Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + 3*Tr[Dagger[Yd].Yd]
	  + Tr[Dagger[Ye].Ye]
          )*Ye
          );

BetaY\[CapitalDelta][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Y\[CapitalDelta].(
          + Dagger[Ye].Ye
	  + 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
          )
          + (
          + Transpose[Ye].Conjugate[Ye]
	  + 3 * Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
          ).Y\[CapitalDelta]
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          )*Y\[CapitalDelta]
          );



Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=1/(16*Pi^2) * (
	\[Kappa].(
        + Dagger[Ye].Ye
	+ 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
        )
	+ (
        + Transpose[Ye].Conjugate[Ye]
        + 3 * Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
        ).\[Kappa]
        + (
	- (6/5)*g1^2
	- 6*g2^2
        + 6 * Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	+ 6*Tr[Dagger[Yu].Yu]
        )*\[Kappa]
        );


Beta\[CapitalLambda]u[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + 7*Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + 6 * Tr[Dagger[Yu].Yu]
          ) * \[CapitalLambda]u
          );


Beta\[CapitalLambda]d[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + 7*Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + 6 * Tr[Dagger[Yd].Yd]
	  + 2 * Tr[Dagger[Ye].Ye]
	  + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          ) * \[CapitalLambda]d
          );

BetaM\[CapitalDelta][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[CapitalDelta]_,\[Kappa]_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
	  + (
          + Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          - (12/5)*g1^2
	  - 8*g2^2
          ) * M\[CapitalDelta]
          );



	(* transition functions *)

ClearAll[TransMSSM0N];
TransMSSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa]}];
];

ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGEY\[Nu]->{},RGEM\[Nu]r->{}}];
];





ClearAll[TransMSSMTriplet0N];
TransMSSMTriplet0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGE\[CapitalLambda]u->l\[CapitalLambda]u,RGE\[CapitalLambda]d->l\[CapitalLambda]d,RGEM\[CapitalDelta]->lM\[CapitalDelta]}];
];

ClearAll[TransMSSMTriplet];
TransMSSMTriplet[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o \[Nu]*)
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGEY\[Nu]->{},RGEM\[Nu]r->{},RGE\[CapitalLambda]u->l\[CapitalLambda]u,RGE\[CapitalLambda]d->l\[CapitalLambda]d,RGEM\[CapitalDelta]->lM\[CapitalLambda]}];
];





ClearAll[Trans2HDM0N];
Trans2HDM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa],lsb,lcb,lTosb,lTocb,lu,le,ld,l\[Nu],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

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
(* changed transition function to 2HDM such that D-terms are correctly taken into account, although there is less flexibility *)
(*,RGE\[Lambda]1->(RGE\[Lambda]1/.Options[lToOpts,RGE\[Lambda]1]),RGE\[Lambda]2->(RGE\[Lambda]2/.Options[lToOpts,RGE\[Lambda]2]),RGE\[Lambda]3->(RGE\[Lambda]3/.Options[lToOpts,RGE\[Lambda]3]),RGE\[Lambda]4->(RGE\[Lambda]4/.Options[lToOpts,RGE\[Lambda]4]),RGE\[Lambda]5->(RGE\[Lambda]5/.Options[lToOpts,RGE\[Lambda]5]) *)

ClearAll[Trans2HDM];
Trans2HDM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lTosb,lTocb,lu,le,ld,l\[Nu],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
	lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=lsb/((RGEzu/.Options[lToOpts,RGEzu]).{lTocb,lTosb});
	ld=lcb/((RGEzd/.Options[lToOpts,RGEzd]).{lTocb,lTosb});
	l\[Nu]=lsb/((RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]]).{lTocb,lTosb});
	le=lcb/lTocb;

        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
	
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,
	RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,
	RGE\[Kappa]1->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[1]]*l\[Kappa]*l\[Nu]^2,
	RGE\[Kappa]2->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[2]]*l\[Kappa]*l\[Nu]^2,
	RGEM\[Nu]r->{},RGE\[Lambda]1->(lg2^2+lg1^2)/2.,RGE\[Lambda]2->(lg2^2+lg1^2)/2.,RGE\[Lambda]3->(lg2^2-lg1^2)/4,RGE\[Lambda]4->-(lg2^2)/2,RGE\[Lambda]5->0.}];
];
(* changed transition function to 2HDM such that D-terms are correctly taken into account, although there is less flexibility *)
(*,RGE\[Lambda]1->(RGE\[Lambda]1/.Options[lToOpts,RGE\[Lambda]1]),RGE\[Lambda]2->(RGE\[Lambda]2/.Options[lToOpts,RGE\[Lambda]2]),RGE\[Lambda]3->(RGE\[Lambda]3/.Options[lToOpts,RGE\[Lambda]3]),RGE\[Lambda]4->(RGE\[Lambda]4/.Options[lToOpts,RGE\[Lambda]4]),RGE\[Lambda]5->(RGE\[Lambda]5/.Options[lToOpts,RGE\[Lambda]5]) *)


ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];

        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}]; (* !!! Parameters !!!*)
];


ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->{},RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}]; (* !!! Parameters !!!*)
];




ClearAll[TransSMTriplet0N];
TransSMTriplet0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM w/o \[Nu] *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
        Print["The transition from MSSMTriplet0N to SMTriplet0N is not fully implemented. The functions to define the parameters in the Higgs potential are not fully implemented."];

        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[CapitalDelta]->lY\[CapitalDelta]*(lsb)^2,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]]),RGE\[CapitalLambda]6->l\[CapitalLambda]u*lM\[CapitalDelta],RGEM\[CapitalDelta]2->lM\[CapitalDelta]^2}]; (* !!! Parameters !!!*)
];


ClearAll[TransSMTriplet];
TransSMTriplet[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,l\[Kappa], lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM w/o \[Nu] *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)

(* exceptions: try to add new particles --> CanNotAddNewParticles
*)
        Print["The transition from MSSMTriplet0N to SMTriplet is not fully implemented. The functions to define the parameters in the Higgs potential are not fully implemented."];

        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        (* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[CapitalDelta],l\[Kappa],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];

	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[CapitalDelta]->lY\[CapitalDelta]*(lsb)^2,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->{},RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]]),RGE\[CapitalLambda]6->l\[CapitalLambda]u*lM\[CapitalDelta],RGEM\[CapitalDelta]2->lM\[CapitalDelta]^2}]; (* !!! Parameters !!!*)
];




(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[CapitalDelta][pScale],\[Kappa][pScale],\[CapitalLambda]u[pScale],\[CapitalLambda]d[pScale],M\[CapitalDelta][pScale]};

ClearAll[SetIntial];
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
			\[CapitalLambda]u[pBoundary]==RGE\[CapitalLambda]u,
			\[CapitalLambda]d[pBoundary]==RGE\[CapitalLambda]d,
			M\[CapitalDelta][pBoundary]==RGEM\[CapitalDelta]
			}//.pInitial
			];
];


End[]; (* end of `Private`*)


EndPackage[];
