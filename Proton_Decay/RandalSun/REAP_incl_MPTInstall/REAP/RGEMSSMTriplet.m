(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)





BeginPackage["REAP`RGEMSSMTriplet`",{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGEMSSM`","REAP`RGESM`","REAP`RGEMSSM0N`","REAP`RGESM0N`","REAP`RGE2HDM`","REAP`RGE2HDM0N`",
"REAP`RGESMTriplet`","REAP`RGEMSSMTriplet0N`","REAP`RGESMTriplet0N`"(* transtions to these models are possible *)
}];

(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! package depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
(* dependencies are marked by !!! in the file *)


(* register MSSMTriplet *)
RGERegisterModel["MSSMTriplet","REAP`RGEMSSMTriplet`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGE\[Epsilon]1Max->`Private`Get\[Epsilon]1Max,RGE\[Kappa]->`Private`Get\[Kappa],RGEMu->`Private`GetMu,RGEMe->`Private`GetMe,RGEM\[Nu]r->`Private`GetM\[Nu]r,RGERawY\[CapitalDelta]->`Private`GetRawY\[CapitalDelta],RGECoupling->`Private`GetCoupling,RGEYd->`Private`GetYd,RGEMixingParameters->`Private`GetMixingParameters,RGEPoleMTop->`Private`GetPoleMTop,RGEM1Tilde->`Private`GetM1Tilde,RGE\[Alpha]->`Private`Get\[Alpha],RGERawM\[Nu]r->`Private`GetRawM\[Nu]r,RGEYe->`Private`GetYe,RGEYu->`Private`GetYu,RGE\[Epsilon]->`Private`Get\[Epsilon],RGE\[Epsilon]1->`Private`Get\[Epsilon]1,RGEM\[Nu]->`Private`GetM\[Nu],RGETwistingParameters->`Private`GetTwistingParameters,RGERaw->`Private`GetRawSolution,RGEAll->`Private`GetSolution,RGE\[Epsilon]Max->`Private`Get\[Epsilon]Max,RGEMd->`Private`GetMd,RGEY\[Nu]->`Private`GetY\[Nu],RGEM\[CapitalDelta]->`Private`GetM\[CapitalDelta],RGERawY\[Nu]->`Private`GetRawY\[Nu]},
{{"SMTriplet",`Private`TransSMTriplet},{"SM0NTriplet",`Private`TransSMTriplet0N},{"SM0N",`Private`TransSM0N},{"MSSMTriplet0N",`Private`TransMSSMTriplet0N},{"SM",`Private`TransSM},{"2HDM",`Private`Trans2HDM},{"MSSM",`Private`TransMSSM},{"MSSMTriplet",`Private`TransMSSMTriplet},{"MSSM0N",`Private`TransMSSM0N},{"2HDM0N",`Private`Trans2HDM0N}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGEMSSM`","REAP`RGESM`","REAP`RGEMSSM0N`","REAP`RGESM0N`","REAP`RGE2HDM`","REAP`RGE2HDM0N`", "REAP`RGESMTriplet`","REAP`RGEMSSMTriplet0N`","REAP`RGESMTriplet0N`"}];

ModelName="MSSMTriplet";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

ClearAll[GetRawM\[Nu]r];
GetRawM\[Nu]r[pScale_,pSolution_,pOpts___]:=Block[{lM},
(* returns the mass matrix of the heavy neutrinos *)
	lM=(M\[Nu]r[pScale]/.pSolution)[[1]];
        Return[lM];
];

ClearAll[GetRawY\[Nu]];
GetRawY\[Nu][pScale_,pSolution_,pOpts___]:=Block[{},
(* returns the Yukawa coupling matrix of the neutrinos *)
    Return[(Y\[Nu][pScale]/.pSolution)[[1]]];
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

(* Get\[Lambda] is not a function of MSSMTriplet *)
(* GetM\[CapitalDelta]2 is not a function of MSSMTriplet *)
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

(* Get\[Kappa]1 is not a function of MSSMTriplet *)
(* Get\[Kappa]2 is not a function of MSSMTriplet *)
ClearAll[GetM\[Nu]r];
GetM\[Nu]r[pScale_,pSolution_,pOpts___]:=Block[{lM,lIntegratedOut,lLenM,lMhigh,lf,lg,lMEval,lCutoff,lUp},
(* returns the mass matrix of the heavy neutrinos *)
        lOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lIntegratedOut=(RGEIntegratedOut/.Options[lOpts]);
	lM=(M\[Nu]r[pScale]/.pSolution)[[1]];
        Catch[
          If[lIntegratedOut>0,
            lCutoff=RGEGetCutoff[Exp[pScale],1];
	    lUp=RGEGetRange[][[2]];
	    If[lCutoff>Exp[lUp],
	    	lMInitial=(RGEM\[Nu]r/.RGEGetInitial[][[2]]);
		lLenM=Length[lMInitial];
		If[lLenM>Length[lM],
	            lMEval=Sqrt[Sort[Abs[Eigenvalues[Dagger[lMInitial].lMInitial]],Greater]];
		    lM=Table[
				If[lf==lg, If[lf<=lLenM-lIntegratedOut,lM[[lf,lg]],lMEval[[lLenM-lf+1]] ],
					   If[(lf<=lLenM-lIntegratedOut)&&(lg<=lLenM-lIntegratedOut),lM[[lf,lg]],0]
				],
				{lf,lLenM},{lg,lLenM}
			];
		];
		Throw[lCutoff,RGEScaleTooBig]
	    ];
            lMhigh=RGEGetSolution[lCutoff,RGEM\[Nu]r,1];
            lLenM=Length[lMhigh];
            lMEval=Sqrt[Sort[Abs[Eigenvalues[Dagger[lMhigh].lMhigh]],Greater]];
            lM=Table[
              If[lf==lg,If[lf<=lLenM-lIntegratedOut,lM[[lf,lg]],lMEval[[lLenM-lf+1]] ],
                        If[(lf<=lLenM-lIntegratedOut)&&(lg<=lLenM-lIntegratedOut),lM[[lf,lg]],0]
              ],
              {lf,lLenM},{lg,lLenM}];
            ];
        ,RGEScaleTooBig];
        Return[lM];
];

ClearAll[GetM\[Nu]];
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa],lY\[Nu],lM,lY\[CapitalDelta],l\[CapitalLambda]u,lM\[CapitalDelta],lvu},
(* returns the mass matrix of the neutrinos *)
	{l\[Kappa],lY\[Nu],lM,lY\[CapitalDelta],l\[CapitalLambda]u,lM\[CapitalDelta]}=({\[Kappa][pScale],Y\[Nu][pScale],M\[Nu]r[pScale],Y\[CapitalDelta][pScale],\[CapitalLambda]u[pScale],M\[CapitalDelta][pScale]}/.pSolution)[[1]];
	lOpts;
	Options[lOpts]=Options[RGEOptions];
	SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lvu=RGEv\[Nu]/.Options[lOpts];
	If[MatrixConditionNumber[lM]>2*Precision[lM],Print["GetM\[Nu]: The
	matrix M=", MatrixForm[lM]], " is ill-conditioned and the condition
	number is ", MatrixConditionNumber[lM]];

	Return[-lvu^2*(1/2)^2*(10^9*l\[Kappa]+2*Transpose[lY\[Nu]].(Inverse[lM]*10^9).lY\[Nu]-2 *l\[CapitalLambda]u/(lM\[CapitalDelta]/10^9)*lY\[CapitalDelta])];
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
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lM,lY\[Nu],l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2},
(* returns all parameters of the SM *)
        {lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2}=(ParametersFunc[pScale]/.pSolution)[[1]];
	lM=GetM\[Nu]r[pScale,pSolution,pOpts];
	lY\[Nu]=GetY\[Nu][pScale,pSolution,pOpts];
        Return[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]2}];
];

ClearAll[GetM1Tilde];
GetM1Tilde[pScale_,pSolution_,pOpts___]:=Block[{lv,lM,lM\[Nu],lm1},
(* returns m1 Tilde *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lv=(RGEv\[Nu]/Sqrt[2])/.Options[lOpts];
   lM=(M\[Nu]r[pScale]/.pSolution)[[1]];
   lM\[Nu]=lv*(Y\[Nu][pScale]/.pSolution)[[1]];
   {lM,lM\[Nu]}=RGERotateM[lM,lM\[Nu]];
   lM\[Nu]=lM\[Nu][[1]];
   lm1=(lM\[Nu].Conjugate[lM\[Nu]])/lM[[1,1]];
   Return[RGEFloor[Re[lm1]][[1]]*10^9];
];

(* GetGWCond is not a function of MSSMTriplet *)
(* GetGWConditions is not a function of MSSMTriplet *)
(* GetVEVratio is not a function of MSSMTriplet *)
(* GetVEVratios is not a function of MSSMTriplet *)
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
SolveModel[{pUp_?NumericQ,pUpModel_,pUpOptions_},{pDown_?NumericQ,pDownModel_,pDownOptions_},pDirection_?NumericQ,pBoundary_?NumericQ,pInitial_,pNDSolveOpts_,pOpts___]:=Block[{lSolution, (* contains the solution *)
lInitial, (* initial values in the format needed by NDSolve *)
lODE, (* RGE and lInitial *)
lAddedModels, (* number of added models *)
lIntegratedOut, (* # of neutrinos integrated out *)
lAlreadyIntegratedOut, (* # of neutrinos already inetgrated out *)
lTransition, (* transition to be added *)
lIndexModel, (* variant of the model *)
lM,lY,l\[Kappa] (* needed for check of too high masses *)
},
(* solves the model; returns the solution, the new scale and the number of added models *)
(* exceptions: too many iterations --> RGETooManyIterations
               \[Nu] mass above cutoff --> RGE\[Nu]MassAboveCutoff 
*)
(* solve differential equation *)
	lOpts;
	lAddModelOpts;
        Options[lOpts]=Options[RGEOptions];
        SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lLogThresholdFactor=Log[(RGEThresholdFactor/.Options[lOpts,RGEThresholdFactor])];
	lInitial;
	Options[lInitial]=pInitial;
	If[Length[(RGEM\[Nu]r/.pInitial)]>0,
		{lM,lY\[Nu],l\[Kappa]}=({RGEM\[Nu]r,RGEY\[Nu],RGE\[Kappa]}/.pInitial);
		If[(RGESearchTransition/.Options[lOpts,RGESearchTransition]), (*CHECK THIS LINE: inverted condition*)
			{lM,lY,l\[Kappa],lIntegratedOut}=RGETestMY\[Nu][Exp[lLogThresholdFactor+pBoundary],lM,lY\[Nu],l\[Kappa]];
			If[lIntegratedOut>0,
				Print["Warning: ",lIntegratedOut, " right-handed neutrinos have a mass above ", Exp[pBoundary], " GeV (",Select[Sqrt[Abs[Eigenvalues[(Dagger[RGEM\[Nu]r].RGEM\[Nu]r/.pInitial)]]],(#1>Exp[pBoundary])&],"). Thus they have been integrated out."];
				SetOptions[lOpts,RGEIntegratedOut->(((RGEIntegratedOut/.Options[lOpts,RGEIntegratedOut])/.{RGEIntegratedOut->0})+lIntegratedOut)];

				SetOptions[lInitial,RGEM\[Nu]r->lM,RGEY\[Nu]->lY,RGE\[Kappa]->l\[Kappa]];
				RGEChangeOptions[Exp[pBoundary],RGEIntegratedOut->(RGEIntegratedOut/.Options[lOpts,RGEIntegratedOut])];
				];
			];
	];
	If[Length[RGEM\[Nu]r/.Options[lInitial]]==0,
		lTransition=RGEGetCutoff[(Exp[pUp]+Exp[pDown])/2];
		RGEDelModel[(Exp[pUp]+Exp[pDown])/2];
		Options[lAddModelOpts]=RGEFilterOptions[RGEOptions,pOpts];
		SetOptions[lAddModelOpts,RGEFilterOptions[lAddModelOpts,RGEAutoGenerated->ModelName]];
		Options[lAddModelOpts]=Options[lAddModelOpts]~Union~{RGEAutoGenerated->ModelName};
		RGEAddEFT[ModelName0N, RGECutoff->lTransition,RGEFilterOptions[RGEOptions,Options[lAddModelOpts]]];
                RGESetModelOptions[ModelName0N,Options[RGEOptions]];
		{lSolution,lTransition,lAddedModels}=(RGEGetSolveModel[ModelName0N])[{pUp,pUpModel,pUpOptions},{pDown,pDownModel,pDownOptions},pDirection,pBoundary,Options[lInitial],pNDSolveOpts,Sequence[pOpts]];
		Return[{lSolution,lTransition,lAddedModels}];
	];
		
	
	lNDSolveOpts;
	Options[lNDSolveOpts]=Options[NDSolve];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,Options[RGEOptions]]];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,pOpts]];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,Sequence[pNDSolveOpts]]];
	
	lInitial=SetInitial[pBoundary,Options[lInitial]];
	lIndexModel=Flatten[ Position[ModelVariants,(RGEModelVariant/.Options[lOpts])] ][[ 1 ]];
	lODE=RGE[[lIndexModel]]/.Options[lOpts];

        lSolution=NDSolve[lODE ~Join~ lInitial, Parameters,{t,pDown,pUp}, Sequence[Options[lNDSolveOpts]]];

     
(* search transitions *)
        lAddedModels=0;
        lAlreadyIntegratedOut=RGEIntegratedOut/.Options[lOpts,RGEIntegratedOut] /. {RGEIntegratedOut->0};
        lPrecision=RGEPrecision/.Options[lOpts,RGEPrecision];
        lIntegratedOut=0;
(* check whether transitions should not be searched for *)
        If[!(RGESearchTransition/.Options[lOpts,RGESearchTransition]),Return[{lSolution,pDown,lAddedModels}]];
        lDownOpts;
        If[pDownModel!="",
            Options[lDownOpts]=RGEGetModelOptions[pDownModel][[1,2]],
	    Options[lDownOpts]={RGEIntegratedOut->0}
        ];
        SetOptions[lDownOpts,RGEFilterOptions[lDownOpts,pDownOptions]];
(* check whether there can not be any transitions *)
        lLengthM\[Nu]=Length[(M\[Nu]r[pUp]/.lSolution)[[1]] ];

(* to do *)
(*        If[(pDownModel==ModelName0N),If[(lLengthM\[Nu]-(RGEIntegratedOut/.Options[lOpts,RGEIntegratedOut])<=1),Return[{lSolution,pDown,lAddedModels}]]];
        If[(pDownModel==ModelName),If[(RGEIntegratedOut/.Options[lDownOpts,RGEIntegratedOut])-(RGEIntegratedOut/.Options[lOpts,RGEIntegratedOut])<=1,Return[{lSolution,pDown,lAddedModels}]]];
*)
         lMTransitions=RGESearchTransitions[(GetRawM\[Nu]r[#1,lSolution])&, pUp,pUp,pDown,Sequence[Options[lOpts]]];
          
         If[Length[lMTransitions]>0,
             lThreshold=N[lLogThresholdFactor+Max[lMTransitions],lPrecision];
             lMTransitions=N[(Select[lMTransitions,(#1>=lThreshold)&])+lLogThresholdFactor,lPrecision];
             While[Length[lMTransitions]>0,
                 lTransition=First[lMTransitions];
                 lDegeneracy=Length[Select[lMTransitions,(Chop[#1-lTransition,10^(-lPrecision)]==0)&]];
                 lMTransitions=Sort[Select[lMTransitions ,(Chop[#1-lTransition,10^(-lPrecision)]!=0)&],Greater];
                 lIntegratedOut+=lDegeneracy;		 
		 If[lIntegratedOut<lLengthM\[Nu],
			Options[lAddModelOpts]=RGEFilterOptions[RGEOptions,pOpts];
			SetOptions[lAddModelOpts,RGEFilterOptions[lAddModelOpts,RGEIntegratedOut->lIntegratedOut+lAlreadyIntegratedOut,RGEAutoGenerated->True]];
			Options[lAddModelOpts]=Options[lAddModelOpts]~Union~{RGEIntegratedOut->lIntegratedOut+lAlreadyIntegratedOut,RGEAutoGenerated->True};
			RGEAddEFT[ModelName,RGECutoff->Exp[lTransition],RGEFilterOptions[RGEOptions,Options[lAddModelOpts]]],
                        RGESetModelOptions[ModelName0N,Options[RGEOptions]];
			Options[lAddModelOpts]=RGEFilterOptions[RGEOptions,pOpts];
			SetOptions[lAddModelOpts,RGEFilterOptions[lAddModelOpts,RGEAutoGenerated->True]];
			Options[lAddModelOpts]=Options[lAddModelOpts]~Union~{RGEAutoGenerated->True};
			RGEAddEFT[ModelName0N,RGECutoff->Exp[lTransition],RGEFilterOptions[RGEOptions,Options[lAddModelOpts]]]
		 ];
		 lAddedModels++;
		],
	lThreshold=pDown;
	];
	
	Return[{lSolution,lThreshold,lAddedModels}];    
];


(* definitions for the Minimal Supersymmetric Standard Model (MSSM) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={   RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
			RGEPrecision->6, (* precision to find transitions *)
                        RGEMaxNumberIterations->20, (* maximum number of iterations in the loops to search transitions *)
                        RGEvEW->246, (* vev of the SM Higgs *)
                        RGEtan\[Beta]->50, (* tan \[Beta]=vu/vd *)  
                        RGEIntegratedOut->0, (* number of the integrated out neutrinos *)
			Method->StiffnessSwitching, (* option of NDSolve *)
			RGESearchTransition->True, (* enables/disables the automatic search for transitions *)
                        RGEThresholdFactor->1 (* neutrinos are integrated out at RGEThresholdFactor*Mass *)
			}; (* options of the model *)

                        
Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[Nu],Y\[CapitalDelta],\[Kappa],M\[Nu]r,\[CapitalLambda]u,\[CapitalLambda]d,M\[CapitalDelta]};
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[Nu],RGEY\[CapitalDelta],RGE\[Kappa],RGEM\[Nu]r,RGE\[CapitalLambda]u,RGE\[CapitalLambda]d,RGEM\[CapitalDelta]};

ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.7044110331165641,
	RGEg2->0.6965468498179075,
	RGEg3->0.6983661130877465,
	RGEYd->RGEGetYd[RGEyd,RGEys,RGEyb,RGEq\[Theta]12,RGEq\[Theta]13,RGEq\[Theta]23,RGEq\[Delta],RGEq\[Delta]e,RGEq\[Delta]\[Mu],RGEq\[Delta]\[Tau],RGEq\[CurlyPhi]1,RGEq\[CurlyPhi]2],
	RGEYu->DiagonalMatrix[{RGEyu,RGEyc,RGEyt}],
	RGE\[Kappa]->0*IdentityMatrix[3],
	RGEY\[CapitalDelta]->0*IdentityMatrix[3],
	RGE\[CapitalLambda]u->0.5,
	RGE\[CapitalLambda]d->0.5,
	RGEM\[CapitalDelta]->10^9,
	RGEM\[Nu]r->RGEGetM[RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu,RGEY\[Nu],RGE\[Kappa],RGEY\[CapitalDelta],RGE\[CapitalLambda]u,RGEM\[CapitalDelta]],
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
	RGE\[Theta]12 -> 27 Degree,
	RGE\[Theta]13 -> 0 Degree, 
	RGE\[Theta]23 -> 45 Degree,
	RGE\[Delta] -> 0 Degree,
	RGE\[Delta]e -> 0 Degree,
	RGE\[Delta]\[Mu] -> 0 Degree,
	RGE\[Delta]\[Tau] -> 0 Degree,
	RGE\[CurlyPhi]1 -> 0 Degree,
	RGE\[CurlyPhi]2 -> 0 Degree, 
	RGEMlightest -> 0.02,
	RGE\[CapitalDelta]m2atm -> 4*10^-3, 
	RGE\[CapitalDelta]m2sol -> 1.6 10^-4,
	RGEY\[Nu]33 -> 0.5,
	RGEY\[Nu]Ratio -> 0.1,
	RGEY\[Nu]->RGEGetY\[Nu][RGEY\[Nu]33,RGEY\[Nu]Ratio],
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGEyu -> 0.00104/RGEvu*Sqrt[2],
	RGEyc -> 0.302/RGEvu*Sqrt[2],
	RGEyt -> 129/RGEvu*Sqrt[2],
	RGEyd -> 0.00133/RGEvd*Sqrt[2],
	RGEys -> 0.0265/RGEvd*Sqrt[2],
	RGEyb ->  1.00/RGEvd*Sqrt[2],
	RGEye -> 0.32502032*10^-3*Sqrt[2]/RGEve,
	RGEy\[Mu] -> 68.59813*10^-3*Sqrt[2]/RGEve,
	RGEy\[Tau] -> 1171.4*10^-3*Sqrt[2]/RGEve
}
}
}; (* a list containing suggestions for initial values *)



ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Y\[Nu][t],t]==BetaY\[Nu][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[Y\[CapitalDelta][t],t]==BetaY\[CapitalDelta][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[M\[Nu]r[t],t]==BetaM\[Nu]r[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[CapitalLambda]u[t],t]==Beta\[CapitalLambda]u[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[\[CapitalLambda]d[t],t]==Beta\[CapitalLambda]d[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]],
		D[M\[CapitalDelta][t],t]==BetaM\[CapitalDelta][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],Y\[CapitalDelta][t],\[Kappa][t],M\[Nu]r[t],\[CapitalLambda]u[t],\[CapitalLambda]d[t],M\[CapitalDelta][t]]
}; (* renormalization group equations of the MSSM ( 1 Loop ) *)



(* Beta functions of the MSSM *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[Nu], BetaY\[CapitalDelta], Beta\[Kappa], BetaM\[Nu]r, Beta\[CapitalLambda]u, Beta\[CapitalLambda]d, BetaM\[CapitalDelta]];

(* 1 loop contributions *)

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(51/5) * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(5) * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=
	(-3) * 1/(16*Pi^2) * g3^3;


BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
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

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Yu.(
          + Dagger[Yd].Yd
	  + 3*Dagger[Yu].Yu
          )
          + (
          - (13/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + 3*Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*Tr[Dagger[Yu].Yu]
          )*Yu
          );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Ye.(
          + 3*Dagger[Ye].Ye
	  + Dagger[Y\[Nu]].Y\[Nu]
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

BetaY\[Nu][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Y\[Nu].(
          + Dagger[Ye].Ye
	  + 3*Dagger[Y\[Nu]].Y\[Nu]
	  + 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
          )
          + (
          - (3/5)*g1^2
	  - 3*g2^2
	  + 3*Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*Tr[Dagger[Yu].Yu]
          )*Y\[Nu]
          );

BetaY\[CapitalDelta][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          Y\[CapitalDelta].(
          + Dagger[Ye].Ye
	  + Dagger[Y\[Nu]].Y\[Nu]
	  + 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
          )
          + (
          + Transpose[Ye].Conjugate[Ye]
	  + Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
	  + 3 * Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
          ).Y\[CapitalDelta]
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          )*Y\[CapitalDelta]
          );



Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] :=1/(16*Pi^2) * (
	\[Kappa].(
        + Dagger[Ye].Ye
        + Dagger[Y\[Nu]].Y\[Nu]
	+ 3 * Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]
        )
	+ (
        + Transpose[Ye].Conjugate[Ye]
	+ Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
        + 3 * Transpose[Y\[CapitalDelta]].Conjugate[Y\[CapitalDelta]]
        ).\[Kappa]
        + (
	- (6/5)*g1^2
	- 6*g2^2
        + 6 * Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	+ 6*Tr[Dagger[Yu].Yu]
        + 2*Tr[Dagger[Y\[Nu]].Y\[Nu]]
        )*\[Kappa]
        );

BetaM\[Nu]r[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          2*M\[Nu]r.Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	  + 2*Y\[Nu].Dagger[Y\[Nu]].M\[Nu]r);


Beta\[CapitalLambda]u[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + 7*Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + 6 * Tr[Dagger[Yu].Yu]
	  + 2 * Tr[Dagger[Y\[Nu]].Y\[Nu]]
          ) * \[CapitalLambda]u
          );


Beta\[CapitalLambda]d[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
          + (
          - (9/5)*g1^2
	  - 7*g2^2
	  + 7*Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + 6 * Tr[Dagger[Yd].Yd]
	  + 2 * Tr[Dagger[Ye].Ye]
	  +  Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          ) * \[CapitalLambda]d
          );

BetaM\[CapitalDelta][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,Y\[CapitalDelta]_,\[Kappa]_,M\[Nu]r_,\[CapitalLambda]u_,\[CapitalLambda]d_,M\[CapitalDelta]_] := 1/(16*Pi^2) * (
	  + (
          + Conjugate[\[CapitalLambda]d]*\[CapitalLambda]d
	  + Conjugate[\[CapitalLambda]u]*\[CapitalLambda]u
	  + Tr[Dagger[Y\[CapitalDelta]].Y\[CapitalDelta]]
          - (12/5)*g1^2
	  - 8*g2^2
          ) * M\[CapitalDelta]
          );



	  

(* transition functions *)

ClearAll[TransMSSMTriplet];
TransMSSMTriplet[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut},
(* make a transition from the MSSM to the MSSM *)

(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEOptions];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
	lToIntegratedOut=RGEIntegratedOut/.Options[lToOpts,RGEIntegratedOut];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
	lFromIntegratedOut=RGEIntegratedOut/.Options[lFromOpts,RGEIntegratedOut];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

	If[lFromIntegratedOut>lToIntegratedOut,
	   (* Print["The model ", RGEGetModel[Exp[pScale]], " at the scale ",
	   Exp[pScale]," pretends to have more particles, but new particles can
	   not be added: ",lFromIntegratedOut,"->",lToIntegratedOut, " Thus the
	   model is changed."]; *)
	   RGEChangeOptions[Exp[pScale],RGEIntegratedOut->lFromIntegratedOut];
	   lToIntegratedOut=lFromIntegratedOut;
	   (*Throw[{lFromIntegratedOut,lToIntegratedOut},RGECanNotAddNewParticles];*)
	];
	If[lToIntegratedOut>lFromIntegratedOut,
		lIntegrateOut=lToIntegratedOut-lFromIntegratedOut;
			l\[Kappa]+=RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,lIntegrateOut];
			lY\[Nu]=RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
			lM\[Nu]r=RGEIntegrateOutM[lM\[Nu]rRotated, lIntegrateOut];
	];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[Nu]->lY\[Nu],RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGEM\[Nu]r->lM\[Nu]r,RGE\[CapitalLambda]u->l\[CapitalLambda]u,RGE\[CapitalLambda]d->l\[CapitalLambda]d,RGEM\[CapitalDelta]->lM\[CapitalDelta]}];
];

ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta],lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut},
(* make a transition from the MSSM to the MSSM *)

(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEOptions];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
	lToIntegratedOut=RGEIntegratedOut/.Options[lToOpts,RGEIntegratedOut];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
	lFromIntegratedOut=RGEIntegratedOut/.Options[lFromOpts,RGEIntegratedOut];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

	If[lFromIntegratedOut>lToIntegratedOut,
	   (* Print["The model ", RGEGetModel[Exp[pScale]], " at the scale ",
	   Exp[pScale]," pretends to have more particles, but new particles can
	   not be added: ",lFromIntegratedOut,"->",lToIntegratedOut, " Thus the
	   model is changed."]; *)
	   RGEChangeOptions[Exp[pScale],RGEIntegratedOut->lFromIntegratedOut];
	   lToIntegratedOut=lFromIntegratedOut;
	   (*Throw[{lFromIntegratedOut,lToIntegratedOut},RGECanNotAddNewParticles];*)
	];
	If[lToIntegratedOut>lFromIntegratedOut,
		lIntegrateOut=lToIntegratedOut-lFromIntegratedOut;
			l\[Kappa]+=RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,lIntegrateOut];
			lY\[Nu]=RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
			lM\[Nu]r=RGEIntegrateOutM[lM\[Nu]rRotated, lIntegrateOut];
	];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[Nu]->lY\[Nu],RGE\[Kappa]->l\[Kappa],RGEM\[Nu]r->lM\[Nu]r}];
];





ClearAll[TransMSSMTriplet0N];
TransMSSMTriplet0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o heavy neutrinos *)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[CapitalDelta]->lY\[CapitalDelta],RGE\[Kappa]->l\[Kappa],RGE\[CapitalLambda]u->l\[CapitalLambda]u,RGE\[CapitalLambda]d->l\[CapitalLambda]d,RGEM\[CapitalDelta]->lM\[CapitalDelta]}];
];

ClearAll[TransMSSM0N];
TransMSSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o heavy neutrinos *)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lY\[CapitalDelta],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa]}];
];




ClearAll[Trans2HDM];
Trans2HDM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut,le,lu,ld,l\[Nu],lTosb,lTocb,lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM *)

(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
	lToIntegratedOut=RGEIntegratedOut/.Options[lToOpts,RGEIntegratedOut];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
	lFromIntegratedOut=RGEIntegratedOut/.Options[lFromOpts,RGEIntegratedOut];

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)
        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

	If[lFromIntegratedOut>lToIntegratedOut,
		   (* Print["The model ", RGEGetModel[Exp[pScale]], " at the
		   scale ", Exp[pScale]," pretends to have more particles, but
		   new particles can not be added:
		   ",lFromIntegratedOut,"->",lToIntegratedOut, " Thus the model
		   is changed."]; *)
		   RGEChangeOptions[Exp[pScale],RGEIntegratedOut->lFromIntegratedOut];
		   lToIntegratedOut=lFromIntegratedOut;
		   (*Throw[{lFromIntegratedOut,lToIntegratedOut},RGECanNotAddNewParticles];*)
	   ];
	If[lToIntegratedOut>lFromIntegratedOut,
		lIntegrateOut=lToIntegratedOut-lFromIntegratedOut;
                
			l\[Kappa]+=RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,lIntegrateOut];
			lY\[Nu]=RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
			lM\[Nu]r=RGEIntegrateOutM[lM\[Nu]rRotated, lIntegrateOut];
	];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=lsb/((RGEzu/.Options[lToOpts,RGEzu]).{lTocb,lTosb});
	ld=lcb/((RGEzd/.Options[lToOpts,RGEzd]).{lTocb,lTosb});
	l\[Nu]=lsb/((RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]]).{lTocb,lTosb});
	le=lcb/lTocb;

	
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,RGEY\[Nu]->lY\[Nu]*l\[Nu],RGE\[Kappa]1->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[1]]*l\[Kappa]*(l\[Nu])^2,RGE\[Kappa]2->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[2]]*l\[Kappa]*(l\[Nu])^2,RGEM\[Nu]r->lM\[Nu]r,RGE\[Lambda]1->(lg2^2+lg1^2)/2.,RGE\[Lambda]2->(lg2^2+lg1^2)/2.,RGE\[Lambda]3->(lg2^2-lg1^2)/4,RGE\[Lambda]4->-(lg2^2)/2,RGE\[Lambda]5->0.}];
];
(* changed transition function to 2HDM such that D-terms are correctly taken into account, although there is less flexibility *)
(*,RGE\[Lambda]1->(RGE\[Lambda]1/.Options[lToOpts,RGE\[Lambda]1]),RGE\[Lambda]2->(RGE\[Lambda]2/.Options[lToOpts,RGE\[Lambda]2]),RGE\[Lambda]3->(RGE\[Lambda]3/.Options[lToOpts,RGE\[Lambda]3]),RGE\[Lambda]4->(RGE\[Lambda]4/.Options[lToOpts,RGE\[Lambda]4]),RGE\[Lambda]5->(RGE\[Lambda]5/.Options[lToOpts,RGE\[Lambda]5]) *)


ClearAll[Trans2HDM0N];
Trans2HDM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lsb,lcb,lTosb,lTocb,lu,ld,le,l\[Nu],lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the MSSM w/o heavy neutrinos *)

        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
        l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lTocb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];
	lTosb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lToOpts,RGEtan\[Beta]];

	lu=lsb/((RGEzu/.Options[lToOpts,RGEzu]).{lTocb,lTosb});
	ld=lcb/((RGEzd/.Options[lToOpts,RGEzd]).{lTocb,lTosb});
	l\[Nu]=lsb/((RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]]).{lTocb,lTosb});
	le=lcb/lTocb;


        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lu,RGEYd->lYd*ld,RGEYe->lYe*le,RGE\[Kappa]1->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[1]]*l\[Kappa]*(l\[Nu])^2,RGE\[Kappa]2->(RGEz\[Nu]/.Options[lToOpts,RGEz\[Nu]])[[2]]*l\[Kappa]*(l\[Nu])^2,RGE\[Lambda]1->(lg2^2+lg1^2)/2.,RGE\[Lambda]2->(lg2^2+lg1^2)/2.,RGE\[Lambda]3->(lg2^2-lg1^2)/4,RGE\[Lambda]4->-(lg2^2)/2,RGE\[Lambda]5->0.}];
];
(* changed transition function to 2HDM such that D-terms are correctly taken into account, although there is less flexibility *)
(*,RGE\[Lambda]1->(RGE\[Lambda]1/.Options[lToOpts,RGE\[Lambda]1]),RGE\[Lambda]2->(RGE\[Lambda]2/.Options[lToOpts,RGE\[Lambda]2]),RGE\[Lambda]3->(RGE\[Lambda]3/.Options[lToOpts,RGE\[Lambda]3]),RGE\[Lambda]4->(RGE\[Lambda]4/.Options[lToOpts,RGE\[Lambda]4]),RGE\[Lambda]5->(RGE\[Lambda]5/.Options[lToOpts,RGE\[Lambda]5]) *)

ClearAll[TransSMTriplet];
TransSMTriplet[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut, lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)


(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        Print["The transition from MSSMTriplet to SMTriplet is not fully implemented. The functions to define the parameters in the Higgs potential are not fully implemented."];

        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
	lToIntegratedOut=RGEIntegratedOut/.Options[lToOpts,RGEIntegratedOut];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
	lFromIntegratedOut=RGEIntegratedOut/.Options[lFromOpts,RGEIntegratedOut];
	
(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

	If[lFromIntegratedOut>lToIntegratedOut,
		   (* Print["The model ", RGEGetModel[Exp[pScale]], " at the
		   scale ", Exp[pScale]," pretends to have more particles, but
		   new particles can not be added:
		   ",lFromIntegratedOut,"->",lToIntegratedOut, " Thus the model
		   is changed."]; *)
		   RGEChangeOptions[Exp[pScale],RGEIntegratedOut->lFromIntegratedOut];
		   lToIntegratedOut=lFromIntegratedOut;
		   (*Throw[{lFromIntegratedOut,lToIntegratedOut},RGECanNotAddNewParticles];*)
	   ];
	If[lToIntegratedOut>lFromIntegratedOut,
		lIntegrateOut=lToIntegratedOut-lFromIntegratedOut;

			l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,lIntegrateOut];
			lY\[Nu]= RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
			lM\[Nu]r= RGEIntegrateOutM[lM\[Nu]rRotated, lIntegrateOut];
	];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[Nu]->lY\[Nu]*lsb,RGEY\[CapitalDelta]->lY\[CapitalDelta]*(lsb)^2,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->lM\[Nu]r,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]]),RGE\[CapitalLambda]6->l\[CapitalLambda]u*lM\[CapitalDelta],RGEM\[CapitalDelta]2->lM\[CapitalDelta]^2}];
];


ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut, lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)


(* evaluate the options *)
(* evaluate IntegratedOut in pToOpts and pFromOpts *)
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
	lToIntegratedOut=RGEIntegratedOut/.Options[lToOpts,RGEIntegratedOut];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
	lFromIntegratedOut=RGEIntegratedOut/.Options[lFromOpts,RGEIntegratedOut];
	
(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];

        If[lFromIntegratedOut>lToIntegratedOut,
		   (* Print["The model ", RGEGetModel[Exp[pScale]], " at the
		   scale ", Exp[pScale]," pretends to have more particles, but
		   new particles can not be added:
		   ",lFromIntegratedOut,"->",lToIntegratedOut, " Thus the model
		   is changed."]; *)
		   RGEChangeOptions[Exp[pScale],RGEIntegratedOut->lFromIntegratedOut];
		   lToIntegratedOut=lFromIntegratedOut;
		   (*Throw[{lFromIntegratedOut,lToIntegratedOut},RGECanNotAddNewParticles];*)
	   ];
	If[lToIntegratedOut>lFromIntegratedOut,
		lIntegrateOut=lToIntegratedOut-lFromIntegratedOut;

			l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,lIntegrateOut];
			lY\[Nu]= RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
			lM\[Nu]r= RGEIntegrateOutM[lM\[Nu]rRotated, lIntegrateOut];
	];
	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[Nu]->lY\[Nu]*lsb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->lM\[Nu]r,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}];
];





ClearAll[TransSMTriplet0N];
TransSMTriplet0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM w/o heavy neutrinos *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)

        Print["The transition from MSSMTriplet to SMTriplet0N is not fully implemented. The functions to define the parameters in the Higgs potential are not fully implemented."];

        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
(* calculate the new parameters *)

	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

	l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[CapitalDelta]->lY\[CapitalDelta]*(lsb)^2,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]]),RGE\[CapitalLambda]6->l\[CapitalLambda]u*lM\[CapitalDelta],RGEM\[CapitalDelta]2->lM\[CapitalDelta]^2}];
];



ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lsb,lcb,lY\[CapitalDelta],l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]},
(* make a transition from the MSSM to the SM w/o heavy neutrinos *)
(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!
the dependencies are marked by !!! *)


        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["SM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
(* calculate the new parameters *)

	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[CapitalDelta],l\[Kappa],lM\[Nu]r,l\[CapitalLambda]u,l\[CapitalLambda]d,lM\[CapitalDelta]}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]-=2 l\[CapitalLambda]u/lM\[CapitalDelta] * lY\[CapitalDelta];
	l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}];
];




(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[Nu][pScale],Y\[CapitalDelta][pScale],\[Kappa][pScale],M\[Nu]r[pScale],\[CapitalLambda]u[pScale],\[CapitalLambda]d[pScale],M\[CapitalDelta][pScale]};

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
			Y\[CapitalDelta][pBoundary]==RGEY\[CapitalDelta],
			\[Kappa][pBoundary]==RGE\[Kappa],
			M\[Nu]r[pBoundary]==RGEM\[Nu]r,
			\[CapitalLambda]u[pBoundary]==RGE\[CapitalLambda]u,
			\[CapitalLambda]d[pBoundary]==RGE\[CapitalLambda]d,
			M\[CapitalDelta][pBoundary]==RGEM\[CapitalDelta]
			}//.pInitial
			];
];

End[]; (* end of `Private`*)


EndPackage[];
