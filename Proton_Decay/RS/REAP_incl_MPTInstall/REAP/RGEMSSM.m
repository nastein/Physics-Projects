(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)





BeginPackage["REAP`RGEMSSM`",{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`",
"REAP`RGESM`","REAP`RGEMSSM0N`","REAP`RGESM0N`","REAP`RGE2HDM`","REAP`RGE2HDM0N`"(* transtions to these models are possible *)
}];

(* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! package depends on RGESM.m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
(* dependencies are marked by !!! in the file *)


(* register MSSM *)
RGERegisterModel["MSSM","REAP`RGEMSSM`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGE\[Alpha]->`Private`Get\[Alpha],RGEYd->`Private`GetYd,RGECoupling->`Private`GetCoupling,RGEPoleMTop->`Private`GetPoleMTop,RGE\[Epsilon]1->`Private`Get\[Epsilon]1,RGEY\[Nu]->`Private`GetY\[Nu],RGE\[Epsilon]Max->`Private`Get\[Epsilon]Max,RGERaw->`Private`GetRawSolution,RGERawM\[Nu]r->`Private`GetRawM\[Nu]r,RGE\[Epsilon]->`Private`Get\[Epsilon],RGE\[Kappa]->`Private`Get\[Kappa],RGEMe->`Private`GetMe,RGETwistingParameters->`Private`GetTwistingParameters,RGEMu->`Private`GetMu,RGEMd->`Private`GetMd,RGEM1Tilde->`Private`GetM1Tilde,RGERawY\[Nu]->`Private`GetRawY\[Nu],RGEYu->`Private`GetYu,RGEM\[Nu]->`Private`GetM\[Nu],RGEMixingParameters->`Private`GetMixingParameters,RGEM\[Nu]r->`Private`GetM\[Nu]r,RGEYe->`Private`GetYe,RGEAll->`Private`GetSolution,RGE\[Epsilon]1Max->`Private`Get\[Epsilon]1Max},
{{"2HDM0N",`Private`Trans2HDM0N},{"2HDM",`Private`Trans2HDM},{"MSSM0N",`Private`TransMSSM0N},{"SM0N",`Private`TransSM0N},{"MSSM",`Private`TransMSSM},{"SM",`Private`TransSM}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`", "REAP`RGESolver`","REAP`RGEUtilities`", "REAP`RGEParameters`","REAP`RGEInitial`","MixingParameterTools`MPT3x3`","REAP`RGESM`","REAP`RGEMSSM0N`","REAP`RGESM0N`","REAP`RGE2HDM`","REAP`RGE2HDM0N`"}];

ModelName="MSSM";
ModelVariants={"1Loop","2Loop"};
RGE={RGE1Loop,RGE2Loop};

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

(* Get\[Lambda] is not a function of MSSM *)
(* GetM\[CapitalDelta]2 is not a function of MSSM *)
(* GetM\[CapitalDelta] is not a function of MSSM *)
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

(* GetRawY\[CapitalDelta] is not a function of MSSM *)
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

(* Get\[Kappa]1 is not a function of MSSM *)
(* Get\[Kappa]2 is not a function of MSSM *)
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
GetM\[Nu][pScale_,pSolution_,pOpts___]:=Block[{l\[Kappa],lY\[Nu],lM,lvu},
(* returns the mass matrix of the neutrinos *)
	{l\[Kappa],lY\[Nu],lM}=({\[Kappa][pScale],Y\[Nu][pScale],M\[Nu]r[pScale]}/.pSolution)[[1]];
	lOpts;
	Options[lOpts]=Options[RGEOptions];
	SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
        lvu=RGEv\[Nu]/.Options[lOpts];
	If[MatrixConditionNumber[lM]>2*Precision[lM],Print["GetM\[Nu]: The
	matrix M=", MatrixForm[lM]], " is ill-conditioned and the condition
	number is ", MatrixConditionNumber[lM]];

	Return[-lvu^2*(1/2)^2*(10^9*l\[Kappa]+2*Transpose[lY\[Nu]].(Inverse[lM]*10^9).lY\[Nu])];
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
        {lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM}=(ParametersFunc[pScale]/.pSolution)[[1]];
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

ClearAll[Get\[Epsilon]1Max];
Get\[Epsilon]1Max[pScale_,pSolution_,pOpts___]:=Block[{lv,lM,lm,lm1Tilde,lm1,lm2,lm3,leps},
(* returns the mass matrix of the charged leptons *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lv=(RGEv\[Nu]/Sqrt[2])/.Options[lOpts];
   lM=(M\[Nu]r[pScale]/.pSolution)[[1]];
   lM\[Nu]=lv*(Y\[Nu][pScale]/.pSolution)[[1]];
   {lM,lM\[Nu]}=RGERotateM[lM,lM\[Nu]];
   lm1Tilde=(lM\[Nu].Dagger[lM\[Nu]])[[1,1]]/lM[[1,1]];
   {lm1,lm2,lm3}=MNSParameters[GetM\[Nu][pScale,pSolution,pOpts],GetYe[pScale,pSolution,pOpts]][[2]]*10^-9;
   leps=3/8/\[Pi]^2* lM[[1,1]] *lm3/lv^2*(1-lm1/lm3 Sqrt[1+(lm3^2-lm1^2)/lm1Tilde^2]);
   Return[RGEFloor[Abs[leps]][[1]]];
];

ClearAll[Get\[Epsilon]1];
Get\[Epsilon]1[pScale_,pSolution_,pOpts___]:=Block[{lv,lM,lm,lY,leps},
(* returns the mass matrix of the charged leptons *)
   lOpts;
   Options[lOpts]=Options[RGEOptions];
   SetOptions[lOpts,RGEFilterOptions[lOpts,pOpts]];
   lv=(RGEv\[Nu]/Sqrt[2])/.Options[lOpts];
   lM=(M\[Nu]r[pScale]/.pSolution)[[1]];
   lY=(Y\[Nu][pScale]/.pSolution)[[1]];
   {lM,lY}=RGERotateM[lM,lY];
   lm=GetM\[Nu][pScale,pSolution,pOpts]*10^-9;
   lY=lY[[1]];
   leps=3/8/\[Pi]*lM[[1,1]]/lv^2*Im[(lY.Conjugate[lm].lY)]/(lY.Conjugate[lY]);
   Return[leps];
];

(* GetGWCond is not a function of MSSM *)
(* GetGWConditions is not a function of MSSM *)
(* GetVEVratio is not a function of MSSM *)
(* GetVEVratios is not a function of MSSM *)
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

                        
Parameters={g1,g2,g3,Yu,Yd,Ye,Y\[Nu],\[Kappa],M\[Nu]r};
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYu,RGEYd,RGEYe,RGEY\[Nu],RGE\[Kappa],RGEM\[Nu]r};

ClearAll[Initial];
Initial={
{"GUT",{
	RGEg1->0.7044110331165641,
	RGEg2->0.6965468498179075,
	RGEg3->0.6983661130877465,
	RGEYd->RGEGetYd[RGEyd,RGEys,RGEyb,RGEq\[Theta]12,RGEq\[Theta]13,RGEq\[Theta]23,RGEq\[Delta],RGEq\[Delta]e,RGEq\[Delta]\[Mu],RGEq\[Delta]\[Tau],RGEq\[CurlyPhi]1,RGEq\[CurlyPhi]2],
	RGEYu->DiagonalMatrix[{RGEyu,RGEyc,RGEyt}],
	RGE\[Kappa]->0*IdentityMatrix[3],
	RGEM\[Nu]r->RGEGetM[RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu,RGEY\[Nu]],
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
},
{"MZ",{
	RGE\[Kappa] -> 0*IdentityMatrix[3], 
	RGEYe->DiagonalMatrix[{RGEye,RGEy\[Mu],RGEy\[Tau]}],
	RGEY\[Nu]33 -> 0.1,
	RGEY\[Nu]Ratio -> 0.5,
	RGEY\[Nu]->RGEGetY\[Nu][RGEY\[Nu]33,RGEY\[Nu]Ratio],
	RGEM\[Nu]r->RGEGetM[RGE\[Theta]12,RGE\[Theta]13,RGE\[Theta]23,RGE\[Delta],RGE\[Delta]e,RGE\[Delta]\[Mu],RGE\[Delta]\[Tau],RGE\[CurlyPhi]1,RGE\[CurlyPhi]2,RGEMlightest,RGE\[CapitalDelta]m2atm,RGE\[CapitalDelta]m2sol,RGEMassHierarchy,RGEvu,RGEY\[Nu]],
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

}; (* a list containing suggestions for initial values *)



ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Y\[Nu][t],t]==BetaY\[Nu][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[M\[Nu]r[t],t]==BetaM\[Nu]r[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]
}; (* renormalization group equations of the MSSM ( 1 Loop ) *)



(* Beta functions of the MSSM *)
ClearAll[Betag1, Betag2, Betag3, BetaYu, BetaYd, BetaYe, BetaY\[Nu], Beta\[Kappa], BetaM\[Nu]r, Beta\[Lambda]];

(* 1 loop contributions *)

Betag1[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=
	(33/5) * 1/(16*Pi^2) * g1^3;

Betag2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=
	(1) * 1/(16*Pi^2) * g2^3;

Betag3[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=
	(-3) * 1/(16*Pi^2) * g3^3;


BetaYd[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] := 1/(16*Pi^2) * (
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

BetaYu[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] := 1/(16*Pi^2) * (
          Yu.(
          + Dagger[Yd].Yd
	  + 3*Dagger[Yu].Yu
          )
          + (
          - (13/15)*g1^2
	  - 3*g2^2
	  - (16/3)*g3^2
	  + Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*Tr[Dagger[Yu].Yu]
          )*Yu
          );

BetaYe[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] := 1/(16*Pi^2) * (
          Ye.(
          + 3*Dagger[Ye].Ye
	  + Dagger[Y\[Nu]].Y\[Nu]
          )
          + (
          - (9/5)*g1^2
	  - 3*g2^2
	  + 3*Tr[Dagger[Yd].Yd]
	  + Tr[Dagger[Ye].Ye]
          )*Ye
          );

BetaY\[Nu][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] := 1/(16*Pi^2) * (
          Y\[Nu].(
          + Dagger[Ye].Ye
	  + 3*Dagger[Y\[Nu]].Y\[Nu]
          )
          + (
          - (3/5)*g1^2
	  - 3*g2^2
	  + Tr[Dagger[Y\[Nu]].Y\[Nu]]
	  + 3*Tr[Dagger[Yu].Yu]
          )*Y\[Nu]
          );

Beta\[Kappa][g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(16*Pi^2) * (
	\[Kappa].(
        + Dagger[Ye].Ye
        + Dagger[Y\[Nu]].Y\[Nu]
        )
	+ (
        + Transpose[Ye].Conjugate[Ye]
	+ Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
        ).\[Kappa]
        + (
	- (6/5)*g1^2
	- 6*g2^2
	+ 6*Tr[Dagger[Yu].Yu]
        + 2*Tr[Dagger[Y\[Nu]].Y\[Nu]]
        )*\[Kappa]
        );

BetaM\[Nu]r[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] := 1/(16*Pi^2) * (
          2*M\[Nu]r.Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	  + 2*Y\[Nu].Dagger[Y\[Nu]].M\[Nu]r);




ClearAll[RGE2Loop];
RGE2Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[g3[t],t]==Betag3[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Yu[t],t]==BetaYu[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+BetaYu2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Yd[t],t]==BetaYd[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+BetaYd2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Ye[t],t]==BetaYe[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+BetaYe2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[Y\[Nu][t],t]==BetaY\[Nu][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+BetaY\[Nu]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[\[Kappa][t],t]==Beta\[Kappa][g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+Beta\[Kappa]2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]],
		D[M\[Nu]r[t],t]==BetaM\[Nu]r[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]+BetaM\[Nu]r2[g1[t],g2[t],g3[t],Yu[t],Yd[t],Ye[t],Y\[Nu][t],\[Kappa][t],M\[Nu]r[t]]
}; (* renormalization group equations of the MSSM ( 2 Loop ) *)


ClearAll[BetaYu2, BetaYd2, BetaYe2, BetaY\[Nu]2, Beta\[Kappa]2, BetaM\[Nu]r2];


(* 2 loop contributions *)

BetaYu2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	 + Yu.(
         - 2*Dagger[Yd].Yd.Dagger[Yd].Yd
	 - 2*Dagger[Yd].Yd.Dagger[Yu].Yu
	 - 4*Dagger[Yu].Yu.Dagger[Yu].Yu
	 - 3*Dagger[Yd].Yd*Tr[Yd.Dagger[Yd]]
	 - Dagger[Yd].Yd*Tr[Ye.Dagger[Ye]]
	 - 9*Dagger[Yu].Yu*Tr[Yu.Dagger[Yu]]
	 - 3*Dagger[Yu].Yu*Tr[Y\[Nu].Dagger[Y\[Nu]]]
         + (2/5)*g1^2*Dagger[Yd].Yd
	 + (2/5)*g1^2*Dagger[Yu].Yu
	 + 6*g2^2*Dagger[Yu].Yu
         )
         + Yu*(
	 - 3*Tr[Dagger[Yu].Yd.Dagger[Yd].Yu]
	 - 9*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
	 - Tr[Y\[Nu].Dagger[Ye].Ye.Dagger[Y\[Nu]]]
	 - 3*Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
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


BetaYd2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	+Yd.(
        - 4*Dagger[Yd].Yd.Dagger[Yd].Yd
	- 2*Dagger[Yu].Yu.Dagger[Yd].Yd
	- 2*Dagger[Yu].Yu.Dagger[Yu].Yu
	- 9*Dagger[Yd].Yd*Tr[Yd.Dagger[Yd]]
	- 3*Dagger[Yd].Yd*Tr[Ye.Dagger[Ye]]
	- Dagger[Yu].Yu*Tr[Y\[Nu].Dagger[Y\[Nu]]]
	- 3*Dagger[Yu].Yu*Tr[Yu.Dagger[Yu]]
	+ 6*Dagger[Yd].Yd*g2^2
	+ (4/5)*g1^2*Dagger[Yd].Yd
	+ (4/5)*g1^2*Dagger[Yu].Yu
        )
	+Yd*(
	- 9*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
	- 3*Tr[Dagger[Yd].Yu.Dagger[Yu].Yd]
	- 3*Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
	- Tr[Ye.Dagger[Y\[Nu]].Y\[Nu].Dagger[Ye]]
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


BetaYe2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	+ Ye.(
        - 4*Dagger[Ye].Ye.Dagger[Ye].Ye
	- 2*Dagger[Y\[Nu]].Y\[Nu].Dagger[Ye].Ye
	- 2*Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]
	- 9*Dagger[Ye].Ye*Tr[Yd.Dagger[Yd]]
	- 3*Dagger[Ye].Ye*Tr[Ye.Dagger[Ye]]
	- Dagger[Y\[Nu]].Y\[Nu]*Tr[Y\[Nu].Dagger[Y\[Nu]]]
	- 3*Dagger[Y\[Nu]].Y\[Nu]*Tr[Yu.Dagger[Yu]]
	+ 6*g2^2*Dagger[Ye].Ye
        )
        +Ye*(
	- 9*Tr[Dagger[Yd].Yd.Dagger[Yd].Yd]
	- 3*Tr[Dagger[Yd].Yu.Dagger[Yu].Yd]
	- 3*Tr[Dagger[Ye].Ye.Dagger[Ye].Ye]
	- Tr[Ye.Dagger[Y\[Nu]].Y\[Nu].Dagger[Ye]]
	+ (6/5)*g1^2*Tr[Dagger[Ye].Ye]
        - (2/5)*g1^2*Tr[Dagger[Yd].Yd]
	+ 16*g3^2*Tr[Dagger[Yd].Yd]
	+ (27/2)*g1^4
        + (9/5)*g1^2*g2^2
	+ (15/2)*g2^4
        )
        );


BetaY\[Nu]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	+ Y\[Nu].(
        - 2*Dagger[Ye].Ye.Dagger[Ye].Ye
	- 2*Dagger[Ye].Ye.Dagger[Y\[Nu]].Y\[Nu]
	- 4*Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]
	- 3*Dagger[Ye].Ye*Tr[Yd.Dagger[Yd]]
	- Dagger[Ye].Ye*Tr[Ye.Dagger[Ye]]
	- 3*Dagger[Y\[Nu]].Y\[Nu]*Tr[Y\[Nu].Dagger[Y\[Nu]]]
	- 9*Dagger[Y\[Nu]].Y\[Nu]*Tr[Yu.Dagger[Yu]]
	+ (6/5)*g1^2*Dagger[Ye].Ye
	+ (6/5)*g1^2*Dagger[Y\[Nu]].Y\[Nu]
	+ 6*g2^2*Dagger[Y\[Nu]].Y\[Nu]
        )
        + Y\[Nu]*(
	- Tr[Y\[Nu].Dagger[Ye].Ye.Dagger[Y\[Nu]]]
	- 3*Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
	- 3*Tr[Dagger[Yu].Yd.Dagger[Yd].Yu]
	- 9*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
        + (4/5)*g1^2*Tr[Dagger[Yu].Yu]
	+ 16*g3^2*Tr[Dagger[Yu].Yu]
	+ (207/50)*g1^4
        + (9/5)*g1^2*g2^2
	+ (15/2)*g2^4
        )
        );


Beta\[Kappa]2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	+ \[Kappa].(
        - 2*Dagger[Ye].Ye.Dagger[Ye].Ye
        - 2*Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]
	- (Tr[Y\[Nu].Dagger[Y\[Nu]]]+3*Tr[Yu.Dagger[Yu]])*Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
	+ ((6/5)*g1^2 - Tr[Ye.Dagger[Ye]]- 3*Tr[Yd.Dagger[Yd]])*Transpose[Ye].Conjugate[Ye]
        )
        + (
	- 2*Transpose[Ye].Conjugate[Ye].Transpose[Ye].Conjugate[Ye]
	- 2*Transpose[Y\[Nu]].Conjugate[Y\[Nu]].Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
	- (Tr[Y\[Nu].Dagger[Y\[Nu]]]+3*Tr[Yu.Dagger[Yu]])*Transpose[Y\[Nu]].Conjugate[Y\[Nu]]
	+ ((6/5)*g1^2 - Tr[Ye.Dagger[Ye]]- 3*Tr[Yd.Dagger[Yd]])*Transpose[Ye].Conjugate[Ye]
        ).\[Kappa]
	+ \[Kappa]*(
	- 6*Tr[Dagger[Yu].Yd.Dagger[Yd].Yu]
	- 18*Tr[Dagger[Yu].Yu.Dagger[Yu].Yu]
	- 2*Tr[Y\[Nu].Dagger[Ye].Ye.Dagger[Y\[Nu]]]
	- 6*Tr[Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]].Y\[Nu]]
        + (8/5)*g1^2*Tr[Dagger[Yu].Yu]
	+ 32*g3^2*Tr[Dagger[Yu].Yu]
        + (207/25)*g1^4
	+ (18/5)*g1^2*g2^2
	+ 15*g2^4
        )
        );


BetaM\[Nu]r2[g1_,g2_,g3_,Yu_,Yd_,Ye_,Y\[Nu]_,\[Kappa]_,M\[Nu]r_] :=1/(4*Pi)^4 * (
	 + M\[Nu]r.(
         - 2*Conjugate[Y\[Nu]].Transpose[Ye].Conjugate[Ye].Transpose[Y\[Nu]]
	 - 2*Conjugate[Y\[Nu]].Transpose[Y\[Nu]].Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	 - 6*Tr[Yu.Dagger[Yu]]*Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	 - 2*Tr[Y\[Nu].Dagger[Y\[Nu]]]*Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	 + (6/5)*g1^2*Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	 + 6*g2^2*Conjugate[Y\[Nu]].Transpose[Y\[Nu]]
	 )
         + (
         - 2*Y\[Nu].Dagger[Ye].Ye.Dagger[Y\[Nu]]
	 - 2*Y\[Nu].Dagger[Y\[Nu]].Y\[Nu].Dagger[Y\[Nu]]
	 - 6*Tr[Yu.Dagger[Yu]]*Y\[Nu].Dagger[Y\[Nu]]
	 - 2*Tr[Y\[Nu].Dagger[Y\[Nu]]]*Y\[Nu].Dagger[Y\[Nu]]
	 + (6/5)*g1^2*Y\[Nu].Dagger[Y\[Nu]]
	 + 6*g2^2*Y\[Nu].Dagger[Y\[Nu]]
         ).M\[Nu]r
	 );

	  

(* transition functions *)

ClearAll[TransMSSM];
TransMSSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
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

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGEY\[Nu]->lY\[Nu],RGE\[Kappa]->l\[Kappa],RGEM\[Nu]r->lM\[Nu]r}];
];



ClearAll[TransMSSM0N];
TransMSSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM},
(* make a transition from the MSSM to the MSSM w/o heavy neutrinos *)

(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

        l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa]}];
];


ClearAll[Trans2HDM];
Trans2HDM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut,le,lu,ld,l\[Nu],lTosb,lTocb,lsb,lcb},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
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


ClearAll[Trans2HDM0N];
Trans2HDM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lsb,lcb,lTosb,lTocb,lu,ld,le,l\[Nu]},
(* make a transition from the MSSM to the MSSM w/o heavy neutrinos *)

        lToOpts;
        Options[lToOpts]=Options[RGEGetModelOptions["2HDM0N"][[1,2]]];
        SetOptions[lToOpts,RGEFilterOptions[lToOpts,pToOpts]];
        lFromOpts;
        Options[lFromOpts]=Options[RGEOptions];
        SetOptions[lFromOpts,RGEFilterOptions[lFromOpts,pFromOpts]];
(* calculate the new parameters *)
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

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


ClearAll[TransSM];
TransSM[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated,lIntegrateOut, lUforM,lToIntegratedOut,lFromIntegratedOut, lsb,lcb},
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
	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
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
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGEY\[Nu]->lY\[Nu]*lsb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGEM\[Nu]r->lM\[Nu]r,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}];
];



ClearAll[TransSM0N];
TransSM0N[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],lY\[Nu]Rotated,l\[Kappa],lM\[Nu]r,lM\[Nu]rRotated, lUforM,lsb,lcb},
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

	{lg1,lg2,lg3,lYu,lYd,lYe,lY\[Nu],l\[Kappa],lM\[Nu]r}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
	{lM\[Nu]rRotated,lY\[Nu]Rotated}=RGERotateM[ lM\[Nu]r,lY\[Nu] ]; (*rotation matrix for lM\[Nu]r*)

	l\[Kappa]+= RGEKappaMatching[lM\[Nu]rRotated,lY\[Nu]Rotated,Length[lM\[Nu]r]];

	lcb=1/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	lsb=RGEtan\[Beta]/Sqrt[1+RGEtan\[Beta]^2]/.Options[lFromOpts,RGEtan\[Beta]];
	Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu*lsb,RGEYd->lYd*lcb,RGEYe->lYe*lcb,RGE\[Kappa]->l\[Kappa]*(lsb)^2,RGE\[Lambda]->(RGE\[Lambda]/.Options[lToOpts,RGE\[Lambda]])}];
];




(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yu[pScale],Yd[pScale],Ye[pScale],Y\[Nu][pScale],\[Kappa][pScale],M\[Nu]r[pScale]};

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
			\[Kappa][pBoundary]==RGE\[Kappa],
			M\[Nu]r[pBoundary]==RGEM\[Nu]r
			}//.pInitial
			];
];

End[]; (* end of `Private`*)


EndPackage[];
