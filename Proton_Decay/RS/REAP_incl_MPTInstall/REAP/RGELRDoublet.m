(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)



BeginPackage["REAP`RGELRDoublet`",{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEParameters`","REAP`RGEUtilities`","MixingParameterTools`MPT3x3`","REAP`RGEInitial`","REAP`RGE2HDMDirac`"}];


(* register LRDoublet *)
RGERegisterModel["LRDoublet","REAP`RGELRDoublet`",
	`Private`GetParameters,
        `Private`SolveModel,
        {RGEGWConditions->`Private`GetGWConditions,RGECoupling->`Private`GetCoupling,RGEVEVratio->`Private`GetVEVratio,RGEAll->`Private`GetRawSolution,RGE\[Alpha]->`Private`Get\[Alpha],RGEGWCondition->`Private`GetGWCond,RGEVEVratios->`Private`GetVEVratios},
{{"2HDMDirac",`Private`Trans2HDMDirac}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
         ];


Begin["`Private`"];
Map[Needs,{"REAP`RGESymbol`","REAP`RGESolver`","REAP`RGEParameters`","REAP`RGEUtilities`","MixingParameterTools`MPT3x3`","REAP`RGEInitial`","REAP`RGE2HDMDirac`"}];

ModelName="LRDoublet";
ModelVariants={"1Loop"};
RGE={RGE1Loop};

ClearAll[GetRawSolution];
GetRawSolution[pScale_,pSolution_,pOpts___]:=Block[{},
(* returns all parameters of the SM *)
        Return[(ParametersFunc[pScale]/.pSolution)[[1]]];
];

(* GetRawM\[Nu]r is not a function of LRDoublet *)
(* GetRawY\[Nu] is not a function of LRDoublet *)
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

(* Get\[Lambda] is not a function of LRDoublet *)
(* GetM\[CapitalDelta]2 is not a function of LRDoublet *)
(* GetM\[CapitalDelta] is not a function of LRDoublet *)
(* GetYe is not a function of LRDoublet *)
(* GetYu is not a function of LRDoublet *)
(* GetYd is not a function of LRDoublet *)
(* GetRawY\[CapitalDelta] is not a function of LRDoublet *)
(* GetY\[Nu] is not a function of LRDoublet *)
(* Get\[Kappa] is not a function of LRDoublet *)
(* Get\[Kappa]1 is not a function of LRDoublet *)
(* Get\[Kappa]2 is not a function of LRDoublet *)
(* GetM\[Nu]r is not a function of LRDoublet *)
(* GetM\[Nu] is not a function of LRDoublet *)
(* GetMu is not a function of LRDoublet *)
(* GetPoleMTop is not a function of LRDoublet *)
(* GetMd is not a function of LRDoublet *)
(* GetMe is not a function of LRDoublet *)
(* GetSolution is not a function of LRDoublet *)
(* GetMixingParameters is not a function of LRDoublet *)
(* GetTwistingParameters is not a function of LRDoublet *)
(* GetM1Tilde is not a function of LRDoublet *)
(* Get\[Epsilon]1Max is not a function of LRDoublet *)
(* Get\[Epsilon]1 is not a function of LRDoublet *)
ClearAll[GetGWCond];
GetGWCond[pScale_,pSolution_,pOpts___]:=Block[{lk1,lk2,lb1,ll1,ll2,ll3,lf1},
(* returns the fine structure constants *)
    {lk1,lk2,lb1,ll1,ll2,ll3,lf1}=Flatten[{\[Kappa]1[pScale],\[Kappa]2[pScale],\[Beta]1[pScale],\[Lambda]1[pScale],\[Lambda]2,\[Lambda]3,f1[pScale]}/.pSolution];
    Return[(lk1+lk2-(lf1-2 lb1)^2/8/ll1)];
];

ClearAll[GetGWConditions];
GetGWConditions[pScale_,pSolution_,pOpts___]:=Block[{lk1,lk2,lb1,ll1,ll2,ll3,lf1},
(* returns the fine structure constants *)
    {lk1,lk2,lb1,ll1,ll2,ll3,lf1}=Flatten[{\[Kappa]1[pScale],\[Kappa]2[pScale],\[Beta]1[pScale],\[Lambda]1[pScale],\[Lambda]2[pScale],\[Lambda]3[pScale],f1[pScale]}/.pSolution];
    Return[{
{
lk1+lk2-lb1^2/2/(ll1+4ll2)+lf1^2/32/ll2,
lk1-lb1^2/2/(ll1+4ll2)+lf1^2/32/ll2,
lk1+lk2-lb1^2/2/(ll1-4ll3)-lf1^2/32/ll3,
lk1-lb1^2/2/(ll1-4ll3)-lf1^2/32/ll3,
ll1+4ll2,
ll1-4ll3
},
{
lk1+lk2-(lf1-2 lb1)^2/8/ll1,
lk1-(lf1-2lb1)^2/8/ll1,
lk1+lk2-(2lb1+lf1)^2/8/ll1,
lk1-(2lb1+lf1)^2/8/ll1,
lk1+lk2,
lk1,
ll1,
ll1}}];
];

ClearAll[GetVEVratio];
GetVEVratio[pScale_,pSolution_,pOpts___]:=Block[{lb1,ll1,lf1},
(* returns the fine structure constants *)
    {lb1,ll1,lf1}=Flatten[{\[Beta]1[pScale],\[Lambda]1[pScale],f1[pScale]}/.pSolution];
    Return[((lf1-2 lb1)/4/ll1)];
];

ClearAll[GetVEVratios];
GetVEVratios[pScale_,pSolution_,pOpts___]:=Block[{lb1,ll1,lf1},
(* returns the fine structure constants *)
    {lb1,ll1,lf1}=Flatten[{\[Beta]1[pScale],\[Lambda]1[pScale],f1[pScale]}/.pSolution];
    Return[((lf1-2 lb1)/4/ll1)];
];

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


(* definitions for the LR Doublet Model (LRDoublet) *)

ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={		  RGEModelVariant->"1Loop", (* different variation of the model *)
			  RGEAutoGenerated->False, (* used to find automatically generated entries *)
				  RGEvEW->246, (* vev for the electroweak transition *)
				  Method->StiffnessSwitching  (* option of NDSolve *)
};

Parameters={g1,g2,g3,Yqp,Yqm,Ylp,Ylm, \[Lambda]1,\[Lambda]2, \[Lambda]3, \[Lambda]4,\[Kappa]1,\[Kappa]2,\[Beta]1,\[Beta]2,\[Beta]3,f1}; 
ParameterSymbols={RGEg1,RGEg2,RGEg3,RGEYqp,RGEYqm,RGEYlp,RGEYlm, RGE\[Lambda]1,RGE\[Lambda]2, RGE\[Lambda]3, RGE\[Lambda]4,RGE\[Kappa]1,RGE\[Kappa]2,RGE\[Beta]1,RGE\[Beta]2,RGE\[Beta]3,RGEf1}; 


ClearAll[Initial];
Initial={
{"MPl",{ RGEg1 -> 0.7044
	, RGEg2 -> 0.6965
	, RGEg3 -> 0.6984
	, RGEYqp -> {{5.980 10^-6, 0, 0},{0, 0.001736, 0},{0, 0, 0.7417}}
	, RGEYqm -> {{5.980 10^-6, 0, 0},{0, 0.001736, 0},{0, 0, 0.7417}}
	, RGEYlp -> {{0.00009344, 0, 0},{0, 0.01972, 0},{0, 0, 0.33678}}
	, RGEYlm -> {{0.00009344, 0, 0},{0, 0.01972, 0},{0, 0, 0.33678}}
	, RGE\[Lambda]1 -> 1
	, RGE\[Lambda]2 -> 1
	, RGE\[Lambda]3 -> -0.1
	, RGE\[Lambda]4 -> 0
	, RGE\[Kappa]1 -> 1
	, RGE\[Kappa]2 -> -0.5
	, RGE\[Beta]1 -> 0.1
	, RGE\[Beta]2 -> 0
	, RGE\[Beta]3 -> 0
	, RGEf1 -> 0.01
	}
}
};


ClearAll[RGE1Loop];
RGE1Loop:={	D[g1[t],t]==Betag1[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[g2[t],t]==Betag2[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[g3[t],t]==Betag2[g3[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[Yqp[t],t]==BetaYqp[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[Yqm[t],t]==BetaYqm[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[Ylp[t],t]==BetaYlp[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[Ylm[t],t]==BetaYlm[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Lambda]1[t],t]==Beta\[Lambda]1[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Lambda]2[t],t]==Beta\[Lambda]2[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Lambda]3[t],t]==Beta\[Lambda]3[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Lambda]4[t],t]==Beta\[Lambda]4[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Kappa]1[t],t]==Beta\[Kappa]1[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Kappa]2[t],t]==Beta\[Kappa]2[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Beta]1[t],t]==Beta\[Beta]1[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Beta]2[t],t]==Beta\[Beta]2[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[\[Beta]3[t],t]==Beta\[Beta]3[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
		, D[f1[t],t]==Betaf1[g1[t],g2[t],g3[t],Yqp[t],Yqm[t],Ylp[t],Ylm[t], \[Lambda]1[t],\[Lambda]2[t], \[Lambda]3[t], \[Lambda]4[t],\[Kappa]1[t],\[Kappa]2[t],\[Beta]1[t],\[Beta]2[t],\[Beta]3[t],f1[t]]
};

(* Beta Functions of the LR Doublet model *)
ClearAll[Betag1, Betag2, Betag3, BetaYqp, BetaYqm, BetaYlp, BetaYlm, Beta\[Lambda]1,Beta\[Lambda]2, Beta\[Lambda]3, Beta\[Lambda]4,Beta\[Kappa]1,Beta\[Kappa]2,Beta\[Beta]1,Beta\[Beta]2,Beta\[Beta]3,Betaf1 ];

Betag1[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := g1^3 (3)/16 /\[Pi]^2;
Betag2[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := -(17/6)/16 /\[Pi]^2 g2^3; 
Betag3[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := -(7)/16 /\[Pi]^2 g3^3;

BetaYqp[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/64/\[Pi]^2 * (
(
-(2/9)g1^2
-9 g2^2 
-32 g3^2 
+ 3 ( Tr[Yqp.Yqp]+ Tr[Yqm.Yqm])
+  (Tr[Ylp.Ylp] + Tr[Ylm.Ylm])
) Yqp
+ (6 Tr[Yqp.Yqm]+2 Tr[Ylp.Ylm])* Yqm
+4 Yqp.Yqp.Yqp 
-2 Yqp.Yqm.Yqm
-2 Yqm.Yqm.Yqp
);


BetaYqm[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/64/\[Pi]^2 * (
(
-(2/9)g1^2
-9 g2^2 
-32 g3^2 
+ 3 (Tr[Yqp.Yqp]+Tr[ Yqm.Yqm])
+  (Tr[Ylp.Ylp]+ Tr[Ylm.Ylm])
) Yqm
+ (6 Tr[Yqp.Yqm] + 2 Tr[Ylp.Ylm])* Yqp
+4 Yqm.Yqm.Yqm 
-2 Yqm.Yqp.Yqp
-2 Yqp.Yqp.Yqm
);

BetaYlp[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/64/\[Pi]^2 * (
(
-6g1^2
-9 g2^2 
+ 3 (Tr[Yqp.Yqp]+ Tr[Yqm.Yqm])
+  (Tr[Ylp.Ylp] + Tr[Ylm.Ylm])
) Ylp
+ (6 Tr[ Yqp.Yqm] + 2 Tr[Ylp.Ylm]) Ylm
+4Ylp.Ylp.Ylp 
-2 Ylp.Ylm.Ylm
-2 Ylm.Ylm.Ylp
);

BetaYlm[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/64/\[Pi]^2 * (
(
-6g1^2
-9 g2^2 
+ 3 (Tr[Yqp.Yqp] + Tr[Yqm.Yqm])
+  (Tr[Ylp.Ylp] + Tr[Ylm.Ylm])
) Ylm
+ (6 Tr[Yqp.Yqm] + 2 Tr[Ylp.Ylm]) Ylp
+4Ylm.Ylm.Ylm 
-2 Ylm.Ylp.Ylp
-2 Ylp.Ylp.Ylm
);

Beta\[Lambda]1[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(128 \[Pi]^2) * 
(
+ 32 b1^2 +8 f1^2 + 9 g2^4 - 72 g2^2 l1 
+  8 ( 
   + 32 l1^2 + 128 l2^2 + 128 l3^2 + 48 l4^2 
   - 2 l4 ( 3 Tr[Yqp.Yqm] + Tr[Ylp.Ylm])
   + l1 (32 l2 - 32 l3 + 3 (Tr[Yqp.Yqp]+ Tr[Yqm.Yqm]) 
   + ( Tr[Ylp.Ylp] + Tr[Ylm.Ylm])
   )
) 
- 12 (Tr[Yqp.Yqp.Yqp.Yqp]+Tr[Yqm.Yqm.Yqm.Yqm]) 
-4( Tr[Ylp.Ylp.Ylp.Ylp]+Tr[Ylm.Ylm.Ylm.Ylm]));

Beta\[Lambda]2[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(512 \[Pi]^2) *
(
+128 b2^2 -8 f1^2 + 3 g2^4 - 288 g2^2 l2 + 768 l1 l2 + 3072 l2^2 + 1024 l2 l3 
+ 384 l4^2 
+ 32 l2  * (
  +3 (Tr[Yqp.Yqp] +  Tr[Yqm.Yqm])+ Tr[Ylp.Ylp]+ Tr[Ylm.Ylm]) 
  - 16 l4 (3 Tr[Yqp.Yqm] +Tr[Ylp.Ylm]) 
  + (1/2)* (
    (
      + Tr[Ylp.Ylp.Ylp.Ylp] 
      + Tr[Ylm.Ylm.Ylm.Ylm] 
      - 2 Tr[Ylp.Ylm.Ylp.Ylm] 
      -4 Tr[Ylp.Ylp.Ylm.Ylm]
      )
  +3*(
	+ Tr[Yqp.Yqp.Yqp.Yqp] 
  	+ Tr[Yqm.Yqm.Yqm.Yqm] 
	- 2 Tr[Yqp.Yqm.Yqp.Yqm] 
	-4 Tr[Yqp.Yqp.Yqm.Yqm])
 ) 
);

Beta\[Lambda]3[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/(256 \[Pi]^2) * 
( 
+ 32 b3^2 + 4 f1^2 - 3 g2^4 - 144 g2^2 l3 + 384 l1 l3 -512 l2 l3 - 1536 l3^2 
+  16 l3  ((3 Tr[Yqp.Yqp] + Tr[Yqm.Yqm])+Tr[Ylp.Ylp]+Tr[Ylm.Ylm]) 
- 3 (
  + Tr[Yqp.Yqp.Yqp.Yqp]
  + Tr[Yqm.Yqm.Yqm.Yqm]
  - 4 Tr[Yqp.Yqp.Yqm.Yqm] 
  + 2 Tr[Yqp.Yqm.Yqp.Yqm] 
  )
-  (
  + Tr[Ylp.Ylp.Ylp.Ylp]
  + Tr[Ylm.Ylm.Ylm.Ylm]
  - 4 Tr[Ylp.Ylp.Ylm.Ylm] 
  + 2 Tr[Ylp.Ylm.Ylp.Ylm] 
  ) 
);

Beta\[Lambda]4[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/(64 \[Pi]^2) * 
( 
+ l4 (
  -36 g2^2 +192 l1 +768 l2 
  + 4 (Tr[Ylp.Ylp]+Tr[Ylm.Ylm]+ 3(Tr[Yqp.Yqp]+Tr[Yqm.Yqm]))
  )
  - 4 (3Tr[ Yqp.Yqm] + Tr[Ylp.Ylm]) (2 l1 + 4 l2) 
  + 32 b1 b2 
  +2 (
  + Tr[Ylm.Ylp.Ylp.Ylp]
  + Tr[Ylp.Ylm.Ylm.Ylm] 
  + 3 (Tr[Yqm.Yqp.Yqp.Yqp]+Tr[Yqp.Yqm.Yqm.Yqm])) 
);

Beta\[Kappa]1[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(512 \[Pi]^2) (
 + 256 b1^2 + 1024 b2^2 + 128 f1^2 + 24 g1^4 + 12 g1^2 g2^2 + 9 g2^4 
 - 96 g1^2 k1 - 144 g2^2 k1 + 576 k1^2 + 384 k1 k2 + 192 k2^2
 );

Beta\[Kappa]2[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(512 \[Pi]^2) (
 -1024 b3^2 + 128 f1^2 + 12 g1^2 g2^2 + 9 g2^4 
 -  96 g1^2 k2 - 144 g2^2 k2 + 512 k1 k2 + 384 k2^2
);

Beta\[Beta]1[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(256 \[Pi]^2) (
 + 32 b1^2 + 128 b2^2 - 128 b3^2 + 24 f1^2 + 9 g2^4 
 + 16 b2 (48  l4 -   6 Tr[Yqp.Yqm] + 2 Tr[Ylp.Ylm]) 
 - 4 b1 (
   + 6 g1^2 + 27 g2^2 
   - 2 (
      +	20 k1 + 4 k2 + 40 l1 + 32 l2 - 32 l3 
      +  3 (Tr[ Yqp.Yqp]+ Tr[Yqm.Yqm]) +  (Tr[Ylp.Ylp] + Tr[Ylm.Ylm])
      )
    )
);

Beta\[Beta]2[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
 1/(64 \[Pi]^2) (
 + b2 (
    -6 g1^2-27 g2^2 
    +2 (
       + 20 k1+4 k2+8 l1+160 l2+32 l3+ 8 b1 
       + 3 (Tr[Yqp.Yqp]+Tr[Yqm.Yqm])
       + (Tr[Ylp.Ylp]+Tr[Ylm.Ylm])
       )
    +48 l4 b1 - 2 b1 (Tr[Ylp.Ylm]+3 Tr[Yqp.Yqm])
));

Beta\[Beta]3[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/(64 \[Pi]^2) * ( b3 (
      + 16 b1 - 6 g1^2 - 27 g2^2 + 8 k1 + 40 k2 
      + 16 (l1 - 4 l2) - 320 l3 
      +  6 (Tr[Yqp.Yqp]+ Tr[Yqm.Yqm]) + 2 ( Tr[Ylp.Ylp]+Tr[Ylm.Ylm])
));

Betaf1[g1_,g2_,g3_,Yqp_,Yqm_,Ylp_,Ylm_, l1_, l2_, l3_, l4_,k1_,k2_,b1_, b2_, b3_, f1_] := 
1/(64\[Pi]^2) * ( f1 (
      + 16 b1 - 6 g1^2 - 27 g2^2 + 8 k1 + 8 k2 
      + 16 (l1 - 4 l2) + 64 l3 
      +  6 (Tr[Yqp.Yqp]+ Tr[Yqm.Yqm]) + 2 ( Tr[Ylp.Ylp]+Tr[Ylm.Ylm])
));

    
ClearAll[Trans2HDMDirac];
Trans2HDMDirac[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{lg1,lg2,lg3,lYqp,lYqm,lYlp,lYlm,ll1,ll2,ll3,ll4,lk1,lk2,lb1,lb2,lb3,lf1},
(* exceptions: try to add new particles --> CanNotAddNewParticles
*)

	Print["transition function has to be implemented"];
(* calculate the new parameters *)
	{lg1,lg2,lg3,lYqp,lYqm,lYlp,lYlm,ll1,ll2,ll3,ll4,lk1,lk2,lb1,lb2,lb3,lf1}=(ParametersFunc[ pScale ]/.pSolution)[[1]];
        Return[{RGEg1->lg1,RGEg2->lg2,RGEg3->lg3,RGEYu->lYu,RGEYd->lYd,RGEYe->lYe,RGE\[Kappa]->l\[Kappa],RGE\[Lambda]->l\[Lambda]}];
];

(* internal functions *)

ClearAll[ParametersFunc];
ParametersFunc[pScale_]:={g1[pScale],g2[pScale],g3[pScale],Yqp[pScale],Yqm[pScale],Ylp[pScale],Ylm[pScale], \[Lambda]1[pScale],\[Lambda]2[pScale], \[Lambda]3[pScale], \[Lambda]4[pScale],\[Kappa]1[pScale],\[Kappa]2[pScale],\[Beta]1[pScale],\[Beta]2[pScale],\[Beta]3[pScale],f1[pScale]};



ClearAll[SetInitial];
SetInitial[pBoundary_?NumericQ,pInitial_]:=Block[{},
(* sets the initial values *)
   Return[	{g1[pBoundary]==RGEg1
   		,g2[pBoundary]==RGEg2
		,g3[pBoundary]==RGEg3
		,Yqp[pBoundary]==RGEYqp
		,Yqm[pBoundary]==RGEYqm
		,Ylp[pBoundary]==RGEYlp
		,Ylm[pBoundary]==RGEYlm
		, \[Lambda]1[pBoundary]==RGE\[Lambda]1
		,\[Lambda]2[pBoundary]==RGE\[Lambda]2
		, \[Lambda]3[pBoundary]==RGE\[Lambda]3
		, \[Lambda]4[pBoundary]==RGE\[Lambda]4
		,\[Kappa]1[pBoundary]==RGE\[Kappa]1
		,\[Kappa]2[pBoundary]==RGE\[Kappa]2
		,\[Beta]1[pBoundary]==RGE\[Beta]1
		,\[Beta]2[pBoundary]==RGE\[Beta]2
		,\[Beta]3[pBoundary]==RGE\[Beta]3
		,f1[pBoundary]==RGEf1
}//.pInitial];
];

End[]; (* end of `Private` *)


EndPackage[];
