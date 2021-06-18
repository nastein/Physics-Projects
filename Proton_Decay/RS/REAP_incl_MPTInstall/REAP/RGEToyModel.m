(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEToyModel`",{"REAP`RGESolver`","REAP`RGESymbol`","REAP`RGEUtilities`"}];

(* register several toy models *)
RGERegisterModel["Toy","REAP`Toy`",
	`Private`GetParameters,
	`Private`SolveModel,
	{RGEAll->`Private`GetSolution},
	{{"Toy",`Private`Transition}},
        `Private`GetInitial,
        `Private`ModelSetOptions,
        `Private`ModelGetOptions
];

Begin["`Private`"];
Map[Needs,{"REAP`RGESolver`","REAP`RGESymbol`","REAP`RGEUtilities`"}];

(* toy model *)

(* options *)
ClearAll[RGEOptions];
RGEOptions;
Options[RGEOptions]={
                       Method->StiffnessSwitching
};

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

(* parameters *)
Parameters={\[Lambda]};

ClearAll[GetParameters];
GetParameters[]:= Block[{},
(* returns the parameters of the model *)
   Return[Parameters];
];


(* initial *)
ClearAll[SetInitial]
SetInitial[pBoundary_?NumericQ,pInitial_]:={
	\[Lambda][pBoundary]==(RGE\[Lambda]/.pInitial)
	};

ClearAll[GetInitial];
GetInitial[pOpts___:{}]:=Block[{lInitial,l\[Lambda]},
	lInitial={RGE\[Lambda]->0.1};
	l\[Lambda]=(RGE\[Lambda]/.pOpts)/.lInitial;
	Return[{RGE\[Lambda]->l\[Lambda]}];
];


(* return solution *)
ClearAll[GetSolution];
GetSolution[pScale_,pSolution_,pOpts___]:=Block[{},
	Return[({\[Lambda][pScale]}/.pSolution)[[1]]];
];


(* RGE *)
ClearAll[RGE];
RGE:={	D[\[Lambda][t],t]==Beta\[Lambda][\[Lambda][t]] };

ClearAll[Beta\[Lambda]];
Beta\[Lambda][\[Lambda]_] :=-7 * 1/(16*Pi^2) * \[Lambda]^3;


(* transition functions *)
ClearAll[Transition];
Transition[pScale_?NumericQ,pDirection_?NumericQ,pSolution_,pToOpts_,pFromOpts_]:=Block[{},
       Return[({RGE\[Lambda]->\[Lambda][pScale]}/.pSolution)[[1]]];
];

(* solve model *)
ClearAll[SolveModel];
SolveModel[{pUp_,pUpModel_,pUpOptions_},{pDown_,pDownModel_,pDownOptions_},pDirection_,pBoundary_,pInitial_,pOpts___]:=Block[{lSolution,lNewScale,lInitial},
	lNDSolveOpt;
	Options[lNDSolveOpts]=Options[NDSolve];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,Options[RGEOptions]]];
	SetOptions[lNDSolveOpts,RGEFilterOptions[NDSolve,pOpts]];
	lInitial=SetInitial[pBoundary,pInitial];
	lSolution=NDSolve[RGE ~Join~ lInitial, Parameters,{t,pDown,pUp}, Sequence[Options[lNDSolveOpts]]];
	If[lDirection>0,lNewScale=pUp,lNewScale=pDown];
	Return[{lSolution,lNewScale,0}];
];


End[]; (* end of toy model *)

EndPackage[];
