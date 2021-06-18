(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)


BeginPackage["REAP`RGEUtilities`",{"REAP`RGESolver`","REAP`RGESymbol`"}];


ClearAll[RGEFilterOptions];
If[$VersionNumber < 6.,
  Needs["Utilities`FilterOptions`"];
  RGEFilterOptions[what_,list___]:=FilterOptions[what,list],
  RGEFilterOptions[what_, list_List] := Sequence @@ FilterRules[list, Options@what];
  RGEFilterOptions[what_, list___] := Sequence @@ FilterRules[{list}, Options@what];
  ];

ClearAll[RGEIntegrateOutM];

RGEIntegrateOutM::usage="RGEIntegrateOutM[M,Num] integrates out the heaviest <Num> degrees of freedom of M.";


ClearAll[RGEIntegrateOutY\[Nu]];

RGEIntegrateOutY\[Nu]::usage="RGEIntegrateOutY\[Nu][Y\[Nu],Num] integrates out <Num> degrees of freedom in Y\[Nu]";


ClearAll[RGETestMY\[Nu]];

RGETestMY\[Nu]::usage="RGETestMY\[Nu][scale,M,Y\[Nu],\[Kappa]] integrates out the degrees of freedom in M which are above scale and returns {M,Y\[Nu],\[Kappa],number of dof which are integrated out} ";


ClearAll[RGEKappaMatching];

RGEKappaMatching::usage="RGEKappaMatching[M,Y\[Nu],IntOut]  does the \[Kappa] matching at transitions where a right-handed neutrino is integrated out. The matrices have to be in a basis where M is diagonal.";


ClearAll[RGEGetNeutrinoMasses];

RGEGetNeutrinoMasses::usage="RGEGetNeutrinoMasses[MassHierarchy,\[CapitalDelta]m2atm,\[CapitalDelta]m2sol,Mlightest]
converts its arguments into a list containing the neutrino mass eigenvalues.";


ClearAll[RGERotateM];

RGERotateM::usage="RGERotateM[M,Y] changes the basis to the eigenstates of M and returns M and Y in that basis";


ClearAll[RGESearchTransitions];

RGESearchTransitions::usage="RGESearchTransitions[Mass,LogScale,MaxLogScale,Min,Options] returns a list of transitions which are found by integrating out degrees of freedom. Mass is a function returning the mass matrix at a given scale and LogScale, MaxLogScale and MinLogScale are the starting point, the maximum and the minimum respectively.";
Options[RGESearchTransitions]={
    RGEPrecision->6,
    RGEMaxNumberIterations->20,
    RGEThresholdFactor->1
};


ClearAll[RGEFloor];

RGEFloor::usage="RGEFloor[value,opts] returns <value> rounded to RGEPrecision which is given as an option <opts>.";
Options[RGEFloor]={RGEPrecision->6};


ClearAll[RGEGetRightHanded\[Nu]Masses];

RGEGetRightHanded\[Nu]Masses::usage="RGEGetRightHanded\[Nu]Masses[Scale] returns the
right-handed neutrino masses in a list, which is ordered by increasing
mass. <Scale> is the scale at which the returned right-handed neutrino masses
are defined.";



 Begin["`Private`"];

 Map[Needs,{"REAP`RGESolver`","REAP`RGESymbol`"}];
(*shortcuts*)
ClearAll[Dagger];
Dagger[x_] := Transpose[Conjugate[x]];


RGEIntegrateOutM[pM_,pNum_]:=Block[{lf,lg,lLen},
(* Integrates out the number of degrees of freedom in pM given in pNum *)
   lLen=Length[pM]-pNum;
   Return[Table[pM[[lf,lg]],{lf,lLen},{lg,lLen}]];
];


RGEIntegrateOutY\[Nu][pY\[Nu]_,pNum_]:=Block[{lf,lg,lx,ly},
(* Integrates out the number of degrees of freedom in pY\[Nu] given in pNum *)
   {lRow,lColumn}=Dimensions[pY\[Nu]];
   Return[Table[pY\[Nu][[lf,lg]],{lf,lRow-pNum},{lg,lColumn}]];
];


RGETestMY\[Nu][pScale_,pM_,pY\[Nu]_,p\[Kappa]_]:=Block[{lM\[Nu]Rotated,lY\[Nu]Rotated,lIntegrateOut,lY\[Nu],lM,l\[Kappa],lEvalues,lLogScale},
	l\[Kappa]=p\[Kappa];
	lY\[Nu]=pY\[Nu];
	lM=pM;
	lLogScale=Log[pScale];
	lEvalues=Map[Log[Abs[#]]/2&,Eigenvalues[N[Dagger[lM].lM]]];
	lIntegrateOut=Length[Select[lEvalues,(#1>lLogScale)&]];
	If[lIntegrateOut>0,
		{lM\[Nu]Rotated,lY\[Nu]Rotated}=RGERotateM[ lM,lY\[Nu] ];
		l\[Kappa]+= RGEKappaMatching[lM\[Nu]Rotated,lY\[Nu]Rotated,lIntegrateOut];
		lY\[Nu]= RGEIntegrateOutY\[Nu][lY\[Nu]Rotated, lIntegrateOut];
		lM= RGEIntegrateOutM[lM\[Nu]Rotated, lIntegrateOut];
	];

	Return[{lM,lY\[Nu],l\[Kappa],lIntegrateOut}];
];



RGEKappaMatching[pM_,pY\[Nu]_,pNum_]:=Block[{li,lIntOut},
(* does the \[Kappa] matching at the transitions *)
	li=Length[pM];
        lIntOut=Range[li-pNum+1,li];

	If[MatrixConditionNumber[pM[[lIntOut,lIntOut]]]>2*Precision[pM[[lIntOut,lIntOut]]],
		Print["RGEKappaMatching: The matrix M=", MatrixForm[ pM[[lIntOut,lIntOut]] ]," is ill-conditioned and the condition number is ",MatrixConditionNumber[ pM[[lIntOut,lIntOut]] ] ]];
	Return[2*Transpose[pY\[Nu][[lIntOut]]].Inverse[pM[[lIntOut,lIntOut]]].pY\[Nu][[lIntOut]]];
];


RGEGetNeutrinoMasses[pMassHierarchy_,p\[CapitalDelta]atm_,p\[CapitalDelta]sol_,pMlightest_]:=Block[{lM1,lM2,lM3,lSolM2},
    Which[pMassHierarchy=="i",
             lM3 = pMlightest;
	     lM2 = lM3^2 + p\[CapitalDelta]atm;
	     lM1 = Sqrt[lM2 - p\[CapitalDelta]sol];
	     lM2 = Sqrt[lM2],
           pMassHierarchy=="r"||pMassHierarchy=="n",
             lM1 = pMlightest;
	     lM2 = lM1^2 + p\[CapitalDelta]sol;
	     lM3 = Sqrt[lM2 + p\[CapitalDelta]atm];
	     lM2 = Sqrt[lM2],
           True, Throw[pMassHierarchy,RGENotAValidMassHierarchy]
    ];
    Return[{lM1,lM2,lM3}];
];


RGERotateM[pM_,pY_] := Block[{lU, lV, lEvalues, lEvectors},
(* returns pM and pY in the basis where pM is diagonal *)
    {lEvalues, lEvectors} = Eigensystem[N[ConjugateTranspose[pM].pM]];
    lU = Transpose[Reverse[lEvectors]];
    lV = lU.DiagonalMatrix[Exp[-\[ImaginaryI] Arg[Diagonal[Transpose[lU].pM.lU]]/2]];
    Return[{DiagonalMatrix[Sqrt[Reverse[lEvalues]]],Transpose[lV].pY}];
];


RGESearchTransitions[pMass_,pLogScale_?NumericQ,pMaxLogScale_?NumericQ,pMinLogScale_?NumericQ,pOpts___]:=Block[{
lTransitions, (* list of found transitions *)
lNewLogScale, (* used to search transitions *)
lOldLogScale, (* used to search transitions *)
lPrecision, (* precision used in searching transitions *)
lCount, (* number of iterations *)
lCountMax, (* maximum number of iterations *)
lLengthM, (* number of particles *)
lM, (* mass matrix *)
lAbsEVal, (* absolut value of the eigenvalues of neutrino mass matrix *)
li, (* variable used in loops *)
lLogThresholdFactor, (* neutrinos are integrated out below their mass; this is the shift *)
lBelowMinimum, (* Signal to indicate if the minimum is reached *)
lMinLogScale (* lower bound *)
},
(* search for transitions *)
(* exceptions: mass out of range --> RGEOutOfRange
*)
    lBelowMinimum;
    lMinLogScale=pMinLogScale;
    lTransitions={};
    lPrecision=RGEPrecision/.{RGEFilterOptions[RGESearchTransitions,pOpts]}/.Options[RGESearchTransitions,RGEPrecision];
    lCountMax=RGEMaxNumberIterations/.{RGEFilterOptions[RGESearchTransitions,pOpts]}/.Options[RGESearchTransitions,RGEMaxNumberIterations];
    lOldLogScale=pLogScale;
    lM=pMass[lOldLogScale];
    lLengthM=Length[lM];
    lAbsEVal=Abs[Eigenvalues[Dagger[lM].lM]];
    lNewLogScale=Log[Max[lAbsEVal]]/2;
    Catch[
    If[lNewLogScale>pMaxLogScale,
       lNewLogScale=RGEFloor[lNewLogScale,RGEPrecision->lPrecision];
       If[lNewLogScale>pMaxLogScale, Throw[{lNewLogScale,pMaxLogScale},RGEOutOfRange]];
       ];
    If[lNewLogScale<lMinLogScale,Throw[lNewLogScale,lBelowMinimum]];
       (* search transitions *)
    For[li=1,li<=lLengthM,li++,
         lCount=0;
         lOldLogScale=lNewLogScale;
         lM=pMass[lOldLogScale];
         lLengthM=Length[lM];
         lAbsEVal=Abs[Eigenvalues[Dagger[lM].lM]];
         lNewLogScale=Log[Sort[lAbsEVal,Greater][[li]]]/2;
         If[lNewLogScale>pMaxLogScale,
            lNewLogScale=RGEFloor[lNewLogScale,RGEPrecision->lPrecision];
            If[lNewLogScale>pMaxLogScale,Throw[0,RGEOutOfRange]];
         ];
         If[lNewLogScale<lMinLogScale,Throw[lNewLogScale,lBelowMinimum]];
         While[N[Abs[lOldLogScale-lNewLogScale],lPrecision]>0,
            lOldLogScale=lNewLogScale;
            lM=pMass[lOldLogScale];
            lAbsEVal=Abs[Eigenvalues[Dagger[lM].lM]];
            lNewLogScale=Log[Sort[lAbsEVal,Greater][[li]]]/2;
            lCount++;
            If[lNewLogScale>pMaxLogScale,
               lNewLogScale=RGEFloor[lNewLogScale,RGEPrecision->lPrecision];
               If[lNewLogScale>pMaxLogScale,Throw[0,RGEOutOfRange]];
            ];
            If[lNewLogScale<lMinLogScale,Throw[lNewLogScale,lBelowMinimum]];
            If[lCount>lCountMax,
                Print["RGESearchTransitions: algorithm to search transitions does not converge. There have been ",lCount," iterations so far. Returning: ",N[Sort[lTransitions,Greater],lPrecision]];
                Return[N[Sort[lTransitions,Greater],lPrecision]];
(*               Throw[N[Sort[lTransitions,Greater],lPrecision],RGETooManyIterations]];*)
            ]];
         lTransitions=Append[lTransitions,lNewLogScale];
         If[li==1,
              lMinLogScale=First[lTransitions]+Log[RGEThresholdFactor]/.{RGEFilterOptions[RGESearchTransitions,pOpts]}/.Options[RGESearchTransitions,RGEThresholdFactor]
         ];
   ],
   lBelowMinimum];
   Return[N[Sort[lTransitions,Greater],lPrecision]];
];


RGEFloor[pValue_,pOpts___]:=Block[{lValue,lPrecision},
    {lMantissa,lExponent}=MantissaExponent[pValue,10];
    lPrecision=RGEPrecision/.{RGEFilterOptions[RGEFloor,pOpts]}/.Options[RGEFloor,RGEPrecision];
    Return[N[lMantissa-10^-{lPrecision},lPrecision]*10^lExponent];
];


RGEGetRightHanded\[Nu]Masses[pScale_]:=Block[{lM},
	lM=RGEGetSolution[pScale,RGEM\[Nu]r];
	Return[Sort[Sqrt[Eigenvalues[Dagger[lM].lM]], Less]];
];


 End[];

 Protect[RGEIntegrateOutM,RGEIntegrateOutY\[Nu],RGETestMY\[Nu],RGEKappaMatching,RGEGetNeutrinoMasses,RGERotateM,RGESearchTransitions,RGEFloor,RGEGetRightHanded\[Nu]Masses];

 EndPackage[];
