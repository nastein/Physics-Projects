(* The package `REAP' is written for Mathematica 7 and is distributed under the
terms of GNU Public License http://www.gnu.org/copyleft/gpl.html *)




BeginPackage["REAP`RGEInitial`",{"MixingParameterTools`MPT3x3`","REAP`RGEUtilities`"}];

ClearAll[RGEGetYe];

RGEGetYe::usage="RGEGetYe[YukawaTauGUT] returns Ye at the GUT scale where Ye[[3,3]]=<YukawaTauGUT>";


ClearAll[RGEGetY\[Nu]];

RGEGetY\[Nu]::usage="RGEGetY\[Nu][Y\[Nu]33,Y\[Nu]Ratio] returns Y\[Nu] at the GUT scale where Y\[Nu]=<Y\[Nu]33>*MatrixForm[DiagonalMatrix[{<Y\[Nu]Ratio>^2,<Y\[Nu]Ratio>,1}]]";


ClearAll[RGEGetM];

RGEGetM::usage="RGEGetM[\[Theta]12,\[Theta]13,\[Theta]23,\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2,Mlightest,\[CapitalDelta]m2atm,\[CapitalDelta]m2sol,MassHierarchy,
vev of the Higgs coupling to the neutrinos  \!\(v\_u\) (GeV),Y\[Nu]] returns M\[Nu]r.";


ClearAll[RGEGet\[Kappa]];

RGEGet\[Kappa]::usage="RGEGet\[Kappa][\[Theta]12,\[Theta]13,\[Theta]23,\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2,Mlightest,\[CapitalDelta]m2atm,\[CapitalDelta]m2sol,MassHierarchy,
vev of the Higgs coupling
to the neutrinos \!\(v\_u\) (GeV)] returns \[Kappa] in \!\(GeV\^\(-1\)\).";


ClearAll[RGEGetDiracY\[Nu]];

RGEGetDiracY\[Nu]::usage="RGEGetDiracY\[Nu][\[Theta]12,\[Theta]13,\[Theta]23,\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2,Mlightest,\[CapitalDelta]m2atm,\[CapitalDelta]m2sol,MassHierarchy, vev of the Higgs coupling
to the neutrinos \!\(v\_u\) (GeV)] returns Y\[Nu].";


ClearAll[RGEGetYd];

RGEGetYd::usage="RGEGetYd[y1,y2,y3,\[Theta]12,\[Theta]13,\[Theta]23,\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2] returns Yd.";



 Begin["`Private`"];

 Map[Needs,{"MixingParameterTools`MPT3x3`","REAP`RGEUtilities`"}];
ClearAll[Dagger];
Dagger[x_]:=Transpose[Conjugate[x]];


RGEGetYe[pYukawaTauGUT_?NumericQ]:=pYukawaTauGUT*DiagonalMatrix[{1/3*0.001, 2/3*0.01, 1}];


RGEGetY\[Nu][p\[Nu]Factor_?NumericQ,pRatio_?NumericQ]:=Block[{lY\[Nu]},
    lY\[Nu]=p\[Nu]Factor*DiagonalMatrix[{pRatio^2,pRatio,1}];
    Return[lY\[Nu]];
];


RGEGetM[p\[Theta]12_?NumericQ,p\[Theta]13_?NumericQ,p\[Theta]23_?NumericQ,pDiracCP_?NumericQ,p\[Delta]e_?NumericQ,p\[Delta]\[Mu]_?NumericQ,p\[Delta]\[Tau]_?NumericQ,p\[Phi]1_?NumericQ,p\[Phi]2_?NumericQ,pMlightest_?NumericQ,p\[CapitalDelta]atm_?NumericQ,p\[CapitalDelta]sol_?NumericQ,pMassHierarchy_,pvu_?NumericQ,pY\[Nu]_?MPTNumericMatrixQ]:=Block[{l\[Kappa],lM\[Nu]r},
    l\[Kappa] = RGEGet\[Kappa][p\[Theta]12, p\[Theta]13, p\[Theta]23, pDiracCP, p\[Delta]e, p\[Delta]\[Mu], p\[Delta]\[Tau], p\[Phi]1, p\[Phi]2, pMlightest,p\[CapitalDelta]atm,p\[CapitalDelta]sol,pMassHierarchy,pvu];
    If[pMlightest==0,Print["Warning: One Neutrino is massless. Thus \[Kappa] is singular and can not be inverted. This results in large numerical errors. The returned right-handed neutrino mass matrix is calculated by \!\(M = 2*\((\(\((Y\_\[Nu]\^T)\)\^\(-1\)\) \[Kappa]\ \ Y\_\[Nu]\^\(-1\))\)\^\(-1\)\) instead of \!\(M = 2*\(Y\_\[Nu]\) \[Kappa]\^\(-1\)\ Y\_\[Nu]\^T\) which is used otherwise. Be careful in the interpretation of the results."];
			      lM\[Nu]r = 2*Inverse[Inverse[Transpose[pY\[Nu]]].l\[Kappa].Inverse[pY\[Nu]]],
			      lM\[Nu]r = 2*pY\[Nu].Inverse[l\[Kappa]].Transpose[pY\[Nu]]
    ];
			      
    Return[lM\[Nu]r];
];
RGEGetM[p\[Theta]12_?NumericQ,p\[Theta]13_?NumericQ,p\[Theta]23_?NumericQ,pDiracCP_?NumericQ,p\[Delta]e_?NumericQ,p\[Delta]\[Mu]_?NumericQ,p\[Delta]\[Tau]_?NumericQ,p\[Phi]1_?NumericQ,p\[Phi]2_?NumericQ,pMlightest_?NumericQ,p\[CapitalDelta]atm_?NumericQ,p\[CapitalDelta]sol_?NumericQ,pMassHierarchy_,pvu_?NumericQ,pY\[Nu]_?MPTNumericMatrixQ,p\[Kappa]_?MPTNumericMatrixQ,pY\[CapitalDelta]_?MPTNumericMatrixQ,p\[CapitalLambda]6_?NumericQ,pM2\[CapitalDelta]_?NumericQ]:=Block[{l\[Kappa],lM\[Nu]r},
    l\[Kappa] = RGEGet\[Kappa][p\[Theta]12, p\[Theta]13, p\[Theta]23, pDiracCP, p\[Delta]e, p\[Delta]\[Mu], p\[Delta]\[Tau], p\[Phi]1, p\[Phi]2, pMlightest,p\[CapitalDelta]atm,p\[CapitalDelta]sol,pMassHierarchy,pvu];
    l\[Kappa]=l\[Kappa]-(p\[Kappa]-2 p\[CapitalLambda]6/pM2\[CapitalDelta]*pY\[CapitalDelta]);

    If[pMlightest==0,Print["Warning: One Neutrino is massless. Thus \[Kappa] is singular and can not be inverted. This results in large numerical errors. The returned right-handed neutrino mass matrix is calculated by \!\(M = 2*\((\(\((Y\_\[Nu]\^T)\)\^\(-1\)\) \[Kappa]\ \ Y\_\[Nu]\^\(-1\))\)\^\(-1\)\) instead of \!\(M = 2*\(Y\_\[Nu]\) \[Kappa]\^\(-1\)\ Y\_\[Nu]\^T\) which is used otherwise. Be careful in the interpretation of the results."];
			      lM\[Nu]r = 2*Inverse[Inverse[Transpose[pY\[Nu]]].l\[Kappa].Inverse[pY\[Nu]]],
			      lM\[Nu]r = 2*pY\[Nu].Inverse[l\[Kappa]].Transpose[pY\[Nu]]
    ];
			      
    Return[lM\[Nu]r];
];


RGEGet\[Kappa][p\[Theta]12_?NumericQ,p\[Theta]13_?NumericQ,p\[Theta]23_?NumericQ,pDiracCP_?NumericQ,p\[Delta]e_?NumericQ,p\[Delta]\[Mu]_?NumericQ,p\[Delta]\[Tau]_?NumericQ,p\[Phi]1_?NumericQ,p\[Phi]2_?NumericQ,pMlightest_?NumericQ,p\[CapitalDelta]atm_?NumericQ,p\[CapitalDelta]sol_?NumericQ,pMassHierarchy_,pvu_?NumericQ]:=Block[{lV, lM, l\[Kappa]},
    lV = MPT3x3UnitaryMatrix[p\[Theta]12, p\[Theta]13, p\[Theta]23, pDiracCP, p\[Delta]e, p\[Delta]\[Mu], p\[Delta]\[Tau], p\[Phi]1, p\[Phi]2];

    lM = DiagonalMatrix[RGEGetNeutrinoMasses[pMassHierarchy,p\[CapitalDelta]atm,p\[CapitalDelta]sol,pMlightest]];
    l\[Kappa] = - 10^-9 *4/pvu^2*Conjugate[lV].lM.Dagger[lV];
    Return[l\[Kappa]];
];


RGEGetDiracY\[Nu][p\[Theta]12_?NumericQ,p\[Theta]13_?NumericQ,p\[Theta]23_?NumericQ,pDiracCP_?NumericQ,p\[Delta]e_?NumericQ,p\[Delta]\[Mu]_?NumericQ,p\[Delta]\[Tau]_?NumericQ,p\[Phi]1_?NumericQ,p\[Phi]2_?NumericQ,pMlightest_?NumericQ,p\[CapitalDelta]atm_?NumericQ,p\[CapitalDelta]sol_?NumericQ,pMassHierarchy_,pvu_?NumericQ]:=Block[{lV, lM, l\[Kappa]},
    lV = MPT3x3UnitaryMatrix[p\[Theta]12, p\[Theta]13, p\[Theta]23, pDiracCP, p\[Delta]e, p\[Delta]\[Mu], p\[Delta]\[Tau], p\[Phi]1, p\[Phi]2];

    lM = DiagonalMatrix[RGEGetNeutrinoMasses[pMassHierarchy,p\[CapitalDelta]atm,p\[CapitalDelta]sol,pMlightest]];
    lY\[Nu] = 10^-9 Sqrt[2]/pvu*lV.lM.Dagger[lV];
    Return[lY\[Nu]];
];


RGEGetYd[py1_?NumericQ,py2_?NumericQ,py3_?NumericQ,p\[Theta]12_?NumericQ,p\[Theta]13_?NumericQ,p\[Theta]23_?NumericQ,p\[Delta]_?NumericQ,p\[Delta]e_?NumericQ,p\[Delta]\[Mu]_?NumericQ,p\[Delta]\[Tau]_?NumericQ,p\[Phi]1_?NumericQ,p\[Phi]2_?NumericQ]:=Block[{lV},
    lV = MPT3x3UnitaryMatrix[p\[Theta]12, p\[Theta]13, p\[Theta]23, p\[Delta], p\[Delta]e, p\[Delta]\[Mu], p\[Delta]\[Tau], p\[Phi]1, p\[Phi]2];
    Return[lV.DiagonalMatrix[{py1,py2,py3}].Dagger[lV]];
];


 End[];

 Protect[RGEGetYe,RGEGetY\[Nu],RGEGetM,RGEGet\[Kappa],RGEGetDiracY\[Nu],RGEGetYd];

 EndPackage[];