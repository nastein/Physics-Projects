BeginPackage["MixingParameterTools`MPT3x3`"]

If[$VersionNumber < 6,
Needs["LinearAlgebra`MatrixManipulation`"];
Needs["LinearAlgebra`Orthogonalization`"];
];

ClearAll[MPT3x3OrthogonalMatrix];
MPT3x3OrthogonalMatrix::usage = "MPT3x3OrthogonalMatrix[\!\(\[Theta]\_12\),\!\(\[Theta]\_13\), \
\!\(\[Theta]\_23\),\[Delta]] \
returns the standard parametrized othogonal 3x3 matrix. \[Delta] can be \
omitted."

ClearAll[MPT3x3UnitaryMatrix];
MPT3x3UnitaryMatrix::usage = "MPT3x3UnitaryMatrix[\!\(\[Theta]\_12\),\!\(\[Theta]\_13\),\
\!\(\[Theta]\_23\),\[Delta],\!\(\[Delta]\_e\),\!\(\[Delta]\_\[Mu]\),\
\!\(\[Delta]\_\[Tau]\),\!\(\[Phi]\_1\),\!\(\[Phi]\_2\)] returns the \
unitary 3x3 matrix parametrized by \
{\!\(\[Theta]\_13\),\!\(\[Theta]\_12\),\!\(\[Theta]\_23\),\[Delta],\!\(\[Delt\
a]\_e\),\!\(\[Delta]\_\[Mu]\),\!\(\[Delta]\_\[Tau]\),\!\(\[Phi]\_1\),\
\!\(\[Phi]\_2\)} in standard parametrization."
 
ClearAll[MPT3x3MixingMatrixL];
MPT3x3MixingMatrixL::usage = "MPT3x3MixingMatrixL[M] \
returns the matrix \!\(U\_L\) which is used for diagonalizing a general \
complex matrix M, i.e. \
\!\(U\_R\^\[Dagger]\)\[CenterDot]M\[CenterDot]\!\(U\_L\)=diag(\!\(M\_1\),\!\(\
M\_2\),\!\(M\_3\)) and \
|\!\(M\_1\)|\[LessEqual]|\!\(M\_2\)|\[LessEqual]|\!\(M\_3\)|."

Options[MPT3x3MixingMatrixL] = {MPTTolerance -> 10^(-6)};
 
ClearAll[MPT3x3MixingMatrixR];
MPT3x3MixingMatrixR::usage = "MPT3x3MixingMatrixR[M] \
returns the matrix \!\(U\_R\) which is used for diagonalizing a general \
complex matrix M, i.e. \
\!\(U\_R\^\[Dagger]\)\[CenterDot]M\[CenterDot]\!\(U\_L\)=diag(\!\(M\_1\),\!\(\
M\_2\),\!\(M\_3\)) and \
|\!\(M\_1\)|\[LessEqual]|\!\(M\_2\)|\[LessEqual]|\!\(M\_3\)|."

Options[MPT3x3MixingMatrixR] = {MPTTolerance -> 10^(-6)};

ClearAll[MPT3x3NeutrinoMixingMatrix];
MPT3x3NeutrinoMixingMatrix::usage = "\
MPT3x3NeutrinoMixingMatrix[m] returns the matrix U which is used for diagonalizing a \
general complex symmetric Matrix m, i.e. \
\!\(U\^T\)\[CenterDot]m\[CenterDot]U=diag(\!\(m\_1\),\!\(m\_2\),\!\(m\_3\)) \
and \
|\!\(\[CapitalDelta]m\_32\^2\)|\[GreaterEqual]|\!\(\[CapitalDelta]m\_21\^2\)| \
and |\!\(U\_12\)|\[LessEqual]|\!\(U\_11\)|.\
ComplexNeutrinoMixingMatrix[M,S] \
does the same only that the mass hierarchy can be fixed to be inverted by \
setting S='i'. Any other S leads to a fixed regular mass hierarchy."

Options[MPT3x3NeutrinoMixingMatrix] = {MPTTolerance -> 10^(-6)};

ClearAll[MPT3x3MixingParameters];
MPT3x3MixingParameters::usage = "\
MPT3x3MixingParameters[U] returns the mixing angles and phases \
{\!\(\[Theta]\_12\),\!\(\[Theta]\_13\),\!\(\[Theta]\_23\),\[Delta],\!\(\[Delt\
a]\_e\),\!\(\[Delta]\_\[Mu]\),\!\(\[Delta]\_\[Tau]\),\!\(\[Phi]\_1\),\!\(\[Ph\
i]\_2\)} where 0\[LessEqual]\!\(\[Theta]\_ij\)\[LessEqual]\!\(\[Pi]\/2\) \
holds."

ClearAll[MNSParameters];
MNSParameters::usage = "MNSParameters[\!\(m\),\!\(Ye\)] returns the MNS mixing 
and mass \
parameters {{\!\(\[Theta]\_12\),\!\(\[Theta]\_13\),\!\(\[Theta]\_23\),\[Delta],\
\!\(\[Delta]\_e\),\!\(\[Delta]\_\[Mu]\),\!\(\[Delta]\_\[Tau]\),\!\(\[Phi]\_1\),\
\!\(\[Phi]\_2\)},{\!\(m\_1\),\!\(m\_2\),\!\(m\_3\)},{\!\(y\_e\),\!\(y\_\[Mu]\),
\!\(y\_\[Tau]\)}} for a Majorana neutrino matrix \!\(m\) and a Yukawa coupling \
matrix \!\(Ye\). The returned parameters obey the conventions \
0\[LessEqual]\!\(\[Theta]\_12\)\[LessEqual]\!\(\[Pi]\/4\), \
0\[LessEqual]\!\(\[Theta]\_13\),\!\(\[Theta]\_23\)\[LessEqual]\!\(\[Pi]\/2\), \
and all other parameters range from 0 to \!\(2\[Pi]\). Furthermore, the \
neutrino masses \!\(m\_i\) are such that \
\!\(|m\_3\^2-m\_2\^2|>|m\_2\^2-m\_1\^2|\). Note that the input \
matrices \!\(m\) and \!\(Ye\) have to be numeric."

ClearAll[MNSMatrix];
MNSMatrix::usage = "MNSMatrix[\!\(m\),\!\(Ye\)] returns the MNS mixing \
matrix, i.e. the matrix \!\(U\_MNS\) which diagonalizes \!\(m\) in the basis \
where \!\(Ye\) is diagonal. Note that the input matrices \!\(m\) and \!\(Ye\) \
have to be numeric."

ClearAll[DiracMNSMatrix];
DiracMNSMatrix::usage = "DiracMNSMatrix[\!\(Y\[Nu]\),\!\(Ye\)] returns the \
MNS matrix for Dirac neutrinos with Yukawa coupling \!\(Y\[Nu]\)."

ClearAll[DiracMNSParameters];
DiracMNSParameters::usage = "DiracMNSParameters[\!\(Y\[Nu]\),\!\(Ye\)] returns \
the MNS parameters {\!\(\[Theta]\_12\),\!\(\[Theta]\_13\),\
\!\(\[Theta]\_23\),\[Delta]}, {\!\(y\_1\),\!\(y\_2\),\!\(y\_3\)} \
and  {\!\(y\_e\),\!\(y\_\[Mu]\),\!\(y\_\[Tau]\)} \
for Dirac neutrinos."

ClearAll[CKMParameters];
CKMParameters::usage = "CKMParameters[\!\(Yu\),\!\(Yd\)] returns \
the CKM mixing parameters {\!\(\[Theta]\_12\),\!\(\[Theta]\_13\),\
\!\(\[Theta]\_23\),\[Delta]} for up- and down-type Yukawa matrices \
\!\(Yu\) and \!\(Yd\), as well as the Yukawa couplings \
{\!\(y\_u\),\!\(y\_c\),\!\(y\_t\)} and {\!\(y\_d\),\!\(y\_s\),\!\(y\_b\)}."

ClearAll[CKMMatrix];
CKMMatrix::usage = "CKMMatrix[\!\(Yu\),\!\(Yd\)] returns the CKM mixing \
matrix, i.e. the matrix \!\(U\_CKM\) which diagonalizes \!\(Yd\) in the basis \
where \!\(Yu\) is diagonal."

ClearAll[CKMReplacementRules];
CKMReplacementRules::usage = "CKMReplacementRules[\!\(Yu\),\!\(Yd\)] returns \
a list of replacement rules \
{MPT\!\(\[Theta]12\) -> \!\(\[Theta]\_12\),\
MPT\!\(\[Theta]13\) -> \!\(\[Theta]\_13\),\
MPT\!\(\[Theta]23\) -> \!\(\[Theta]\_23\),\
MPT\!\(\[Delta]\) -> \!\(\[Delta]\),\
MPT\!\(yu\) -> \!\(y\_u\), MPT\!\(yc\) -> \!\(y\_c\), MPT\!\(yt\) -> \!\(y\_t\),\
MPT\!\(yd\) -> \!\(y\_d\), MPT\!\(ys\) -> \!\(y\_s\), MPT\!\(yb\) -> \!\(y\_b\)},
suitable to replace the parameters left of `->' by the numerical value \
calculated from the input matrices \!\(Yu\) and \!\(Yd\)."

ClearAll[MPTNumericMatrixQ];
MPTNumericMatrixQ::usage = "MPTNumericMatrixQ[M] returns \
True if M is a numeric matrix."

ClearAll[MPTDebug];
MPTDebug::usage = "MPTDebug : for testing only."

ClearAll[MPT\[Theta]12,MPT\[Theta]13,MPT\[Theta]23,MPT\[Delta],
MPTyu,MPTyc,MPTyt,MPTyd,MPTys,MPTyb];


Begin["`Private`"];

ClearAll[MPTDebugPrint,MPTSecureArg,MPTAvoidAmbiguity,NumericMatrixQ,MPT3x3InvSort];

MPTDebugPrint[Expr__] := If[MPTDebug == True, Print[Expr]]

MPTSecureArg[x_] := Which[x != 0, Arg[x], x == 0, 0]

MPTAvoidAmbiguity[x_] = 2 * Pi * (x/(2*Pi)-Floor[x/(2*Pi)]);

MPTNumericMatrixQ[x_] := MatrixQ[x, NumericQ];

MPT3x3InvSort[m_?MPTNumericMatrixQ]:=Block[{v1,v2,v3},
	v1 = Transpose[m][[3]];
	v2 = Transpose[m][[2]]; 
	v3 = Transpose[m][[1]];
	Return[Transpose[{v1, v2, v3}]]];



MPT3x3OrthogonalMatrix[\[Theta]12_, \[Theta]13_, \[Theta]23_] := 
    Block[{CKMt, \[Delta]}, \[Delta] = 0; 
      CKMt = {{Cos[\[Theta]13]*Cos[\[Theta]12], Sin[\[Theta]12]*
          Cos[\[Theta]13], Sin[\[Theta]13]*Exp[(-I)*\[Delta]]}, 
        {(-Sin[\[Theta]12])*Cos[\[Theta]23] - Cos[\[Theta]12]*Sin[\[Theta]23]*
           Sin[\[Theta]13]*Exp[I*\[Delta]], Cos[\[Theta]12]*Cos[\[Theta]23] - 
          Sin[\[Theta]13]*Sin[\[Theta]23]*Sin[\[Theta]12]*Exp[I*\[Delta]], 
         Sin[\[Theta]23]*Cos[\[Theta]13]}, {Sin[\[Theta]12]*Sin[\[Theta]23] - 
          Cos[\[Theta]12]*Cos[\[Theta]23]*Sin[\[Theta]13]*Exp[I*\[Delta]], 
         (-Cos[\[Theta]12])*Sin[\[Theta]23] - Sin[\[Theta]12]*Cos[\[Theta]23]*
           Sin[\[Theta]13]*Exp[I*\[Delta]], Cos[\[Theta]23]*
          Cos[\[Theta]13]}}; CKMt]
 
MPT3x3OrthogonalMatrix[\[Theta]12_, \[Theta]13_, \[Theta]23_, \[Delta]_] := 
    Block[{CKMt}, CKMt = {{Cos[\[Theta]13]*Cos[\[Theta]12], 
         Sin[\[Theta]12]*Cos[\[Theta]13], Sin[\[Theta]13]*
          Exp[(-I)*\[Delta]]}, {(-Sin[\[Theta]12])*Cos[\[Theta]23] - 
          Cos[\[Theta]12]*Sin[\[Theta]23]*Sin[\[Theta]13]*Exp[I*\[Delta]], 
         Cos[\[Theta]12]*Cos[\[Theta]23] - Sin[\[Theta]13]*Sin[\[Theta]23]*
           Sin[\[Theta]12]*Exp[I*\[Delta]], Sin[\[Theta]23]*Cos[\[Theta]13]}, 
        {Sin[\[Theta]12]*Sin[\[Theta]23] - Cos[\[Theta]12]*Cos[\[Theta]23]*
           Sin[\[Theta]13]*Exp[I*\[Delta]], (-Cos[\[Theta]12])*
           Sin[\[Theta]23] - Sin[\[Theta]12]*Cos[\[Theta]23]*Sin[\[Theta]13]*
           Exp[I*\[Delta]], Cos[\[Theta]23]*Cos[\[Theta]13]}}; CKMt]
 
MPT3x3UnitaryMatrix[\[Theta]12_, \[Theta]13_, \[Theta]23_, \[Delta]_, 
     \[Delta]e_, \[Delta]\[Mu]_, \[Delta]\[Tau]_, \[Phi]1_, \[Phi]2_] := 
    Block[{tmpU}, tmpU = DiagonalMatrix[{Exp[I*\[Delta]e], 
          Exp[I*\[Delta]\[Mu]], Exp[I*\[Delta]\[Tau]]}] . 
        MPT3x3OrthogonalMatrix[\[Theta]12, \[Theta]13, \[Theta]23, 
         \[Delta]] . DiagonalMatrix[{Exp[(-I)*(\[Phi]1/2)], 
          Exp[(-I)*(\[Phi]2/2)], 1}]; Return[tmpU]]
 
MPT3x3UnitaryMatrix[MA_List] := Block[{tmpList}, 
     tmpList = Join[MA, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}]; 
      Return[MPT3x3UnitaryMatrix[tmpList[[1]], tmpList[[2]], tmpList[[3]], 
        tmpList[[4]], tmpList[[5]], tmpList[[6]], tmpList[[7]], tmpList[[8]], 
        tmpList[[9]]]]]


MPT3x3MixingMatrixL[M_?MPTNumericMatrixQ,pOpts___Rule] := 
	Block[{tmpU, evalues, evectors, iLauf, iLen, 
      tmpDev, tmpV, tmpPrecision, tmpVec, tmpTolerance, tmpSVD},
	 tmpTolerance = MPTTolerance /. {pOpts} /. Options[MPT3x3MixingMatrixL,
	  MPTTolerance];
     tmpPrecision = $MachinePrecision + 10; 
	 iLen = Length[M[[1]]]; 
	 If[iLen != 3, Throw["MPT3x3MixingMatrixL::Encountered non-3x3 matrix!"]];
	 evalues = Sort[SingularValueList[N[M],Tolerance->0]];
	 tmpU = MPT3x3InvSort[SingularValueDecomposition[N[M]][[3]]];
	 tmpV = MPT3x3InvSort[SingularValueDecomposition[N[M]][[1]]];
	 tmpDev = Abs[MatrixNorm[Conjugate[Transpose[tmpV]].M.tmpU
	 	- DiagonalMatrix[evalues]]];
     If[tmpDev > Abs[tmpTolerance * Tr[DiagonalMatrix[evalues]]], 
	  Throw["MPT3x3MixingMatrixL::\!\(U\_L\) could not be determined!"]]; 
	 Return[tmpU]];

MPT3x3MixingMatrixR[M_?MPTNumericMatrixQ,pOpts___Rule] := 
	Block[{tmpU, evalues, evectors, iLauf, iLen, 
      tmpDev, tmpV, tmpPrecision, tmpVec,tmpTolerance},
	 tmpTolerance = MPTTolerance /. {pOpts} /. Options[MPT3x3MixingMatrixR,
	  MPTTolerance];
     tmpPrecision = $MachinePrecision + 10; 
	 iLen = Length[M[[1]]]; 
	 If[iLen != 3, Throw["MPT3x3MixingMatrixR::Encountered non-3x3 matrix!"]];
	 evalues = Sort[SingularValueList[N[M],Tolerance->0]];
	 tmpU = MPT3x3InvSort[SingularValueDecomposition[N[M]][[3]]];
	 tmpV = MPT3x3InvSort[SingularValueDecomposition[N[M]][[1]]];
	 tmpDev = Abs[MatrixNorm[Conjugate[Transpose[tmpV]].M.tmpU
	 	- DiagonalMatrix[evalues]]];
     If[tmpDev > Abs[tmpTolerance * Tr[DiagonalMatrix[evalues]]], 
	  Throw["MPT3x3MixingMatrixL::\!\(U\_L\) could not be determined!"]]; 
	 Return[tmpV]];
 
MPT3x3NeutrinoMixingMatrix[M_?MPTNumericMatrixQ,pOpts___Rule] := 
	Block[{tmpU, tmpVec, evalues, evectors, iLauf, jLauf, iLen, 
	 tmpV, tmp\[Phi]1, tmp\[Phi]2, tmp\[Phi]3, tmpTolerance, tmpSDV},
	 tmpTolerance = MPTTolerance /. {pOpts} /. 
	 	Options[MPT3x3NeutrinoMixingMatrix, MPTTolerance]; 
     If[MatrixNorm[M - Transpose[M]] > tmpTolerance, 
       Throw["MPT3x3NeutrinoMixingMatrix::Encountered non-symmetric matrix!"]]; 
	 iLen = Length[M[[1]]]; 
	 If[iLen != 3, Throw["MPT3x3NeutrinoMixingMatrix::Encountered non-3x3 matrix!"]];
	 tmpU = {}; 
     {evalues, evectors} = Eigensystem[N[Conjugate[Transpose[M]] . M]]; 
	 MPTDebugPrint[Eigensystem[N[Conjugate[Transpose[M]] . M]]];
     evalues = Sqrt[evalues];	 
	 For[iLauf = iLen, iLauf >= 1, 
       tmpVec = evectors[[iLauf]]; 
	   tmpVec = tmpVec/Abs[Sqrt[Conjugate[tmpVec] . tmpVec]]; 
       If[iLauf < iLen, tmpVec = evectors[[iLauf]] - 
           Sum[Conjugate[evectors[[jLauf]]] . evectors[[iLauf]]*
             evectors[[jLauf]], {jLauf, 1, iLauf - 1}]]; 
       tmpVec = tmpVec/Abs[Sqrt[Conjugate[tmpVec] . tmpVec]]; 
       evectors[[iLauf]] = tmpVec; tmpU = Append[tmpU, {evalues[[iLauf]], 
           tmpVec}]; 
	   MPTDebugPrint[iLauf, tmpVec, ",", evalues[[iLauf]]]; 
      iLauf--];	 
	 tmpU = Sort[tmpU, Abs[#1[[1]]] < Abs[#2[[1]]] & ]; 
     If[Abs[Abs[tmpU[[2]][[1]]]^2 - Abs[tmpU[[1]][[1]]]^2] > 
        Abs[Abs[tmpU[[3]][[1]]]^2 - Abs[tmpU[[2]][[1]]]^2], 
      If[Abs[tmpU[[3]][[2]][[1]]] < Abs[tmpU[[2]][[2]][[1]]], 
         tmpU = {tmpU[[2]], tmpU[[3]], tmpU[[1]]}, 
         tmpU = {tmpU[[3]], tmpU[[2]], tmpU[[1]]}; ], 
      If[Abs[tmpU[[1]][[2]][[1]]] < Abs[tmpU[[2]][[2]][[1]]], 
         tmpU = {tmpU[[2]], tmpU[[1]], tmpU[[3]]}; ]; ]; 
	 For[iLauf = 1, iLauf <= iLen, 
      evectors[[iLauf]] = tmpU[[iLauf]][[2]]; iLauf++]; 
	 evalues = {tmpU[[1]][[1]], tmpU[[2]][[1]], tmpU[[3]][[1]]}; 
	 tmpU = Transpose[evectors]; 
	 tmpV = Transpose[tmpU] . M . tmpU; 
	 tmp\[Phi]1 = (-(1/2))*MPTSecureArg[tmpV[[1,1]]]; 
     tmp\[Phi]2 = (-(1/2))*MPTSecureArg[tmpV[[2,2]]]; 
     tmp\[Phi]3 = (-(1/2))*MPTSecureArg[tmpV[[3,3]]]; 
     tmpU = tmpU . DiagonalMatrix[{Exp[I*tmp\[Phi]1], Exp[I*tmp\[Phi]2], 
          Exp[I*tmp\[Phi]3]}]; 
     MPTDebugPrint["\!\(U\) = ",tmpU]; 
	 MPTDebugPrint["\!\(U\^T\)\[CenterDot]M\[CenterDot]U = ", 
      Transpose[tmpU] . M . tmpU]; 
	 MPTDebugPrint["Diag[Eigenvalues[M]] = ", Abs[DiagonalMatrix[evalues]]];
	 If[evalues[[1]] == evalues[[2]] || evalues[[2]] == evalues[[3]]
	  || evalues[[1]] == evalues[[3]],
	  Print["MPT3x3NeutrinoMixingMatrix::Warning: degenerate eigenvalues"]];
	 If[MatrixNorm[N[Transpose[tmpU] . M . tmpU - 
            Abs[DiagonalMatrix[evalues]]]] > tmpTolerance,
      Throw["MPT3x3NeutrinoMixingMatrix::U could not be determined"]];
	 Return[tmpU]];
 
MPT3x3NeutrinoMixingMatrix[M_?MPTNumericMatrixQ, S_,pOpts___Rule] := 
	Block[{tmpU, tmpVec, evalues, evectors, iLauf, jLauf, iLen, tmpV, 
	  tmp\[Phi]1, tmp\[Phi]2, tmp\[Phi]3, tmpTolerance, tmpSDV},
	 tmpTolerance = MPTTolerance /. {pOpts} /. 
	 	Options[MPT3x3NeutrinoMixingMatrix, MPTTolerance];  
	 If[MatrixNorm[M - Transpose[M]] > 0.01, 
       Throw["MPT3x3NeutrinoMixingMatrix::Encountered non-symmetric matrix!"]]; iLen = Length[M[[1]]]; tmpU = {}; 
     {evalues, evectors} = Eigensystem[N[Conjugate[Transpose[M]] . M]]; 
     evalues = Sqrt[evalues];
	 MPTDebugPrint[MatrixForm[evectors . M . Transpose[evectors]]]; 
	 MPTDebugPrint[MatrixForm[evectors . M . Transpose[evectors]]]; 
	 For[iLauf = iLen, iLauf >= 1, 
       tmpVec = evectors[[iLauf]]; tmpVec = 
         tmpVec/Abs[Sqrt[Conjugate[tmpVec] . tmpVec]]; 
        If[iLauf < iLen, tmpVec = evectors[[iLauf]] - 
           Sum[Conjugate[evectors[[jLauf]]] . evectors[[iLauf]]*
             evectors[[jLauf]], {jLauf, 1, iLauf - 1}]]; 
        tmpVec = tmpVec/Abs[Sqrt[Conjugate[tmpVec] . tmpVec]]; 
        evectors[[iLauf]] = tmpVec; 
		tmpU = Append[tmpU, {evalues[[iLauf]], tmpVec}]; 
		MPTDebugPrint[iLauf, tmpVec, ",", evalues[[iLauf]]]; 
        iLauf--]; tmpU = Sort[tmpU, Abs[#1[[1]]] < Abs[#2[[1]]] & ]; 
      MPTDebugPrint["After first sort ", tmpU]; 
	  If[S == "i", MPTDebugPrint["Encountered inverted mass hierarchy"]; 
         tmpU = {tmpU[[2]], tmpU[[3]], tmpU[[1]]}]; 
      MPTDebugPrint["After second sort ", tmpU]; 
	  MPTDebugPrint["After final sort ", tmpU]; 
	  For[iLauf = 1, iLauf <= iLen, 
       evectors[[iLauf]] = tmpU[[iLauf]][[2]]; 
	   MPTDebugPrint[iLauf, tmpU[[iLauf]]]; iLauf++]; 
	  evalues = {tmpU[[1]][[1]], tmpU[[2]][[1]], tmpU[[3]][[1]]}; 
	  tmpU = Transpose[evectors]; 
      MPTDebugPrint["U is now ", tmpU]; tmpV = Transpose[tmpU] . M . tmpU; 
      MPTDebugPrint["\!\(U\^T\)\[CenterDot]M\[CenterDot]U = ", 
       MatrixForm[tmpV]]; tmp\[Phi]1 = (-(1/2))*MPTSecureArg[tmpV[[1,1]]]; 
      tmp\[Phi]2 = (-(1/2))*MPTSecureArg[tmpV[[2,2]]]; 
      tmp\[Phi]3 = (-(1/2))*MPTSecureArg[tmpV[[3,3]]]; 
      tmpU = tmpU . DiagonalMatrix[{Exp[I*tmp\[Phi]1], Exp[I*tmp\[Phi]2], 
          Exp[I*tmp\[Phi]3]}]; 
	  MPTDebugPrint["The following should be real:"]; 
      MPTDebugPrint["\!\(U\^T\)\[CenterDot]M\[CenterDot]U = ", Transpose[tmpU] . M . tmpU]; 
	  MPTDebugPrint["Diag[Eigenvalues[M]] = ", Abs[DiagonalMatrix[evalues]]]; 
	  MPTDebugPrint["MatrixNorm of the difference = ", MatrixNorm[N[Transpose[tmpU] . M . tmpU - Abs[DiagonalMatrix[evalues]]]]]; 
      If[MatrixNorm[N[Transpose[tmpU] . M . tmpU - Abs[DiagonalMatrix[evalues]]]] > 0.01, 
        Throw["MPT3x3NeutrinoMixingMatrix::Mixing matrix could not be determined!"];
		MPTDebugPrint["Warning: Deviation ", MatrixNorm[
           N[Transpose[tmpU] . M . tmpU - DiagonalMatrix[evalues]]]]; 
          MPTDebugPrint["\!\(U\^T\)\[CenterDot]M\[CenterDot]U=", MatrixForm[N[Transpose[tmpU] . M . tmpU]]]; 
		  MPTDebugPrint["Eigenvalues: ", MatrixForm[DiagonalMatrix[evalues]]]; 
		  MPTDebugPrint["U= ", MatrixForm[N[tmpU]]]];
	  MPTDebugPrint[evalues]; 
	  Return[tmpU]];
  
MPT3x3MixingParameters[U_?MPTNumericMatrixQ] := Block[{\[Theta]13, \[Theta]12, \[Theta]23, 
      \[Delta], \[Delta]\[Mu], \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2, 
      \[CapitalDelta]\[Phi], tmpDev, tmpV, c, s}, 
	 If[Abs[U[[1,1]]] == 1 && Abs[U[[2,2]]] == 1 && Abs[U[[3,3]]] == 1,
	  Return[{0,0,0,0,Arg[U[[1,1]]],Arg[U[[2,2]]],Arg[U[[3,3]]],0,0}]];
     \[Delta]\[Mu] = MPTSecureArg[U[[2,3]]]; 
	 \[Delta]\[Tau] = MPTSecureArg[U[[3,3]]]; 
	 \[Theta]13 = ArcSin[Abs[U[[1]][[3]]]];
	 (*several cases for U_{13}=0 or \ne 0*)
     Which[
	  Abs[U[[1,3]]] != 1, 
	   Which[
	    Abs[U[[1,1]]] > Abs[U[[1,2]]] && Abs[U[[1,2]]] > 0, 
         \[Theta]12 = ArcTan[Abs[U[[1,2]]/U[[1,1]]]]; 
         Which[
		  U[[3]][[3]] != 0, 
		   \[Theta]23 = ArcTan[Abs[U[[2,3]]/U[[3,3]]]], 
		  U[[3]][[3]] == 0, 
           \[Theta]23 = N[Pi/2]], 
	    Abs[U[[1,1]]] > Abs[U[[1,2]]] && Abs[U[[1,2]]] == 0, 
         \[Theta]12 = 0; 
         Which[
		  U[[3,3]] != 0, 
		   \[Theta]23 = ArcTan[Abs[U[[2,3]]/U[[3,3]]]], 
		  U[[3,3]] == 0, 
           \[Theta]23 = N[Pi/2]], 
		Abs[U[[1,1]]] <= Abs[U[[1,2]]] && Abs[U[[1,1]]] > 0, 
         \[Theta]12 = ArcCot[Abs[U[[1,1]]/U[[1,2]]]]; 
         Which[
		  U[[3,3]] != 0, 
		   \[Theta]23 = ArcTan[Abs[U[[2,3]]/U[[3,3]]]], 
		  U[[3,3]] == 0, 
           \[Theta]23 = N[Pi/2]], 
		Abs[U[[1,1]]] == 0 && Abs[U[[1,2]]] > 0,
		 \[Theta]12 = N[Pi/2];
		 (*\[Theta]13 = ArcCos[Abs[U[[1,2]]]];*)
         Which[
		  U[[3]][[3]] != 0, 
		   \[Theta]23 = ArcTan[Abs[U[[2]][[3]]/U[[3]][[3]]]], 
		  U[[3]][[3]] == 0, 
           \[Theta]23 = N[Pi/2]],
		Abs[U[[1,1]]] > 0 && Abs[U[[1,2]]] == 0,
		 \[Theta]12 = 0;
		 (*\[Theta]13 = ArcSin[Abs[U[[1,1]]]];*)
		 Which[
		  U[[3]][[3]] != 0, 
		   \[Theta]23 = ArcTan[Abs[U[[2]][[3]]/U[[3]][[3]]]], 
		  U[[3]][[3]] == 0, 
           \[Theta]23 = N[Pi/2]],
		Abs[U[[1,1]]] == 0 && Abs[U[[1,2]]] == 0,
		 Throw["MPT3x3MixingParameters::\!\(\[Theta]\_13=\[Pi]/2\) encountered"];	 
		MPTDebugPrint[U];
		MPTDebugPrint[MPT3x3UnitaryMatrix[\[Theta]12,\[Theta]13,\[Theta]23,
		\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2]];
		], 
	  Abs[U[[1]][[3]]] == 1, 
	  	\[Delta] = 0;
			\[Delta]e = MPTSecureArg[U[[1,3]]];
			\[Theta]12 = 0;
			If[ U[[3,1]] != 0,
			  \[Theta]23 = ArcTan[Abs[U[[2,1]]/U[[3,1]]]];
			  If[ \[Theta]23 != 0,
			    \[Delta]\[Mu] = MPTSecureArg[-U[[2,1]]];
			    \[Phi]2 = 2*(\[Delta]\[Mu] - Arg[U[[2,2]]]);
			  ,
			    \[Delta]\[Mu] = MPTSecureArg[U[[2,2]]];
			    \[Phi]2 = 0;
			  ];
			  \[Delta]\[Tau] = MPTSecureArg[-U[[3,1]]];
			,
			  \[Theta]23 = N[Pi/2];
			  \[Delta]\[Mu] = MPTSecureArg[-U[[2,1]]];
			  \[Phi]2 = 0;
			  \[Delta]\[Tau] = MPTSecureArg[-U[[3,2]]];
			];
		\[Phi]1 = 0;
		Print["MPT3x3MixingParameters::Warning: \!\(\[Theta]\_13=\[Pi]/2\)
		encountered! Various parameters are not uniquely determined!"];
		If[MatrixNorm[U-MPT3x3UnitaryMatrix[\[Theta]12,\[Theta]13,\[Theta]23,
			\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2]]<0.1,
			Return[{\[Theta]12,\[Theta]13,\[Theta]23,
			\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2}],
			MPTDebugPrint["U=", MatrixForm[U]]; 
			MPTDebugPrint[{\[Theta]12,\[Theta]13,\[Theta]23,
			\[Delta],\[Delta]e,\[Delta]\[Mu],\[Delta]\[Tau],\[Phi]1,\[Phi]2}]; 
			MPTDebugPrint["Parameters reinserted: ", 
		 	MatrixForm[MPT3x3UnitaryMatrix[\[Theta]12, \[Theta]13, \[Theta]23, 
           \[Delta], \[Delta]e, \[Delta]\[Mu], \[Delta]\[Tau], \[Phi]1, 
           \[Phi]2]]]; 
		Throw["MPT3x3MixingParameters::\!\(\[Theta]\_13=\[Pi]/2\) encountered"];];	 
	  ]; 
	 (*end of several cases*)
     \[Delta] = Which[
	  Sin[2*\[Theta]12]*Sin[2*\[Theta]13]*Sin[2*\[Theta]23]*Cos[\[Theta]13] != 0, 
	  -MPTSecureArg[((Conjugate[U[[1,1]]]*U[[1,3]]*U[[3,1]]*Conjugate[U[[3,
                3]]])/(Cos[\[Theta]12]*Cos[\[Theta]13]^2*Cos[\[Theta]23]*
              Sin[\[Theta]13]) + Cos[\[Theta]12]*Cos[\[Theta]23]*
             Sin[\[Theta]13])/(Sin[\[Theta]12]*Sin[\[Theta]23])],
	  Sin[2*\[Theta]12]*Sin[2*\[Theta]13]*Sin[2*\[Theta]23]*Cos[\[Theta]13] == 0,
	   0]; 
	  Which[
  		Abs[U[[1,3]]] != 0 && Abs[U[[1,3]]] != 1,
    	 \[Delta]e = MPTSecureArg[Exp[I*\[Delta]]*U[[1,3]]];
    	 \[Phi]1 = 2*MPTSecureArg[Exp[I*\[Delta]e]*Conjugate[U[[1,1]]]];
    	 \[Phi]2 = 2*MPTSecureArg[Exp[I*\[Delta]e]*Conjugate[U[[1,2]]]];
	     If[Abs[U[[2,3]]] == 0,
	      \[Delta]\[Mu] = MPTSecureArg[Exp[I*\[Phi]2/2]*U[[2,2]]];
	     ];
	     If[Abs[U[[3,3]]] == 0,
	      \[Delta]\[Tau] = MPTSecureArg[-Exp[I*\[Phi]2/2]*U[[3,2]]];
	     ];
	     If[Abs[U[[1,2]]] == 0,
	      \[Phi]2 = 2*MPTSecureArg[Exp[I*\[Delta]\[Mu]]*Conjugate[U[[2,2]]]];
	     ];
	    ,
  		Abs[U[[1,3]]] == 0,
    	 \[Delta] = 0;
    	 \[Phi]2 = 2*MPTSecureArg[Exp[I*\[Delta]\[Mu]]*Conjugate[U[[2,2]]]];
    	 \[Delta]e = MPTSecureArg[Exp[I*\[Phi]2/2]*U[[1,2]]];
    	 \[Phi]1 = 2*MPTSecureArg[Exp[I*\[Delta]e]*Conjugate[U[[1,1]]]];
    	 If[Abs[U[[1,2]]] == 0,
      	  \[Phi]1 = 0;
      	  \[Delta]e = MPTSecureArg[U[[1,1]]];
         ];
         If[Abs[U[[2,3]]] == 0,
          \[Theta]23 = 0;
          \[Phi]1 = 0;
          \[Delta]e = MPTSecureArg[U[[1,1]]];
          \[Phi]2 = 2*MPTSecureArg[Exp[I*\[Delta]e]*Conjugate[U[[1,2]]]];
          \[Delta]\[Mu] = MPTSecureArg[-U[[2,1]]];
         ];
        If[Abs[U[[3,3]]] == 0,
         \[Theta]23 = Pi/2;
         \[Phi]1 = 0;
         \[Delta]e = MPTSecureArg[U[[1,1]]];
         \[Phi]2 = 2*MPTSecureArg[Exp[I*\[Delta]e]*Conjugate[U[[1,2]]]];
         \[Delta]\[Tau] = MPTSecureArg[-Exp[I*\[Phi]2/2]*U[[3,2]]];
        ];
	   ];
	  \[Theta]12 = MPTAvoidAmbiguity[\[Theta]12];
	  \[Theta]13 = MPTAvoidAmbiguity[\[Theta]13];
	  \[Theta]23 = MPTAvoidAmbiguity[\[Theta]23];
	  \[Delta] = MPTAvoidAmbiguity[\[Delta]];
	  \[Delta]e = MPTAvoidAmbiguity[\[Delta]e];
	  \[Delta]\[Mu] = MPTAvoidAmbiguity[\[Delta]\[Mu]];
	  \[Delta]\[Tau] = MPTAvoidAmbiguity[\[Delta]\[Tau]];
	  \[Phi]1 = 2 * MPTAvoidAmbiguity[(1/2)*\[Phi]1];
	  \[Phi]2 = 2 * MPTAvoidAmbiguity[(1/2)*\[Phi]2];
      tmpDev = Abs[MatrixNorm[N[Abs[MPT3x3UnitaryMatrix[\[Theta]12, 
             \[Theta]13, \[Theta]23, \[Delta], \[Delta]e, \[Delta]\[Mu], 
             \[Delta]\[Tau], \[Phi]1, \[Phi]2] - U]]]]; 
      MPTDebugPrint["MatrixNorm of the deviation : ", tmpDev]; 
      If[tmpDev > 0.1, Print["Warning : Deviation = ", tmpDev]; 
	    Throw["MPT3x3MixingParameters::Inconsistent result"];
        MPTDebugPrint["U=", MatrixForm[U]]; 
		MPTDebugPrint["Parameters reinserted: ", 
		 MatrixForm[MPT3x3UnitaryMatrix[\[Theta]12, \[Theta]13, \[Theta]23, 
           \[Delta], \[Delta]e, \[Delta]\[Mu], \[Delta]\[Tau], \[Phi]1, 
           \[Phi]2]]]; ]; 
		MPTDebugPrint["Accuracy = ", tmpDev, 
       " : Proportionality-constant from \!\(J\_CP\) : ", 
       Sin[2*\[Theta]12]*Sin[2*\[Theta]13]*Sin[2*\[Theta]23]*
        Cos[\[Theta]13]]; 
	  Return[{\[Theta]12, \[Theta]13, \[Theta]23, \[Delta], \[Delta]e, 
	  \[Delta]\[Mu], \[Delta]\[Tau], \[Phi]1, \[Phi]2}]];


MNSParameters[m_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ] := 
	Block[{\[Theta]13, \[Theta]12, \[Theta]23, 
      \[Delta], \[Delta]\[Mu], \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2, 
      tmpUL, tmpU, tmpmdiag},
	  If[MatrixNorm[m - Transpose[m]] > 0.01, 
	   Throw["MNSParameters::Encountered non-symmetric matrix!"]];
	  tmpUL = MPT3x3MixingMatrixL[Ye];
	  tmpU = MPT3x3NeutrinoMixingMatrix[Transpose[tmpUL] . m . tmpUL];	  
	  tmpmdiag = Chop[Transpose[tmpU] . Transpose[tmpUL] . m . tmpUL . tmpU];
	  {\[Theta]13, \[Theta]12, \[Theta]23, \[Delta], \[Delta]\[Mu], 
	  \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2} = 
	  MPT3x3MixingParameters[tmpU];
	  \[Phi]1 = MPTAvoidAmbiguity[\[Phi]1];
	  \[Phi]2 = MPTAvoidAmbiguity[\[Phi]2];
	  Return[{{\[Theta]13, \[Theta]12, \[Theta]23, \[Delta], \[Delta]\[Mu], 
	  \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2},
	  {tmpmdiag[[1,1]],tmpmdiag[[2,2]],tmpmdiag[[3,3]]},
	  Sort[SingularValueList[N[Ye]]]}]];

MNSParameters[m_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ,S_String] := 
	Block[{\[Theta]13, \[Theta]12, \[Theta]23, 
      \[Delta], \[Delta]\[Mu], \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2, 
      tmpUL, tmpU, tmpmdiag},
	  If[MatrixNorm[m - Transpose[m]] > 0.01, 
	   Throw["MNSParameters::Encountered non-symmetric matrix!"]];
	  tmpUL = MPT3x3MixingMatrixL[Ye];
	  tmpU = MPT3x3NeutrinoMixingMatrix[Transpose[tmpUL] . m . tmpUL, S];
	  tmpmdiag = Chop[Transpose[tmpU] . Transpose[tmpUL] . m . tmpUL . tmpU];
	  {\[Theta]13, \[Theta]12, \[Theta]23, \[Delta], \[Delta]\[Mu], 
	  \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2} = 
	  MPT3x3MixingParameters[tmpU];
	  \[Phi]1 = MPTAvoidAmbiguity[\[Phi]1];
	  \[Phi]2 = MPTAvoidAmbiguity[\[Phi]2];
	  Return[{{\[Theta]13, \[Theta]12, \[Theta]23, \[Delta], \[Delta]\[Mu], 
	  \[Delta]\[Tau], \[Delta]e, \[Phi]1, \[Phi]2},
  	  {tmpmdiag[[1,1]],tmpmdiag[[2,2]],tmpmdiag[[3,3]]},
	  Sort[SingularValueList[N[Ye]]]}]];

MNSMatrix[m_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ] := 
	MPT3x3UnitaryMatrix[MNSParameters[m,Ye][[1]]];

MNSMatrix[m_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ,S_String] := 
	MPT3x3UnitaryMatrix[MNSParameters[m,Ye,S][[1]]];

DiracMNSMatrix[Yn_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ] := 
	Block[{tmpUL, tmpYnPrime, tmpU, tmpYYdiag},
  	 tmpUL = Catch[MPT3x3MixingMatrixL[Ye]];
	 tmpYnPrime = Yn . tmpUL;
	 tmpU = Catch[MPT3x3MixingMatrixL[tmpYnPrime]];
	 tmpYYdiag =		
		Abs[Conjugate[Transpose[tmpU]].Conjugate[Transpose[tmpYnPrime]].tmpYnPrime.tmpU];
	 If[tmpYYdiag[[3,3]] - tmpYYdiag[[2,2]] < tmpYYdiag[[2,2]] - tmpYYdiag[[1,1]],
	  tmpU = 
	  Transpose[{Transpose[tmpU][[2]],Transpose[tmpU][[3]],Transpose[tmpU][[1]]}]];
	  (*inverted hierarchy*)
	 If[Abs[tmpU[[1,1]]] < Abs[tmpU[[1,2]]],
	  tmpU = 
	  Transpose[{Transpose[tmpU][[2]],Transpose[tmpU][[1]],Transpose[tmpU][[3]]}]];	 
	  (* fix theta12*)
	 (*Print["U Y Y U=",
		Chop[Conjugate[Transpose[tmpU]].Conjugate[Transpose[tmpYnPrime]].tmpYnPrime.tmpU]];*)
	 Return[tmpU];
	];

DiracMNSParameters[Yn_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ] := 
	Block[{tmpUL, tmpYnPrime, tmpUMNS, tmpList, tmpmDD},
  	 tmpUL = Catch[MPT3x3MixingMatrixL[Ye]];
	 tmpYnPrime = Yn . tmpUL;
	 tmpUMNS = DiracMNSMatrix[Yn,Ye];
	 tmpmDD = 
	 	Conjugate[Transpose[tmpUMNS]].Conjugate[Transpose[tmpYnPrime]].tmpYnPrime.tmpUMNS;
	 tmpList = Catch[MPT3x3MixingParameters[tmpUMNS]];
	 Return[{{tmpList[[1]],tmpList[[2]],tmpList[[3]],tmpList[[4]]},
	  {Sqrt[Abs[tmpmDD[[1,1]]]],Sqrt[Abs[tmpmDD[[2,2]]]],Sqrt[Abs[tmpmDD[[3,3]]]]},
	   Sort[SingularValueList[Ye]]}]
	];

CKMParameters[Yu_?MPTNumericMatrixQ,Yd_?MPTNumericMatrixQ] := 
     Block[{tmpList},
	  tmpList = Catch[MPT3x3MixingParameters[CKMMatrix[Yu,Yd]]];
	  Return[{{tmpList[[1]],tmpList[[2]],tmpList[[3]],tmpList[[4]]},
	   Sort[SingularValueList[Yu]],Sort[SingularValueList[Yd]]}]
	 ];

CKMReplacementRules[Yu_?MPTNumericMatrixQ,Yd_?MPTNumericMatrixQ] :=
	Block[{tmpCKM},
	 tmpCKM = Catch[CKMParameters[Yu,Yd]];
	 Return[{MPT\[Theta]12 -> tmpCKM[[1,1]], MPT\[Theta]13 -> tmpCKM[[1,2]], 
	  MPT\[Theta]23 -> tmpCKM[[1,3]], MPT\[Delta] -> tmpCKM[[1,4]],
	  MPTyu -> tmpCKM[[2,1]], MPTyc -> tmpCKM[[2,2]], MPTyt -> tmpCKM[[2,3]],
	  MPTyd -> tmpCKM[[3,1]], MPTys -> tmpCKM[[3,2]], MPTyb -> tmpCKM[[3,3]]}]
	];
	
CKMMatrix[Yu_?MPTNumericMatrixQ,Yd_?MPTNumericMatrixQ] := 
	Block[{tmpUL, tmpUCKM},
  	  tmpUL = Catch[MPT3x3MixingMatrixL[Yu]];
	  tmpUCKM = Catch[MPT3x3MixingMatrixL[Yd . tmpUL]];
	  Return[tmpUCKM]];

End[]; 
EndPackage[];
