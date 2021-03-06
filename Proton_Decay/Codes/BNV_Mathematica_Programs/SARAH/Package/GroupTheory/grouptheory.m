(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* ::Input::Initialization:: *)
epsTensorMatrix[n_]:=Block[{i,j,res,NR},
res=epsTensor@@Table[NR[j],{j,1,n}];
For[i=1,i<=n,
res = Table[res,{NR[i],1,n}];
i++];
Return[res];
];


RM[G_, D_,x_][a__Integer]:=RM[G,-D,x][a] /; D<0;
GINV[G_, D_][a__Integer]:=GINV[G,-D][a] /; D<0;



(* ::Input::Initialization:: *)
GaugeInteractionMatrix[group_,dim_]:=Block[{i,j,pos,res,res2,nr1,x1,x2,stringName,indname},

pos=Position[Gauge,group][[1,1]];
nr1=SA`NrIndices[group,dim];
stringName = ToString[Gauge[[pos,3]]];
res = getGenerator[pos,dim,10,1,2]*(RM[group,dim,x1]@@Table[ToExpression[stringName<>appendIndex[[i]]],{i,1,nr1}] /. subGC[1])(RM[group,dim,x2]@@Table[ToExpression[stringName<>appendIndex[[i]]],{i,1,nr1}] /. subGC[2]); 
resSave=res;
For[i=1,i<=nr1,
res=Hold[Sum[Sum[RES,NAME1],NAME2]] /. NAME1 ->{(ToExpression[stringName<>appendIndex[[i]]]/. subGC[1]),1,Gauge[[pos,2,1]]}/. NAME2 -> {(ToExpression[stringName<>appendIndex[[i]]]/. subGC[2]),1,Gauge[[pos,2,1]]} /. RES->res;
res = ReleaseHold[res];
i++;];

(* If[Gauge[[i,5]]===True, *)
res = res /.( (ToExpression["a"<>StringTake[ToString[Gauge[[pos,3]]],3]<>"10"]):>gen10);
(* ]; *)

res = Sum[res,{gen10,1,Gauge[[pos,2,1]]^2-1}];
res=Table[Table[D[D[res,IR[x1][i]],IR[x2][j]],{i,1,dim}],{j,1,dim}];
Off[Part::"pspec"];
Off[Part::"pkspec1"];
ReleaseHold[Hold[Set[LHS,RHS]] /. LHS -> GINV[group,dim][a__Integer] /. RHS -> (res[[a]])];
On[Part::"pspec"];
On[Part::"pkspec1"];

];

getInvariantMatrix[fields_,coup_]:=Block[{i,j,k,res=1},
For[i=1,i<=Length[Gauge],
If[Gauge[[i,2]]=!=U[1],
Off[Part::"pkspec1"];
res=res*(GenerateInvariantsTensor[Gauge[[i,2]],Gauge[[i,3]],Table[Fields[[Position[ListFields,fields[[j]]][[1,1]],3+i]],{j,1,Length[fields]}]]/.CG[a_,b_]:>InvariantMatrixSusyno[a,b] /. 1[]->1);
On[Part::"pkspec1"];
];
i++;];
Return[res];
];

getRepresentationMatrix[field_,nr_]:=Block[{i,j,ind,pos,res={},indFinal={},k},
pos=Position[ListFields,field][[1,1]];
For[j=1,j<=Length[Gauge],
ind=Select[ListFields[[pos]],(FreeQ[#,Gauge[[j,3]]]==False)&];
If[ind=!={},
ind=ind[[1]]; 
indFinal={};
For[k=1,k<=Length[ind[[1]]],
If[ind[[2,k,1]]===Gauge[[j,3]] || ind[[2,k,1]]===-Gauge[[j,3]],
indFinal = Join[indFinal,{ind[[1,k]]}];
];
k++;];
If[indFinal=!={},
res=Join[res,{RM[Gauge[[j,2]],Fields[[pos,3+j]],nr]@@indFinal}];,
If[Gauge[[j,2]]===U[1],res=Join[res,{IR[nr][1]}];,res=Join[res,{1}];];
];,
If[Gauge[[j,2]]===U[1],res=Join[res,{IR[nr][1]}];,res=Join[res,{1}];];
];
j++;];
Return[res];
];


GetCubicDynkin[reps_,N_]:=Block[{known,unknown,iterator=1,i,j,repsDyn={}},

For[i=1,i<=Length[reps],
If[reps[[i]]=!=1,
If[Head[reps[[i]]]===List,
repsDyn=Join[repsDyn,{reps[[i,2]]}];,
repsDyn=Join[repsDyn,{getDynkinLabels[reps[[i]],N]}];
];
];
i++;];

For[i=1,i<=Length[repsDyn],
SA`Dynkin3[repsDyn[[i]],N]=TriangularAnomalyValue[{SusynoForm[SU[N]]},{repsDyn[[i]]}][[1]];
i++;];

];

SA`Dynkin3[a_Integer,dim_]:=SA`Dynkin3[getDynkinLabels[a,dim],dim];
SA`Dynkin3[{a_Integer,b_List},dim_]:=SA`Dynkin3[b,dim];
SA`Dynkin3[a_Integer,dim_]:=0 /;(a===1 && dim>1)

GetQuadraticDynkin[reps_,N_]:=Block[{i,j},
SA`Dynkin2[N,N]=1;
SA`Dynkin2[-N,N]=1;
SA`Dynkin2[0,N]=0;

For[i=1,i<=Length[reps],
If[reps[[i]]===1,
SA`Dynkin2[reps[[i]],N]=0;,
SA`Dynkin2[reps[[i]],N]=TestDim[Abs[reps[[i]]],N][[5]];
];
i++;];

];

GenerateQuadraticAndCubicDynkins:=Block[{i,j,k,reps,ggroups},
PrintDebug["Generate cubic Dynkins"];
DynamicCheckAnomalies="Generate cubic Dynkins"!;
ggroups=Table[Gauge[[i,2]],{i,1,Length[Gauge]}];
reps=Table[{},{Length[ggroups]}];

For[i=1,i<=Length[Gauge],
pos=Position[ggroups,Gauge[[i,2]]][[1,1]];
For[j=1,j<=Length[Fields],
reps[[pos]]=Join[reps[[pos]],{Fields[[j,3+i]]}];
j++;];
i++;];

(* reps=Intersection/@DeleteCases[reps,1,5]; *)

For[i=1,i<=Length[reps],
If[Head[ggroups[[i]]]=== U || Head[ggroups[[i]]]===SU,
If[ggroups[[i]]=!=U[1],
GetQuadraticDynkin[reps[[i]],ggroups[[i,1]]];
If[ggroups[[i]]=!=SU[2],
GetCubicDynkin[reps[[i]],ggroups[[i,1]]];,
SA`Dynkin3[a__,2]=0;
];,
For[j=1,j<=Length[reps[[i]]],
SA`Dynkin3[reps[[i,j]],1]=reps[[i,j]]^3;
SA`Dynkin2[reps[[i,j]],1]=reps[[i,j]]^2;
j++;];
];
];
i++;];

];


GetMultiplicites:=Block[{i,j,k},
For[i=1,i<=Length[Gauge],
For[j=1,j<=Length[Fields],
mfac=Fields[[j,2]];
For[k=1,k<=Length[Gauge],
If[k=!=i&&Gauge[[k,2]]=!=U[1],
mfac=mfac*Abs[Fields[[j,3+k]]];
];
k++;];
MultiplicityFactorSF[Fields[[j,3]],Gauge[[i,3]]]=mfac;
j++;];
i++;];

];

getDynkinLabels[dim_,group_]:=Block[{},
If[Abs[dim]===1,Return[{0}]];
Switch[Head[group],
SU,Return[getDynkinLabels[dim,group[[1]]]];,
_,Return[getDynkinLabelsSusyno[SusynoForm[group],dim]];
];
];


getDynkinLabels[dim_,N_Integer]:=Block[{},
If[Abs[dim]===1,Return[{0}]];
(* If[dim\[Equal]-2 && N\[Equal]2,Return[{-1}];]; *)
If[N===2 && dim < 0,
Return[-YoungToDynkin[getYoungTableaux[dim,N],N]];,
Return[YoungToDynkin[getYoungTableaux[dim,N],N]];
];
];

DynkinLabels[dim,SU[N]]:=YoungToDynkin[getYoungTableaux[dim,N],N];
YoungToDynkin[tab_,N_]:=If[Length[DeleteCases[tab,0,2]]>N,NOTAB,Table[If[i<Length[tab],tab[[i]]-tab[[i+1]],If[i==Length[tab],tab[[i]],0]],{i,1,N-1}]];
DynkinToYoung[dyn_]:=DeleteCases[Table[Plus@@Take[dyn,{i,Length[dyn]}],{i,1,Length[dyn]}],0,2];
DynkinToDim[dyn_,N_]:=UseHookFormular[DynkinToYoung[dyn],N];

getDynkinLabelsAdjoint[group_]:=Adjoint[SusynoForm[group]];

getNumberStatesAdjoint[group_]:=Switch[Head[group],
U,1,
SU, group[[1]]^2-1,
_,DimR[SusynoForm[group],Adjoint[SusynoForm[group]]]
];

getDimAdjoint[group_]:=getNumberStatesAdjoint[group];

getDimFundamental[group_]:=Switch[Head[group],
U, 1,
SU, group[[1]],
SO, group[[1]],
F, Switch[group[[1]],4,26],
E,Switch[group[[1]],6,27,7,56,8,248]
];

SA`DynL[conj[a_],b_Integer]:=ConjugatedRep[SA`DynL[a,b],Gauge[[b,2]]];
SA`DynL[bar[a_],b_Integer]:=ConjugatedRep[SA`DynL[a,b],Gauge[[b,2]]];
SA`DynL[conj[a_],b_Symbol]:=ConjugatedRep[SA`DynL[a,b],Gauge[[Position[Gauge,b][[1,1]]]][[2]]]/;(FreeQ[Gauge,b]==False);
SA`DynL[bar[a_],b_Symbol]:=ConjugatedRep[SA`DynL[a,b],Gauge[[Position[Gauge,b][[1,1]]]][[2]]]/;(FreeQ[Gauge,b]==False);
SA`DynL[conj[a_],b_Symbol]:=ConjugatedRep[SA`DynL[a,b],AuxGauge[[Position[AuxGauge,b][[1,1]]]][[2]]]/;(FreeQ[Gauge,b]==True && FreeQ[AuxGauge,b]==False );
SA`DynL[bar[a_],b_Symbol]:=ConjugatedRep[SA`DynL[a,b],AuxGauge[[Position[AuxGauge,b][[1,1]]]][[2]]]/;(FreeQ[Gauge,b]==True &&  FreeQ[AuxGauge,b]==False );


SA`DynL[a_,b_]:=SA`DynL[(a/.diracSub[ALL])[[1]],b]/;(FreeQ[diracFermions[ALL],a]==False);

CG[a_,{}][b___]=1;



ConjugatedRepQ[dyn_,groups_]:=Block[{dim},
dim=DimR[SusynoForm[groups],dyn];
If[dyn==={0,2} && groups===SU[3], Return[True];];
If[dyn==={2,0} && groups===SU[3], Return[False];];
If[FreeQ[getDynkinLabelsSusyno[SusynoForm[groups],dim],dyn],
Return[True];,
Return[False];
];
];

ConjugatedRep[dyn_,groups_]:=Block[{i,l,r},
If[dyn==={0},Return[{0}];];
Switch[Head[groups],
U,Return[-dyn];,
SU,
	l=groups[[1]]-1;
	Return[Table[dyn[[l+1-i]],{i,1,Length[dyn]}]];,
_,Return[ConjugateIrrep[SusynoForm[groups],dyn]];
];
];

getTensorIndizes[dim_,group_]:=Block[{temp},
Switch[Head[group],
SU,
	temp=TestDim[dim,group[[1]]];
	Return[temp[[2]]+temp[[3]]];
	
];
];

InitGaugeGroups:=Block[{i,j,k,l,list,reps,groups,pos,rep,crep},
PrintDebug["   Calculate Lie Group constants"];
DynamicInitGaugeG="Construct gauge group constants";
For[i=1,i<=Length[Gauge],
If[Gauge[[i,2]]=!=U[1],
CG[Gauge[[i,2]],{{0},{0}}][a_,b_]:=1;
SA`KnonwCG=Join[SA`KnonwCG,{CG[Gauge[[i,2]],{{0},{0}}]}];
reps=Join[{getDynkinLabelsAdjoint[Gauge[[i,2]]]},DeleteCases[getDynkinLabels[#,Gauge[[i,2]]]&/@Intersection[Transpose[Fields][[3+i]]],{0}]];
(* reps=Abs[reps]; *)
GenerateDynkinCasimir[Gauge[[i,2]] ,#]&/@reps;
GenerateGeneratorsUnbrokenGroup[i,#]&/@reps;
If[Gauge[[i,5]]==True,
GenerateGeneratorsBrokenGroup[i,#]&/@reps; ,
GenerateGeneratorsUnbrokenGroup[i,#]&/@reps;
];
If[Head[Gauge[[i,2]]]=!=SU || Gauge[[i,2,1]]>3,
GeneratorMatrices[Gauge[[i,2]]]=Normal[RepMatrices[SusynoForm[Gauge[[i,2]]],getDynkinLabels[getDimFundamental[Gauge[[i,2]]],Gauge[[i,2]]]]];
];
];
(* Initialize Deltas for N^* N *)
For[j=1,j<=Length[reps],
rep=reps[[j]];
crep=ConjugatedRep[rep,Gauge[[i,2]]];
If[rep=!=crep && crep=!={},
If[FreeQ[SA`KnonwCG,CG[Gauge[[i,2]],{crep,rep}]],
SA`KnonwCG = Join[SA`KnonwCG,{CG[Gauge[[i,2]],{crep,rep}]}];
CG[Gauge[[i,2]],{crep,rep}][a__]=Delta[a];
];
If[FreeQ[SA`KnonwCG,CG[Gauge[[i,2]],{rep,crep}]],
SA`KnonwCG = Join[SA`KnonwCG,{CG[Gauge[[i,2]],{rep,crep}]}];
CG[Gauge[[i,2]],{rep,crep}][a__]=Delta[a];
];
];
j++;];
i++;];

For[i=1,i<=Length[AuxGauge],
If[AuxGauge[[i,2]]=!=U[1],
CG[AuxGauge[[i,2]],{{0},{0}}][a_,b_]:=1;
SA`KnonwCG=Join[SA`KnonwCG,{CG[AuxGauge[[i,2]],{{0},{0}}]}];
reps=Intersection[Join[{getDynkinLabelsAdjoint[AuxGauge[[i,2]]]},DeleteCases[getDynkinLabels[#,AuxGauge[[i,2]]]&/@Intersection[Select[Flatten[Select[Flatten[AuxDimFields,1],FreeQ[#,{AuxGauge[[i,3]],_}]==False&]],Head[#]===Integer&]],{0}]]];
(* reps=Abs[reps]; *)
GenerateDynkinCasimir[AuxGauge[[i,2]] ,#]&/@reps;
GenerateGeneratorsUnbrokenGroup[i,#,True]&/@reps;
If[Head[AuxGauge[[i,2]]]=!=SU || AuxGauge[[i,2,1]]>3,
GeneratorMatrices[AuxGauge[[i,2]]]=Normal[RepMatrices[SusynoForm[AuxGauge[[i,2]]],getDynkinLabels[getDimFundamental[AuxGauge[[i,2]]],AuxGauge[[i,2]]]]];
];
];
i++;];


If[AuxGaugesPresent===True,
For[i=1,i<=Length[AuxGauge],
If[AuxGauge[[i,2]]=!=U[1],
CG[AuxGauge[[i,2]],{{0},{0}}][a_,b_]:=1;
SA`KnonwCG=Join[SA`KnonwCG,{CG[AuxGauge[[i,2]],{{0},{0}}]}];
reps=Join[{getDynkinLabelsAdjoint[AuxGauge[[i,2]]]},DeleteCases[getDynkinLabels[#,AuxGauge[[i,2]]]&/@Intersection[Transpose[Transpose[AuxDimFields][[2]]][[2]]],{0}]];
(* reps=Abs[reps]; *)
GenerateDynkinCasimir[AuxGauge[[i,2]] ,#]&/@reps;
];
(* Initialize Deltas for N^* N *)
For[j=1,j<=Length[reps],
rep=reps[[j]];
crep=ConjugatedRep[rep,AuxGauge[[i,2]]];
If[rep=!=crep && crep=!={},
If[FreeQ[SA`KnonwCG,CG[AuxGauge[[i,2]],{crep,rep}]],
SA`KnonwCG = Join[SA`KnonwCG,{CG[AuxGauge[[i,2]],{crep,rep}]}];
CG[AuxGauge[[i,2]],{crep,rep}][a__]=Delta[a];
];
If[FreeQ[SA`KnonwCG,CG[AuxGauge[[i,2]],{rep,crep}]],
SA`KnonwCG = Join[SA`KnonwCG,{CG[AuxGauge[[i,2]],{rep,crep}]}];
CG[AuxGauge[[i,2]],{rep,crep}][a__]=Delta[a];
];
];
j++;];
i++;];
];

InitStandardSU2;
InitStandardSU3;

conj[CG[a__]]:=CG[a];
Generator[_,{0}][___]=0;
Generator[_,{0},_][___]=0;

];

(* CG[A_,{a_,b_,a_,b_}][i1_,i2_,i3_,i4_]:=CG[A,{a,b}][i1,i2] CG[A,{a,b}][i3,i4]; *)

InitStandardSU2:=Block[{},
CG[SU[2],{{1},{1}}][a__]=epsTensor[a];
CG[SU[2],{{-1},{-1}}][a__]=epsTensor[a];
CG[SU[2],{{-1},{1}}][a__]=Delta[a];
CG[SU[2],{{1},{-1}}][a__]=Delta[a];
 CG[SU[2],{{1},{1},{1},{1}}][a_,b_,c_,d_]:=epsTensor[a,b] epsTensor[c,d]; 

(*
CG[SU[2],{{1},{2},{1}}][a_Integer,b_Integer,c_Integer]:={{{0,0},{0,1/Sqrt[2]},{-1,0}},{{0,-1},{1/Sqrt[2],0},{0,0}}}[[a,b,c]];
CG[SU[2],{{-1},{2},{1}}][a_Integer,b_Integer,c_Integer]:={{{0,-1},{1/Sqrt[2],0},{0,0}},{{0,0},{0,-(1/Sqrt[2])},{1,0}}}[[a,b,c]];
SA`KnonwCG=Join[SA`KnonwCG,{CG[SU[2],{{1},{2},{1}}],CG[SU[2],{{-1},{2},{1}}]}];
*)

(*
CG[SU[2],{{s1_ },{s2_ },{s3_ },{s4_ }}][a_,b_,c_,d_]:=CG[SU[2],{{s1},{s2}}][a,b]CG[SU[2],{{s3},{s4}}][c,d] /; (Abs[s1]\[Equal]1 && Abs[s2]\[Equal] 1&& Abs[s3]\[Equal]1 && Abs[s4]\[Equal]1);

CG[SU[2],{{2},{2},{2},{2}}][a_,b_,c_,d_]:=InvMat[200][a,b] InvMat[200][c,d];
CG[SU[2],{{1},{1},{2},{2}}][a_,b_,c_,d_]:=epsTensor[a,b] InvMat[200][c,d];
CG[SU[2],{{1},{-1},{2},{2}}][a_,b_,c_,d_]:=Delta[a,b]InvMat[200][c,d];
CG[SU[2],{{-1},{1},{2},{2}}][a_,b_,c_,d_]:=Delta[a,b] InvMat[200][c,d];
*)

];

InitStandardSU3:=Block[{},

CG[SU[3],{{0,1},{0,1},{0,1}}][a_,b_,c_]:=epsTensor[a,b,c];
CG[SU[3],{{1,0},{1,0},{1,0}}][a_,b_,c_]:=epsTensor[a,b,c];

CG[SU[3],{{1,0},{0,1}}][a__]=Delta[a];
CG[SU[3],{{0,1},{1,0}}][a__]=Delta[a];
CG[SU[3],{{1,0},{1,0},{1,0}}]=epsTensor;

(*
CG[SU[3],{{0,1},{1,0},{0,1},{1,0}}][a_,b_,c_,d_]=Delta[a,b]Delta[c,d];
CG[SU[3],{{0,1},{1,0},{1,0},{0,1}}][a_,b_,c_,d_]=Delta[a,b]Delta[c,d];
CG[SU[3],{{1,0},{0,1},{1,0},{0,1}}][a_,b_,c_,d_]=Delta[a,b]Delta[c,d];
CG[SU[3],{{1,0},{0,1},{0,1},{1,0}}][a_,b_,c_,d_]=Delta[a,b]Delta[c,d];
*)

];

GenerateDynkinCasimir[group_,dyn_]:=Block[{i,j,k,casimir,dim,dimAdjoint},
If[dyn==={-1} && group===SU[2],
casimir=Casimir[SusynoForm[group],-dyn];
dim=DimR[SusynoForm[group],-dyn];
dimAdjoint=getDimAdjoint[group];,
casimir=Casimir[SusynoForm[group],dyn];
dim=DimR[SusynoForm[group],dyn];
dimAdjoint=getDimAdjoint[group];
];
SA`Casimir[dyn,group]=casimir;
SA`Dynkin[dyn,group]=casimir*dim/dimAdjoint;

(* Generator *)
];

MakeGenerator[nr_,dyn_,cov_,con_]:=Block[{temp,temp2,i,j,complete,name},
name=Gauge[[nr,3]];
complete=cov+con;
(* dimGauge=Gauge[[gaugeNr,2]]; *)
temp=Plus@@Table[(Hold[TA[DIMGAUGE,genf[lor],IndexName[GAUGE,NR]/.subGC[p1],IndexName[GAUGE,NR]/.subGC[p2]]]/.NR->i/.GAUGE->name/.DIMGAUGE->Gauge[[nr,2]])Product[If[j==i,1,(Hold[Delta[IndexName[GAUGE,NR]/.subGC[p1],IndexName[GAUGE,NR]/.subGC[p2]]]/.NR->j/.GAUGE->name)],{j,1,complete}],{i,1,cov}]-Plus@@Table[(Hold[TA[DIMGAUGE,genf[lor],IndexName[GAUGE,NR]/.subGC[p2],IndexName[GAUGE,NR]/.subGC[p1]]]/.NR->i/.GAUGE->name/.DIMGAUGE->Gauge[[nr,2]])Product[If[j==i,1,(Hold[Delta[IndexName[GAUGE,NR]/.subGC[p1],IndexName[GAUGE,NR]/.subGC[p2]]]/.NR->j/.GAUGE->name)],{j,1,complete}],{i,cov+1,complete}];
Generator[Gauge[[nr,2]],dyn,Gauge[[nr,3]]][lor_,p1_,p2_]=temp;
];


GenerateGeneratorsBrokenGroup[nr_,dyn_]:=Block[{i,j,group,ind1, ind2,dim},
group=Gauge[[nr,2]];
If[Head[group]=!=SU, Print["Only broken SU(N) groups are supported!"]; Return[]];
dim=DimR[SusynoForm[group],Abs[dyn]];
If[ConjugatedRepQ[Abs[dyn],group],dim=-dim;];
If[dim===2 && dyn===-{1},dim=-dim;];
Switch[dim,
group[[1]]^2-1, MakeGenerator[nr,dyn,1,1];,
group[[1]],MakeGenerator[nr,dyn,1,0];,
-group[[1]],MakeGenerator[nr,dyn,0,1];,
_,
ytab=TestDim[{DimR[SusynoForm[group],dyn],dyn},group[[1]]];
MakeGenerator[nr,dyn,ytab[[2]],ytab[[3]]];
];
];

GenerateGeneratorsUnbrokenGroup[nr_,dyn_]:=GenerateGeneratorsUnbrokenGroup[nr,dyn,False];
GenerateGeneratorsUnbrokenGroup[nr_,dyn_,auxgauge_]:=Block[{i,j,group},
If[auxgauge=!=True,
group=Gauge[[nr,2]];,
group=AuxGauge[[nr,2]];
];
If[group==SU[2] && dyn=={1},
Generator[SU[2],{1}][a__]=Sig[a]/2;
Return[];
];
If[group==SU[2] && dyn=={-1},
Generator[SU[2],{-1}][a_,b_,c_]=-Sig[a,c,b]/2;
Return[];
];
If[group==SU[2] && dyn=={2},
(* CG[group,{dyn,dyn}][a_,b_]:= Delta[a,b];   *)
Delta3[a_Integer,b_Integer]:={{0,0,1},{0,-1,0},{1,0,0}}[[a,b]];
(* Generator[SU[2],{2}][a__Integer]:=Normal[RepMatrices[SusynoForm[SU[2]],{2}]][[a]];  *)
Generator[SU[2],{2}][a__Integer]:={{{0,-(1/Sqrt[2]),1/Sqrt[2]},{-(1/Sqrt[2]),0,0},{1/Sqrt[2],0,0}},{{0,-(I/Sqrt[2]),-(I/Sqrt[2])},{I/Sqrt[2],0,0},{I/Sqrt[2],0,0}},{{0,0,0},{0,1,0},{0,0,-1}}}[[a]];
Return[];
];
If[group==SU[3] && dyn=={1,0},
Generator[SU[3],{1,0}][a__]=Lam[a]/2;
Return[];
];
If[group==SU[3] && dyn=={0,1},
Generator[SU[3],{0,1}][a_,b_,c_]=-Lam[a,c,b]/2;
Return[];
];
If[group==SU[3] && dyn=={1,1},
CG[group,{dyn,dyn}][a_,b_]:=Delta[a,b];
(* Generator[SU[3],{1,1}][a__Integer]:=Normal[RepMatrices[SusynoForm[SU[3]],{1,1}]][[a]]; *)
Generator[SU[3],{1,1}][a__]=fSU3[a]/I; (* Check!*)
Return[];
];
If[group==SU[3] || group==SU[2],NonFundamentalSU2SU3=True;]; (* to make sure to use the same matrices as Susyno *)

If[Head[group]===SU,
If[ConjugatedRepQ[dyn,group],
Generator[group,dyn][a_,b_,c_]=-Generator[group,ConjugatedRep[dyn,group]][a,c,b];
Return[];,
Off[Part::"pspec"];
Off[Part::"pkspec1"];
Generator[group,dyn][a__Integer]=Normal[RepMatrices[SusynoForm[group],dyn]][[a]]; 
On[Part::"pspec"];
On[Part::"pkspec1"];
Return[];
];
];

repm=Normal[RepMatrices[SusynoForm[group],dyn]];

If[DimR[SusynoForm[group],dyn]===getDimAdjoin[group],
 CG[group,{dyn,dyn}][a_,b_]:=Delta[a,b];
];

Off[Part::"pspec"];
Off[Part::"pkspec1"];
temp=Hold[SetDelayed[Generator[group,dyn][a___Integer],REP[[a]]]]/. REP->repm;
On[Part::"pspec"];
On[Part::"pkspec1"];
ReleaseHold[temp];

SA`SavedGenerators=Join[SA`SavedGenerators,{{Generator[group,dyn],repm}}];


];



SA`DynL[a_,b_Integer]:=SA`DynL[a,Gauge[[b,3]]];
SA`Casimir[a_,b_Integer]:=SA`Casimir[a,Gauge[[b,3]]];
SA`Dynkin[a_,b_Integer]:=SA`Dynkin[a,Gauge[[b,3]]];
SA`MulFactor[a_,b_Integer]:=SA`MulFactor[a,Gauge[[b,3]]];
SA`DimensionGG[a_,b_Integer]:=SA`DimensionGG[a,Gauge[[b,3]]];
Generator[a_,b_Integer]:=Generator[a,Gauge[[b,3]]];
SA`Casimir[conj[a_],b_Integer]:=SA`Casimir[a,Gauge[[b,3]]];
SA`Dynkin[conj[a_],b_Integer]:=SA`Dynkin[a,Gauge[[b,3]]];
SA`MulFactor[conj[a_],b_Integer]:=SA`MulFactor[a,Gauge[[b,3]]];
Generator[conj[a_],b_Integer]:=Generator[a,Gauge[[b,3]]];

CleanUpGaugeConstants:=Block[{i,j,k},
SA`ClebschGordon=Intersection[SA`ClebschGordon];

];


