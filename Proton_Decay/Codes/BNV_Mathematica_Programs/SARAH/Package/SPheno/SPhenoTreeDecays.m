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
GenerateSPhenoTreeLevelDecays[Eigenstates_]:=Block[{i},
(*
Print["--------------------------------------"];
Print["Writing Two Body Decays for SPheno "];
Print["--------------------------------------"];
*)
Print[StyleForm["Write tree-level decays","Section",FontSize->12]];
Particles[Current]=Particles[Eigenstates];

MakeCouplingLists;


(* $sarahCurrentSPhenoDir=ToFileName[{$sarahCurrentOutputDir,"SPheno"}]; *)
(* CreateDirectory[$sarahCurrentSPhenoDir]; *)
sphenoDecay=OpenWrite[ToFileName[$sarahCurrentSPhenoDir,"TreeLevel_Decays_"<>ModelName<>".f90"]];

WriteHeaderDecays;


savedDecayInfos={};
savedDecayInfos3Body={};
BR2and3={};
BR2={};
All3BodyWidths ={};

MakeWidthList;

Print["  Writing two-body decay: ",Dynamic[DynamicTBDnr],"/",Length[ListDecayParticles],"(",Dynamic[DynamicTBDpart],")"];
For[i=1,i<=Length[ListDecayParticles],
DynamicTBDnr=i;
DynamicTBDpart=ListDecayParticles[[i]];
MakeDecayLists[ListDecayParticles[[i]]];
i++;];
DynamicTBDpart="All Done";

NeededIntegralsComplete={};

(*
Print["--------------------------------------"];
Print["Writing Three Body Decays for SPheno "];
Print["--------------------------------------"];
*)

Print["  Writing three-body decay: ",Dynamic[Dynamic3BDnr],"/",Length[ListDecayParticles3B],"(",Dynamic[Dynamic3BDpart],")"];

If[Length[ListDecayParticles3B]>0,
MakeListPhaseSpaceIntegrals[25000,500000,9,12,10,12,16];
 For[i=1,i<=Length[ListDecayParticles3B],
Dynamic3BDnr=i;
Dynamic3BDpart=ListDecayParticles3B[[i,1]];
MakeDecayLists3Body[ListDecayParticles3B[[i,1]],ListDecayParticles3B[[i,2]]];
i++;];
Dynamic3BDpart="All Done";

];


AppendSourceCode["ScalarToTwoVectorbosons.f90",sphenoDecay];

WriteString[sphenoDecay, "End Module TreeLevelDecays_"<>ModelName<>" \n \n"];

Close[sphenoDecay];

];

MakeWidthList :=Block[{i,particle,type},
SPhenoWidthBR = {};
SPhenoWidthBR1L = {};

For[i=1,i<=Length[ListDecayParticles],
SPhenoWidthBR = Join[SPhenoWidthBR,{ToExpression["gP"<>ToString[ListDecayParticles[[i]]]],ToExpression["gT"<>ToString[ListDecayParticles[[i]]]],ToExpression["BR"<>ToString[ListDecayParticles[[i]]]]}];
realVar= Join[realVar,{ToExpression["gP"<>ToString[ListDecayParticles[[i]]]],ToExpression["gT"<>ToString[ListDecayParticles[[i]]]],ToExpression["BR"<>ToString[ListDecayParticles[[i]]]]}];

SPhenoWidthBR1L = Join[SPhenoWidthBR1L,{ToExpression["gP"<>"1L"<>ToString[ListDecayParticles[[i]]]],ToExpression["gT"<>"1L"<>ToString[ListDecayParticles[[i]]]],ToExpression["BR"<>"1L"<>ToString[ListDecayParticles[[i]]]]}];
realVar= Join[realVar,{ToExpression["gP"<>"1L"<>ToString[ListDecayParticles[[i]]]],ToExpression["gT"<>"1L"<>ToString[ListDecayParticles[[i]]]],ToExpression["BR"<>"1L"<>ToString[ListDecayParticles[[i]]]]}];


i++;];

For[i=1,i<=Length[Particles[Current]],
particle=Particles[Current][[i,1]];
If[FreeQ[ListDecayParticles,particle]==True &&  FreeQ[ListDecayParticles,particle /. diracSubBack1[SPheno`Eigenstates] /. diracSubBack2[SPheno`Eigenstates]],
type=getType[particle];
If[(type===S || type ===V || type === F) && FreeQ[massless,particle]==True && FreeQ[SPhenoParameters,SPhenoWidth[particle]]==True,
realVar=Join[realVar,{SPhenoWidth[particle]}];
 If[getGenSPheno[particle]>1,
SPhenoParameters=Join[SPhenoParameters,{{SPhenoWidth[particle],{generation},{getGenSPheno[particle]}}}];
SPhenoParameters=Join[SPhenoParameters,{{SPhenoWidth1L[particle],{generation},{getGenSPheno[particle]}}}];,
SPhenoParameters=Join[SPhenoParameters,{{SPhenoWidth[particle],{},{}}}];
SPhenoParameters=Join[SPhenoParameters,{{SPhenoWidth1L[particle],{},{}}}];
]; 
];
];
i++;];

SPhenoWidthVP={};

For[i=1,i<=Length[Particles[Current]],
If[getType[Particles[Current][[i,1]]]==V && FreeQ[massless,Particles[Current][[i,1]]] && FreeQ[ListDecayParticles,Particles[Current][[i,1]]],
SPhenoWidthVP = Join[SPhenoWidthVP,{SPhenoWidth[Particles[Current][[i,1]]]}];
realVar=Join[realVar,{SPhenoWidth[SPhenoWidthVP]}];
];
i++;];

];


WriteHeaderDecays:=Block[{},
WriteCopyRight[sphenoDecay];

WriteString[sphenoDecay, "Module TreeLevelDecays_"<>ModelName<>"\n \n"];
WriteString[sphenoDecay, "Use Control \n"];
WriteString[sphenoDecay, "Use DecayFunctions \n"];
WriteString[sphenoDecay, "Use Settings \n"];
WriteString[sphenoDecay, "Use LoopCouplings_"<>ModelName<>" \n"];
WriteString[sphenoDecay, "Use CouplingsForDecays_"<>ModelName<>" \n"];
WriteString[sphenoDecay, "Use Model_Data_"<>ModelName<>" \n"];
WriteString[sphenoDecay, "Use Mathematics, Only: Li2 \n \n "];
WriteString[sphenoDecay, "Contains \n \n  \n"];
];


(* ::Input::Initialization:: *)

MakeDecayLists[particle_]:=Block[{i,temp},
(* Print["Decay of ", particle]; *)

SA`SkipFields=SkipFields=Select[Transpose[PART[S]][[1]],(getGenSPhenoStart[#]>getGenSPheno[#])&];

If[SA`AddOneLoopDecay===True,
temp = TwoBodyDecayAllAllowed[particle];,
temp = TwoBodyDecay[particle]; 
];
ProcessList={};
For[i=1,i<=Length[temp],
If[Select[SA`SkipFields,(FreeQ[temp[[i]] /. conj[x_]->x,#]==False)&]==={},
ProcessList=Join[ProcessList,{temp[[i]]}];
];
i++];


NeededMasses={SPhenoMass[particle]};



If[particle=!=HiggsBoson, 
If[particle=!=PseudoScalar,
NeededCouplings={};
NeededCouplingsInsert={};
FullInformation = {};,
NeededCouplings={{cplPseudoHiggsPP},{cplPseudoHiggsGG}};
NeededCouplingsInsert={};
FullInformation = {{VectorP,VectorP,{{{cplPseudoHiggsPP},{PseudoScalar,VectorP,VectorP}}},1,1},{VectorG,VectorG,{{{cplPseudoHiggsGG},{PseudoScalar,VectorG,VectorG}}},1,1}};
];,
NeededCouplings={{cplHiggsPP},{cplHiggsGG},{cplHiggsZZvirt},{cplHiggsWWvirt}};
NeededCouplingsInsert={};
FullInformation = {{VectorP,VectorP,{{{cplHiggsPP},{HiggsBoson,VectorP,VectorP}}},1,1},{VectorG,VectorG,{{{cplHiggsGG},{HiggsBoson,VectorG,VectorG}}},1,1},{VectorZ,VectorZ,{{{cplHiggsZZvirt},{HiggsBoson,VectorZ,VectorZ}}},1,1},{VectorW,VectorW,{{{cplHiggsWWvirt},{HiggsBoson,VectorW,VectorW}}},1,1}};
 ];


For[i=1,i<=Length[ProcessList],
If[FreeQ[NeededMasses, SPhenoMass[ProcessList[[i,1]]]]==True && FreeQ[massless,getBlank[ProcessList[[i,1]]]]==True,
NeededMasses=Join[NeededMasses,{SPhenoMass[ProcessList[[i,1]]]}];
];
If[FreeQ[NeededMasses, SPhenoMass[ProcessList[[i,2]]]]==True && FreeQ[massless,getBlank[ProcessList[[i,2]]]]==True,
NeededMasses=Join[NeededMasses,{SPhenoMass[ProcessList[[i,2]]]}];
];

If[FreeQ[NeededCouplingsInsert, ProcessList[[i,3]]]==True && ProcessList[[i,3]]=!=LOOP,
NeededCouplingsInsert=Join[NeededCouplingsInsert,{{ProcessList[[i,3]]}}];
NeededCouplings=Join[NeededCouplings,{getSPhenoCoupling[ProcessList[[i,3]]][[1]]}];
];

If[ProcessList[[i,3]]===LOOP,
FullInformation = Join[FullInformation,{{ProcessList[[i,1]],ProcessList[[i,2]],LOOP,ProcessList[[i,4]],ProcessList[[i,5]]}}];,
FullInformation = Join[FullInformation,{{ProcessList[[i,1]],ProcessList[[i,2]],{getSPhenoCoupling[ProcessList[[i,3]]]},ProcessList[[i,4]],ProcessList[[i,5]]}}];
];

i++;];

If[(particle===HiggsBoson || particle===PseudoScalar) && AddSMrunning=!=False && coupAlphaStrong=!={}, 
NeededMasses = Intersection[Join[NeededMasses,SPhenoMass/@Transpose[coupAlphaStrong][[1]]]];
];

savedDecayInfos = Join[savedDecayInfos,{{particle,NeededCouplings,NeededMasses,FullInformation}}];

channels=CountNumberEntries[FullInformation];
(*
channels=0;
For[i=1,i\[LessEqual]Length[FullInformation],
channels += (1+getGenSPheno[FullInformation[[i,1]]]-getGenSPhenoStart[FullInformation[[i,1]]])*(1+getGenSPheno[FullInformation[[i,2]]]-getGenSPhenoStart[FullInformation[[i,2]]]);
i++;];
*)
SPhenoParameters= Join[SPhenoParameters,{{ToExpression["gP"<>ToString[particle]],{generation,generation},{getGenSPheno[particle],channels}}}];
SPhenoParameters= Join[SPhenoParameters,{{ToExpression["gT"<>ToString[particle]],{generation},{getGenSPheno[particle]}}}];
SPhenoParameters = Join[SPhenoParameters,{{ToExpression["BR"<>ToString[particle]],{generation,generation},{getGenSPheno[particle],channels}}}];

SPhenoParameters= Join[SPhenoParameters,{{ToExpression["gP1L"<>ToString[particle]],{generation,generation},{getGenSPheno[particle],channels}}}];
SPhenoParameters= Join[SPhenoParameters,{{ToExpression["gT1L"<>ToString[particle]],{generation},{getGenSPheno[particle]}}}];
SPhenoParameters = Join[SPhenoParameters,{{ToExpression["BR1L"<>ToString[particle]],{generation,generation},{getGenSPheno[particle],channels}}}];


BR2=Join[BR2,{{particle,channels}}];
WriteTwoBodyDecay[particle,NeededCouplings,NeededMasses,FullInformation];

];




WriteTwoBodyDecay[particle_, couplings_,masses_,processes_]:=Block[{i,n2, factor,temp,ind,t1,t2,pt1,pt2},

(* MakeSubroutineTitle[ToString[particle]<>"TwoBodyDecay",Flatten[{couplings,masses}],{"i_in"},{"gPartial","gT","BR"},sphenoDecay]; 

MakeSubroutineTitle[ToString[particle]<>"TwoBodyDecay",Flatten[{listAllParametersAndVEVs,masses}],{"i_in","deltaM"},{"gPartial","gT","BR"},sphenoDecay]; *)

MakeSubroutineTitle[ToString[particle]<>"TwoBodyDecay",Flatten[{NewMassParameters,listAllParametersAndVEVs}],{"i_in","deltaM"},{"gPartial","gT","BR"},sphenoDecay];

WriteString[sphenoDecay, "Implicit None \n \n"];
(* MakeVariableList[Flatten[{couplings,masses}],",Intent(in)",sphenoDecay]; 
MakeVariableList[Flatten[{listAllParametersAndVEVs,masses}],",Intent(in)",sphenoDecay]; *)
MakeVariableList[Flatten[{listAllParametersAndVEVs,NewMassParameters}],",Intent(in)",sphenoDecay];
MakeVariableList[Flatten[couplings],"",sphenoDecay];
WriteString[sphenoDecay,"Integer, Intent(in) :: i_in \n"];
If[getGenSPheno[particle]>1,
WriteString[sphenoDecay, "Real(dp), Intent(inout) :: gPartial(:,:), gT(:) \n"];,
WriteString[sphenoDecay, "Real(dp), Intent(inout) :: gPartial(:,:), gT \n"];
];
WriteString[sphenoDecay, "Real(dp), Intent(in) :: deltaM \n"];
WriteString[sphenoDecay, "Real(dp), Optional, Intent(inout) :: BR(:,:) \n"];
WriteString[sphenoDecay, "Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 \n"];
WriteString[sphenoDecay, "Real(dp) :: gam, m_in, m1out, m2out, coupReal \n"];
WriteString[sphenoDecay, "Complex(dp) :: coupC, coupR, coupL, coup \n \n"];

If[particle === HiggsBoson || particle === PseudoScalar,
WriteString[sphenoDecay, "Real(dp) :: alpha3 \n"];
];




WriteString[sphenoDecay, "Iname = Iname + 1 \n"];
WriteString[sphenoDecay, "NameOfUnit(Iname) = '"<>ToString[particle]<>"TwoBodyDecay'\n \n"];





genIN = ToString[getGenSPheno[particle]];


If[getGenSPheno[particle]>1,
WriteString[sphenoDecay, "If (i_in.Lt.0) Then \n"];
WriteString[sphenoDecay, "  i_start = "<>ToString[getGenSPhenoStart[particle]]<>" \n"];
WriteString[sphenoDecay, "  i_end = "<>genIN <> "\n"];
WriteString[sphenoDecay, "  gT = 0._dp \n"];
WriteString[sphenoDecay, "  gPartial = 0._dp \n"];
WriteString[sphenoDecay, "Else If ( (i_in.Ge.1).And.(i_in.Le."<>genIN<>") ) Then \n"];
WriteString[sphenoDecay, "  i_start = i_in \n"];
WriteString[sphenoDecay, "  i_end = i_in \n"];
WriteString[sphenoDecay, "  gT(i_in) = 0._dp \n"];
WriteString[sphenoDecay, "  gPartial(i_in,:) = 0._dp \n"];
WriteString[sphenoDecay, "Else \n"];
WriteString[sphenoDecay, "  If (ErrorLevel.Ge.-1) Then \n"];
WriteString[sphenoDecay, "     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) \n"];
WriteString[sphenoDecay,  "     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,"<>genIN<>"\n\n"];
WriteString[sphenoDecay, "     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) \n"];
WriteString[sphenoDecay,  "     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,"<>genIN<>"\n\n"];
WriteString[sphenoDecay, "  End If \n"];
WriteString[sphenoDecay, "  If (ErrorLevel.Gt.0) Call TerminateProgram \n"];
WriteString[sphenoDecay, "  If (Present(BR)) BR = 0._dp \n"];
WriteString[sphenoDecay, "  Iname = Iname -1 \n"];
WriteString[sphenoDecay, "  Return \n"];
WriteString[sphenoDecay, "End If \n \n"];,

WriteString[sphenoDecay, "If (i_in.Lt.0) Then \n"];
WriteString[sphenoDecay, "  i_start = 1 \n"];
WriteString[sphenoDecay, "  i_end = 1 \n"];
WriteString[sphenoDecay, "  gT = 0._dp \n"];
WriteString[sphenoDecay, "  gPartial = 0._dp \n"];
WriteString[sphenoDecay, "Else \n"];
WriteString[sphenoDecay, "  If (ErrorLevel.Ge.-1) Then \n"];
WriteString[sphenoDecay, "     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) \n"];
WriteString[sphenoDecay,  "     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,"<>genIN<>"\n\n"];
WriteString[sphenoDecay, "     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) \n"];
WriteString[sphenoDecay,  "     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,"<>genIN<>"\n\n"];
WriteString[sphenoDecay, "  End If \n"];
WriteString[sphenoDecay, "  If (ErrorLevel.Gt.0) Call TerminateProgram \n"];
WriteString[sphenoDecay, "  If (Present(BR)) BR = 0._dp \n"];
WriteString[sphenoDecay, "  Iname = Iname -1 \n"];
WriteString[sphenoDecay, "  Return \n"];
WriteString[sphenoDecay, "End If \n \n"];



];


 If[getGenSPheno[particle]>1,
WriteString[sphenoDecay, "Do i1=i_start,i_end \n"];
WriteString[sphenoDecay,"m_in = "<> SPhenoMass[particle,i1]<>" \n"];
WriteString[sphenoDecay, "If (m_in.Eq.0._dp) Cycle \n"];,
WriteString[sphenoDecay, "i1=1 \n"];
WriteString[sphenoDecay,"m_in = "<> SPhenoMass[particle,i1]<>" \n"];
]; 

MakeCall["CouplingsFor_"<>SPhenoForm[particle]<>"_decays_"<>"2B",Flatten[{NewMassParameters,listAllParametersAndVEVs,couplings}],{"m_in","i1"},{"deltaM"},sphenoDecay];


If[AddSMrunning=!=False, 
If[particle === HiggsBoson || particle === PseudoScalar,
WriteString[sphenoDecay,"!alpha3 = AlphaSDR(m_in,"];

For[i=1,i<=Length[coupAlphaStrong],
WriteString[sphenoDecay,ToString[SPhenoMass[coupAlphaStrong[[i,1]]]]];
If[i!= Length[coupAlphaStrong],
WriteString[sphenoDecay,","];
];
i++;];
WriteString[sphenoDecay,") \n"];
WriteString[sphenoDecay,"alpha3 = g3running**2/(4._dp*Pi) \n"];
];
]; 


WriteString[sphenoDecay, "i_count = 1 \n"];



For[i=1,i<=Length[processes],
pt1=processes[[i,1]];
pt2=processes[[i,2]];
t1=getType[pt1];
t2=getType[pt2];
Which[
t1=== F && t2=== F,p1=pt1;p2=pt2;,
t1=== S && t2=== S,p1=pt1;p2=pt2;,
t1=== F && t2=== S,p1=pt1;p2=pt2;,
t1=== S && t2=== F,p1=pt2;p2=pt1;,
t1=== F && t2=== V,p1=pt1;p2=pt2;,
t1=== V && t2=== F,p1=pt2;p2=pt1;,
t1=== V && t2=== S,p1=pt2;p2=pt1;,
t1=== S && t2=== V,p1=pt1;p2=pt2;,
t1=== V && t2=== V,p1=pt2;p2=pt1;
];

WriteString[sphenoDecay,"\n \n"];
WriteString[sphenoDecay,"! ----------------------------------------------\n"];
WriteString[sphenoDecay,"! " <> ToString[p1] <>", " <>ToString[p2] <>"\n"];
WriteString[sphenoDecay,"! ----------------------------------------------\n"];
WriteString[sphenoDecay,"\n \n"];


If[getGenSPheno[p1]>1,
WriteString[sphenoDecay,"Do gt1= "<>ToString[getGenSPhenoStart[p1]] <> ", "<>ToString[getGenSPheno[p1]] <> "\n"];
];
If[getGenSPheno[p2]>1,
If[p1===p2,
WriteString[sphenoDecay,"  Do gt2= gt1, "<>ToString[getGenSPheno[p2]] <> "\n"];,WriteString[sphenoDecay,"  Do gt2="<>ToString[getGenSPhenoStart[p2]] <> ", "<>ToString[getGenSPheno[p2]] <> "\n"];
];
];

WriteString[sphenoDecay,"m1out = "<>SPhenoMass[p1,gt1] <>"\n"];
WriteString[sphenoDecay,"m2out = "<>SPhenoMass[p2,gt2] <>"\n"];

If[processes[[i,3]]=!=LOOP,
temp=MakeIndicesCoupling[{AntiField[particle],i1},{p1,gt1},{p2,gt2},processes[[i,3,1,2]]];
ind = temp[[1]];
checkHC=temp[[2]];

 Switch[getType[particle],
S,
Switch[VType[getType[particle], getType[p1],getType[p2]],
FFS,
	WriteString[sphenoDecay, "coupL = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay, "coupR = "<>ToString[processes[[i,3,1,1,2]]]<>ind <> "\n"];
	WriteString[sphenoDecay,"Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) \n"];
		If[AddSMrunning=!=False,
	If[(particle === HiggsBoson || particle === PseudoScalar) &&  (getBlank[p1]===TopQuark || getBlank[p1]===BottomQuark),
	WriteString[sphenoDecay,"gam = gam * FFqcd(m1out,m_in,alpha3) \n"];
	];
	];,
SSS,
	WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]] <>ind <> "\n"];
	WriteString[sphenoDecay,"Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) \n"];,
SSV,
	WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay,"Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) \n"];,
SVV,
	Switch[processes[[i,3,1,1,1]],
	cplHiggsPP,
		If[getGen[particle]>1,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplHiggsPP(i1))**2 \n"];,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplHiggsPP)**2 \n"];
		];,
	cplHiggsGG,
		If[getGen[particle]>1,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplHiggsGG(i1))**2 \n"];,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplHiggsGG)**2 \n"];
		];,
	cplPseudoHiggsPP,
		If[getGen[particle]>1,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplPseudoHiggsPP(i1))**2 \n"];,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplPseudoHiggsPP)**2 \n"];
		];,
	cplPseudoHiggsGG,
		If[getGen[particle]>1,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplPseudoHiggsGG(i1))**2 \n"];,
		WriteString[sphenoDecay,"gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplPseudoHiggsGG)**2 \n"];
		];,
	cplHiggsZZvirt,
		WriteString[sphenoDecay,"If (m_in.le.2._dp*m1out) Then \n"];
		WriteString[sphenoDecay, "coupReal = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
		WriteString[sphenoDecay,"Call ScalarToVectorBosonsVR(m_in,m1out,coupReal,gam) \n"];
		WriteString[sphenoDecay, "Else \n"];
		WriteString[sphenoDecay, "gam = 0._dp \n"];
		WriteString[sphenoDecay, "End if \n"];,
	cplHiggsWWvirt,
		WriteString[sphenoDecay,"If (m_in.le.2._dp*m1out) Then \n"];
		WriteString[sphenoDecay, "coupReal = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
		WriteString[sphenoDecay,"Call ScalarToVectorBosonsVR(m_in,m1out,coupReal,gam) \n"];
		WriteString[sphenoDecay, "Else \n"];
		WriteString[sphenoDecay, "gam = 0._dp \n"];
		WriteString[sphenoDecay, "End if \n"];,
	_,
		WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
		WriteString[sphenoDecay,"Call ScalarToTwoVectorBosons(m_in,m1out,m2out,coup,gam) \n"];
	];
];

If[getGenSPheno[p1]>1,
If[p1===p2,
WriteString[sphenoDecay,"If (gt1.ne.gt2) gam = 2._dp*gam \n"];
];
];

factor = processes[[i,4]]*processes[[i,5]];
If[AntiField[particle]===particle &&  processes[[i,3,1,1,1]]=!= cplHiggsWWvirt,
If[(AntiField[p1]=!=p1 || AntiField[p2]=!=p2) && AntiField[p1]=!=p2,
factor = 2*factor;
];
];,
F,
Switch[VType[getType[particle], getType[p1],getType[p2]],
FFS,
	WriteString[sphenoDecay, "coupL = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay, "coupR = "<>ToString[processes[[i,3,1,1,2]]] <>ind <> "\n"];
	WriteString[sphenoDecay,"Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) \n"];,
FFV,
	WriteString[sphenoDecay, "coupL = "<>ToString[processes[[i,3,1,1,1]]] <>ind <> "\n"];
	WriteString[sphenoDecay, "coupR = "<>ToString[processes[[i,3,1,1,2]]] <>ind <> "\n"];
	If[(FreeQ[massless,p2]==True || getType[p2]===F) && (FreeQ[massless,p1]==True || getType[p1]===F),
	WriteString[sphenoDecay,"Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) \n"];,
	WriteString[sphenoDecay,"Call FermionToFermionVectorBosonMassless(m_in,m1out,m2out,coupL,coupR,gam) \n"];
		];
];
factor = processes[[i,4]]*processes[[i,5]];
If[AntiField[particle]===particle && (AntiField[p1]=!=p1 || AntiField[p2]=!=p2) && AntiField[p1]=!=p2,factor = 2*factor;];,

V,
Switch[VType[getType[particle], getType[p1],getType[p2]],
SSV,
	WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay,"Call VectorBosonToTwoScalars(m_in,m1out,m2out,1,coup,gam) \n"];,
VVV,
	WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay,"Call VectorBosonToTwoVectorBosons(m_in,m1out,m2out,coup,gam) \n"];,
FFV,
	WriteString[sphenoDecay, "coupL = "<>ToString[processes[[i,3,1,1,1]]] <>ind <> "\n"];
	WriteString[sphenoDecay, "coupR = "<>ToString[processes[[i,3,1,1,2]]] <>ind <> "\n"];
	WriteString[sphenoDecay,"Call VectorBosonToTwoFermions(m_in,m1out,m2out,1,coupL,coupR,gam) \n"];,
SVV,
	WriteString[sphenoDecay, "coup = "<>ToString[processes[[i,3,1,1,1]]]<>ind <> "\n"];
	WriteString[sphenoDecay,"Call VectorBosonToScalarAndVectorBoson(m_in,m1out,m2out,coup,gam) \n"];
];


If[getGenSPheno[p1]>1,
If[p1===p2,
WriteString[sphenoDecay,"If (gt1.ne.gt2) gam = 2._dp*gam \n"];
];
];

factor = processes[[i,4]]*processes[[i,5]];
If[AntiField[particle]===particle &&  processes[[i,3,1,1,1]]=!= cplHiggsWWvirt,
If[(AntiField[p1]=!=p1 || AntiField[p2]=!=p2) && AntiField[p1]=!=p2,
factor = 2*factor;
];
];

];

(*
factor = processes[[i,4]]*processes[[i,5]];


If[AntiField[particle]===particle,
If[(AntiField[p1]=!=p1 || AntiField[p2]=!=p2) && AntiField[p1]=!=p2,
factor = 2*factor;
];
];
*)


(* If[Head[factor]===Integer,factor=SPhenoForm[1* factor];,factor=SPhenoForm[factor];]; *)

factor=SPhenoForm[factor];,
factor="0._dp";
WriteString[sphenoDecay,"gam = 0._dp \n"];
];

If[getGenSPheno[particle]>1,
WriteString[sphenoDecay, "gPartial(i1,i_count) = "<> factor <> "*gam \n"];
WriteString[sphenoDecay, "gT(i1) = gT(i1) + gPartial(i1,i_count) \n"];,
WriteString[sphenoDecay, "gPartial(1,i_count) = "<> factor <> "*gam \n"];
WriteString[sphenoDecay, "gT = gT + gPartial(1,i_count) \n"];
];

Switch[particle,
HiggsBoson,
	If[RE[p1]===Lepton && RE[p2]===Lepton,
		WriteString[sphenoDecay,"If ((gt1.le.3).and.(gt2.le.3)) Then \n"];
		WriteString[sphenoDecay, "  BR_Hll(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
		WriteString[sphenoDecay,"End if\n"];
	];
	If[p1===HiggsBoson && p2===HiggsBoson,
	If[getGenSPheno[HiggsBoson]==1,
	WriteString[sphenoDecay, "  BRHHH(1,1) = gPartial(1,i_count) \n"];,
	WriteString[sphenoDecay,"If (gt1.eq.gt2) Then \n"];
	WriteString[sphenoDecay, "  BRHHH(i1,gt1) = gPartial(i1,i_count) \n"];
	WriteString[sphenoDecay, "End if \n"];
	WriteString[sphenoDecay, "  BRHHHijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	];
	If[p1===PseudoScalar && p2===PseudoScalar,
	WriteString[sphenoDecay,"If (gt1.eq.gt2) Then \n"];
	WriteString[sphenoDecay, "  BRHAA(i1,gt1) = gPartial(i1,i_count) \n"];
	WriteString[sphenoDecay, "End if \n"]; 
	WriteString[sphenoDecay, "  BRHAAijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===PseudoScalar && p2===HiggsBoson,
	WriteString[sphenoDecay, "  BRHHAijk(i1,gt2,gt1) = gPartial(i1,i_count) \n"];
	];
	If[p2===PseudoScalar && p1===HiggsBoson,
	WriteString[sphenoDecay, "  BRHHAijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===VectorZ && p2===HiggsBoson,
	WriteString[sphenoDecay, "  BRHHZ(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===VectorZ && p2===PseudoScalar,
	WriteString[sphenoDecay, "  BRHAZ(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p2===VectorZ && p1===HiggsBoson,
	WriteString[sphenoDecay, "  BRHHZ(i1,gt1) = gPartial(i1,i_count) \n"];
	];
	If[p2===VectorZ && p1===PseudoScalar,
	WriteString[sphenoDecay, "  BRHAZ(i1,gt1) = gPartial(i1,i_count) \n"];
	];
	If[RE[p1]===VectorW && RE[p2]===ChargedHiggs,
	WriteString[sphenoDecay, "  BRHHpW(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[RE[p2]===VectorW && RE[p1]===ChargedHiggs,
	WriteString[sphenoDecay, "  BRHHpW(i1,gt2) = gPartial(i1,i_count) \n"];
	];,

	PseudoScalar,
	If[RE[p1]===Lepton && RE[p2]===Lepton,
		WriteString[sphenoDecay,"If ((gt1.le.3).and.(gt2.le.3)) Then \n"];
		WriteString[sphenoDecay, "  BR_All(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
		WriteString[sphenoDecay,"End if\n"];
	];
	If[p1===HiggsBoson && p2===HiggsBoson,
	WriteString[sphenoDecay,"If (gt1.eq.gt2) Then \n"];
	WriteString[sphenoDecay, "  BRAHH(i1,gt1) = gPartial(i1,i_count) \n"];
	WriteString[sphenoDecay, "End if \n"];
	WriteString[sphenoDecay, "  BRAHHijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===PseudoScalar && p2===PseudoScalar,
	WriteString[sphenoDecay,"If (gt1.eq.gt2) Then \n"];
	WriteString[sphenoDecay, "  BRAAA(i1,gt1) = gPartial(i1,i_count) \n"];
	WriteString[sphenoDecay, "End if \n"]; 
	WriteString[sphenoDecay, "  BRAAAijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===HiggsBoson && p2===PseudoScalar,
	WriteString[sphenoDecay, "  BRAHAijk(i1,gt1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p2===HiggsBoson && p1===PseudoScalar,
	WriteString[sphenoDecay, "  BRAHAijk(i1,gt2,gt1) = gPartial(i1,i_count) \n"];
	];
	If[p1===VectorZ && p2===HiggsBoson,
	WriteString[sphenoDecay, "  BRAHZ(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p1===VectorZ && p2===PseudoScalar,
	WriteString[sphenoDecay, "  BRAAZ(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[p2===VectorZ && p1===HiggsBoson,
	WriteString[sphenoDecay, "  BRAHZ(i1,gt1) = gPartial(i1,i_count) \n"];
	];
	If[p2===VectorZ && p1===PseudoScalar,
	WriteString[sphenoDecay, "  BRAAZ(i1,gt1) = gPartial(i1,i_count) \n"];
	];
	If[RE[p1]===VectorW && RE[p2]===ChargedHiggs,
	WriteString[sphenoDecay, "  BRAHpW(i1,gt2) = gPartial(i1,i_count) \n"];
	];
	If[RE[p2]===VectorW && RE[p1]===ChargedHiggs,
	WriteString[sphenoDecay, "  BRAHpW(i1,gt2) = gPartial(i1,i_count) \n"];
	];,

	ChargedHiggs,
	Switch[RE[p1],
		VectorW,
			If[p2===VectorZ,
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_HpWZ(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_HpWZ = gPartial(i1,i_count) \n"];
			];
			];
			If[p2===HiggsBoson,
			WriteString[sphenoDecay, "  BR_HpHW(i1,gt2) = gPartial(i1,i_count) \n"];
			];
			If[p2===PseudoScalar,
			WriteString[sphenoDecay, "  BR_HpAW(i1,gt2) = gPartial(i1,i_count) \n"];
			];,
		VectorZ,
			If[RE[p2]===VectorW,
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_HpWZ(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_HpWZ = gPartial(i1,i_count) \n"];
			];
			];,
		HiggsBoson,
			If[RE[p2]===VectorW,
			WriteString[sphenoDecay, "  BR_HpHW(i1,gt1) = gPartial(i1,i_count) \n"];
			];,
		PseudoScalar,
			If[RE[p2]===VectorW,
			WriteString[sphenoDecay, "  BR_HpAW(i1,gt1) = gPartial(i1,i_count) \n"];
			];,
		Neutrino,
			WriteString[sphenoDecay,"If ((gt1.eq.gt2).and.(gt1.eq.3)) Then \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_Htaunu(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Htaunu = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];,
		Lepton,
			WriteString[sphenoDecay,"If ((gt1.eq.gt2).and.(gt1.eq.3)) Then \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_Htaunu(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Htaunu = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];,
		BottomQuark,
			If[RE[p2]===TopQuark,
			WriteString[sphenoDecay,"If ((gt1.eq.3).and.(gt2.eq.2)) Then \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_Hcb(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Hcb = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay,"If ((gt1.eq.2).and.(gt2.eq.2)) Then \n"];
			WriteString[sphenoDecay, "  BR_Hcs(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Hcs = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay,"If ((gt1.eq.3).and.(gt2.eq.3)) Then \n"];
			WriteString[sphenoDecay, "  BR_HpTB(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_HpTB = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			];,
		TopQuark,
			If[RE[p2]===BottomQuark,
			WriteString[sphenoDecay,"If ((gt1.eq.2).and.(gt2.eq.3)) Then \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_Hcb(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Hcb = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			WriteString[sphenoDecay,"If ((gt1.eq.2).and.(gt2.eq.2)) Then \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay, "  BR_Hcs(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_Hcs = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			If[getGen[chargedHiggs]>1,
			WriteString[sphenoDecay,"If ((gt1.eq.3).and.(gt2.eq.3)) Then \n"];
			WriteString[sphenoDecay, "  BR_HpTB(i1) = gPartial(i1,i_count) \n"];,
			WriteString[sphenoDecay, "  BR_HpTB = gPartial(i1,i_count) \n"];
			];
			WriteString[sphenoDecay, "End if \n"];
			];
		];,
	TopQuark,
	Switch[RE[p1],
	VectorW,
		If[RE[p2]===BottomQuark,
			WriteString[sphenoDecay,"If (gt2.eq.3) Then \n"];
			WriteString[sphenoDecay, "  BR_tWb = gPartial(i1,i_count) \n"];
			WriteString[sphenoDecay, "End if \n"];
				];,
	BottomQuark,
		If[RE[p2]===VectorW,
			WriteString[sphenoDecay,"If (gt1.eq.3) Then \n"];
			WriteString[sphenoDecay, "  BR_tWb = gPartial(i1,i_count) \n"];
			WriteString[sphenoDecay, "End if \n"];,
    	If[RE[p2]===ChargedHiggs,
		WriteString[sphenoDecay,"If (gt1.eq.3) Then \n"];
		WriteString[sphenoDecay, "  BR_tHb = gPartial(i1,i_count) \n"];
		WriteString[sphenoDecay, "End if \n"];
		];
		];,
	ChargedHiggs,
		If[RE[p2]===BottomQuark,
		WriteString[sphenoDecay,"If (gt2.eq.3) Then \n"];
		WriteString[sphenoDecay, "  BR_tHb(i1) = gPartial(i1,i_count) \n"];
		WriteString[sphenoDecay, "End if \n"];
				];
	];
];



WriteString[sphenoDecay, "i_count = i_count +1 \n"];

If[getGenSPheno[p1]>1,
WriteString[sphenoDecay,"  End Do \n"];
];
If[getGenSPheno[p2]>1,
WriteString[sphenoDecay,"End Do \n \n"];
];


i++;];

(*
Switch[particle,
HiggsBoson,
	If[getGenSPheno[HiggsBoson]\[Equal]1,
	WriteString[sphenoDecay, "  BRHHH(1,1) = BRHHH(1,1)/gT \n"];,
	WriteString[sphenoDecay, "  BRHHH(i1,:) = BRHHH(i1,:)/gT(i1) \n"];
	WriteString[sphenoDecay, "  BRHAA(i1,:) = BRHAA(i1,:)/gT(i1) \n"];
	];,
PseudoScalar,
	WriteString[sphenoDecay, "  BRAHH(i1,:) = BRAHH(i1,:)/gT(i1) \n"];
	WriteString[sphenoDecay, "  BRAAA(i1,:) = BRAAA(i1,:)/gT(i1) \n"];
];
*)

If[getGenSPheno[particle]>1,
WriteString[sphenoDecay, "If ((Present(BR)).And.(gT(i1).Eq.0)) Then \n"];
WriteString[sphenoDecay, "  BR(i1,:) = 0._dp \n"];
WriteString[sphenoDecay, "Else If (Present(BR)) Then \n"];
WriteString[sphenoDecay, "  BR(i1,:) = gPartial(i1,:)/gT(i1) \n"];
WriteString[sphenoDecay, "End if \n \n"];
WriteString[sphenoDecay, "End Do \n \n"];,
WriteString[sphenoDecay, "If ((Present(BR)).And.(gT.Eq.0)) Then \n"];
WriteString[sphenoDecay, "  BR(1,:) = 0._dp \n"];
WriteString[sphenoDecay, "Else If (Present(BR)) Then \n"];
WriteString[sphenoDecay, "  BR(1,:) = gPartial(1,:)/gT \n"];
WriteString[sphenoDecay, "End if \n \n"];


]; 



WriteString[sphenoDecay, "Iname = Iname - 1 \n \n"];


If[particle === HiggsBoson || particle === PseudoScalar,
WriteString[sphenoDecay, "Contains \n \n"];

If[particle === HiggsBoson ,
AppendSourceCode["FFqcdScalar.f90",sphenoDecay];,
AppendSourceCode["FFqcdPseudoScalar.f90",sphenoDecay];
];

];



WriteString[sphenoDecay, "End Subroutine "<>ToString[particle]<>"TwoBodyDecay"<>"\n \n \n"];

];

CountNumberEntries[FullInformation_]:=Block[{i,j,k,channels},
channels=0;
For[i=1,i<=Length[FullInformation],
If[FullInformation[[i,1]]===FullInformation[[i,2]],
For[j=getGenSPhenoStart[FullInformation[[i,1]]],j<=getGenSPheno[FullInformation[[i,1]]],
For[k=j,k<=getGenSPheno[FullInformation[[i,1]]],
channels++;
k++;];
j++;];,
channels += (1+getGenSPheno[FullInformation[[i,1]]]-getGenSPhenoStart[FullInformation[[i,1]]])*(1+getGenSPheno[FullInformation[[i,2]]]-getGenSPhenoStart[FullInformation[[i,2]]]);
];
i++;];
Return[channels];
];







(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
 


