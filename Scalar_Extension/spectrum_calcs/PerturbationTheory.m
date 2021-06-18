(************************************************************************
 ***********   RAYLEIGH-SCHROEDINGER PERTURBATION THEORY      ***********
 ***********                 version 1.01                     ***********
 ***********  July -93                           D.Nordfors   ***********
 ***********  e-mail: David@Nordfors.se                       ***********
 ************************************************************************
 **  This version manages up to second order perturbaton theory        **
 ************************************************************************
 *** This note may not be removed                                     ***
 ************************************************************************)

BeginPackage["Physics`PerturbationTheory`"]

(*     symbols
===============================================
	ham0  - unperturbed matrix hamiltonian - a diagonal matrix
	ham1  - perturbation hamiltonian
	state - eigenstate number
*)

PTEigenvalue::usage =
	"PTEigenvalue[NthState,H0,H1,NthOrder] calculates the NthState
eigenvalue of the matrix H0+H1 using NthOrder Rayleigh-Schroedinger
perturbation theory. H0 is the unperturbed hamiltonian - a diagonal matrix.
H1 is the perturbation - a matrix of the same dimension as H0. The package
manages up to NthOrder = 2."

PTEigenvector::usage =
	"PTEigenvector[NthState,H0,H1,NthOrder] calculates the NthState
eigenvector of the matrix H0+H1 using NthOrder Rayleigh-Schroedinger
perturbation theory. H0 is the unperturbed hamiltonian - a diagonal matrix.
H1 is the perturbation - a matrix of the same dimension as H0. The package
manages up to NthOrder = 2."


Begin["`Private`"]


(* 
=========================================
    **** EIGENVALUES  *****
=========================================*)

PTEigenvalue[state_,ham0_,ham1_,order_] :=
         Sum[pten[i][state,ham0,ham1],{i,0,order}]


(*  Unperturbed eigenvalues (0th order)   
-------------------*)
pten[0][state_,ham0_,ham1_] := ham0[[state,state]]

(* 1st order correction to eigenvalues
----------------------------------------------------------------------------*)
pten[1][state_,ham0_,ham1_] := ham1[[state,state]]

(* 2nd order correction to eigenvectors
----------------------------------------------------------------------------*)
pten[2][state_,ham0_,ham1_] := 
        Sum[If[state == j, 0,
                ham1[[state,j]]^2/(ham0[[state,state]]-ham0[[j,j]])
              ],
            {j,1,Dimensions[ham1][[1]]}]


(* 
=========================================
    **** EIGENVECTORS *****
=========================================*)

PTEigenvector[state_,ham0_,ham1_,order_] := 
	Sum[ptwf[i][state,ham0,ham1],{i,0,order}]


(*  Unperturbed eigenvectors (0th order) 
----------------------------------------------------------------------------*)
ptwf[0][state_,ham0_,ham1_] := IdentityMatrix[Dimensions[ham0][[1]]][[state]]

(* 1st order correction to eigenvectors
----------------------------------------------------------------------------*)
ptwf[1][state_,ham0_,ham1_] := Block[{e,h},
	h[i_,j_] := ham1[[i,j]];
	e[i_]    := ham0[[i,i]];
	Table[If[state == k, 0,
                 h[state,k]/(e[state]-e[k])
              ],{k,1,Dimensions[ham1][[1]]}]
]

ptwf[2][state_,ham0_,ham1_] := Block[{e,h,dim},
	dim       = Dimensions[ham1][[1]];
	h[i_,j_] := ham1[[i,j]];
	e[i_]    := ham0[[i,i]];
	Table[If[ state == j, 0,
                  Sum[If[k == state, 0,
                        (h[j,k]*h[k,state])/((e[state]-e[k])(e[state]-e[j])) -
                        (h[state,state]*h[j,state])/(e[state]-e[j])^2 
                        ],{k,1,dim}]],
		{j,1,dim} ]
]

End[]
EndPackage[]		

