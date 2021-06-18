(* ::Package:: *)

(* ***
 RunDec: a Mathematica package for running and decoupling of the
 strong coupling and quark masses

 First version: K.G. Chetyrkin, J.H. K"uhn and M. Steinhauser (Jan. 2000)

 v2.0: F. Herren and M. Steinhauser (Jan. 2016)
 v2.1: F. Herren and M. Steinhauser (Jun. 2016)
        
 Improvements and extensions:
 - several minor improvements and bug fixes
 - improved version of mMS2mSI[]: mMS2mSInew[]
 - AsmMSrunexact[] implemented
 - five-loop result for \gamma_m implemented (arXiv:1402.6611)
   preparations for implementation of five-loop beta function 
 - four-loop decoupling terms of \zeta_\alpha_s and \zeta_m implemented
   (only for on-shell heavy quark masses), 
   (hep-ph/0512058, hep-ph/0512060,  arXiv:1502.04719)
 - mMS2OS, mOS2MS, mOS2SI extedned to four loops (arXiv:1502.01030)
 - relations between PS, 1S, RS, RS' and OS, MS, SI masses implemented to N3LO
   (arXiv:1502.01030)

References:

  K.G. Chetyrkin, J.H. Kuhn and M. Steinhauser,
  "RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses"
  Comput. Phys. Commun.  133 (2000) 43,
  arXiv:hep-ph/0004189

  "CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses"
  B. Schmidt, M. Steinhauser
  Comput.Phys.Commun. 183 (2012) 1845-1848,
  arXiv:1201.6149, SFB/CPP-12-03, TTP12-02

*** *)

(* ************************************************************ *)
(* ************************************************************ *)

(*
 The examples shown in the paper can be found in the the file
 RunDec_ex.m
*)

(* ************************************************************ *)

Print["RunDec: a Mathematica package for running and decoupling of the"];
Print["        strong coupling and quark masses"];
Print["by K.G. Chetyrkin, J.H. K\\\"uhn and M. Steinhauser (January 2000)"];
Print["by F. Herren and M. Steinhauser (April 2016, v2.1)"];

(* ************************************************************ *)
(* ************************************************************ *)

(*
 default values of the numerical constants
*)

NumDef = {
    asMz -> 0.118,
    Mz -> 91.18,
    Mt -> 175,
    Mb -> 4.7,
    Mc -> 1.6,
    muc -> 1.2,
    mub -> 3.97,
    Mtau -> 1.777
    };

(* ************************************************************ *)
(* ************************************************************ *)

BeginPackage["RunDec`"];

 fdelm::usage = "Factor in front of 4-loop coefficient; introduced in 
mMS2mOS[], mOS2mMS[] and mOS2mSI[]. This option is also available
in the relations between the MS or SI and thresold masses."

 LamExpl::usage =
    "LamExpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the explicite formulae \\ln(\\mu/\\Lambda) = f(als,mu,nf).";

 LamImpl::usage =
    "LamImpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the implicite formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasLam::usage = 
    "AlphasLam[lam ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) to l-loop accuracy
using the formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasExact::usage = 
    "AlphasExact[als0 ,mu0 ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) integrating
the \beta function numerically. The initial condition \\alpha_s(mu0)=als0
is used.";

 AsmMSrunexactOLD::usage = 
    "AsmMSrunexact[mmu0, asmu0, mu0, mu, nf, l] computes \\alpha_s^(nf)(mu) 
and m_q^(nf)(mu) solving the coupled system of differential equations 
simultaneously. The initial conditions are \\alpha_s(mu0)=als0 
and m_q(mu)=mmu0.";

 mOS2mMS::usage =
    "mOS2mMS[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. {mq} represents a set of light quark masses in the on-shell 
scheme leading to corrections of O(as^2 mq/mOS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mOS2mMS[mOS, nf, l] computes the MS-bar mass at the the scale mOS. 
\\alpha_s^(nf)(M) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mMS2mOS::usage =
    "mMS2mOS[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
{mq} represents a set of light quark masses in the on-shell scheme
leading to corrections of O(as^2 mq/mMS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mMS2mOS[mum, nf, l] computes the on-shell mass where mum is 
defined through mum=m^(nf)(mum).
\\alpha_s^(nf)(mum) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mOS2mMSrun::usage =
    "mOS2mMSrun[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu.
In contrast to mOS2mMS[..] in a first step mMS(mMS) is computed
and only then the conversion to mOS is performed.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mOSrun::usage =
    "mMS2mOSrun[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
In contrast to mMS2mOS[..] in a first step mMS(mMS) is computed
and afterwards mMS(mu) is evaluated.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mOS2mMSit::usage =
    "mOS2mMSit[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. However, in contrast to mOS2mMS[..], the relation
M_OS = m_MS * (1 + ...) is solved iteratively. This has the advantage
that the light quark masses can be evaluated in the MS-bar scheme.";

 mOS2mSI::usage =
    "mOS2mSI[mOS, {mq}, asM, nf, l] computes the mass
\\mu_m = m_MS(\\mu_m). {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS).
asM = \\alpha_s^(nf)(M_OS), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mMS::usage =
    "mMS2mMS[mmu0, asmu0, asmu, nf, l] computes m(mu) from the knowledge
of m(mu0). asmu = \\alpha_s(\\mu), asmu0 = \\alpha_s(\\mu0), nf is the number 
of active flavours and l represents the number of loops.";

 mMS2mSI::usage =
    "mMS2mSI[mMS, asmu, mu, nf, l] computes the scale
invarant mass mMS(mMS). mMS is the MS-bar mass at the scale mu,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mSInew::usage = 
    "mMS2mSInew[mMS, asmu, mu, nf, l] compues mMS(mMS)
without using AlphasLam[] but mMS2mMS[] iteratively";

 mMS2mRI::usage =
    "mMS2mRI[mMS, asmu, nf, l] computes the regularization invariant
mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mRI2mMS::usage =
    "mRI2mMS[mRI, asmu, nf, l] computes the MS-bar quark mass.
mRI is the regularization invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGI::usage =
    "mMS2mRGI[mMS, asmu, nf, l] computes the renormalization group
invariant mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGImod::usage =
    "mMS2mRGImod[mMS, asmu, nf, l] computes the renormalization group
invariant mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
(slightly modified version of mMS2mRGI[], factor '2 beta_0')";

 mRGI2mMS::usage =
    "mRGI2mMS[mRGI, asmu, nf, l] computes the MS-bar quark mass.
mRGI is the renormalization group invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mOS2mPS::usage =
    "mOS2mPS[mOS, mq, asmu, mu, muf, nl, loops] computes the PS mass mPS(muf) at 
the scale mu. {mq} represents a set of light quark masses in the on-shell 
scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nf)(\\mu), nf=nl+1 is the number
of active flavours and loops represents the number of loops.";

 mMS2mPS::usage =
    "mMS2mPS[mMS, mq, asmu, mu, muf, nl, loops, fdelm] 
computes the PS mass mPS(muf) at 
the scale mu. mMS is the MSbar mass at the scale mu. 
asmu = \\alpha_s^(nl)(\\mu), nf=nl+1 is the number
of active flavours and loops represents the number of loops.
fdelm (=1) can be omitted.";

 mPS2mMS::usage =
    "mPS2mMS[mPS, mq, asnlmu, Mu, Muf, nl, loops, fdelm] computes
mMS(Mu) from mPS(Muf). asnmu = \\alpha_s^(nl)(\\mu).
fdelm (=1) can be omitted.";

 mPS2mSI::usage =
    "mPS2mSI[mPS_, mq_, asfct_, Muf_, nl_, loops_, fdelm_] computes
mMS(mMS) from mPS(Muf). asfct head of function which computes
\\alpha_s^{(nl)}.
fdelm (=1) can be omitted.";

 mOS2mRS::usage =
    "mOS2mRS[mOS, mq, asmu, mu, nuf, nl, loops]  computes the RS
mass mRS(nuf) at the scale mu. {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nf)(\\mu), nf=nl+1 is the number
of active flavours and loops represents the number of loops.";

 mMS2mRS::usage =
    "mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, fdelm]
computes mRS(nuf) from mMS(Mu). asmu = \\alpha_s^(nl)(Mu).
fdelm (=1) can be omitted.";

 mRS2mMS::usage =
    "mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, fdelm]
computes mMS(Mu) from mRS(nuf). asnlmu = \\alpha_s^(nl)(Mu).
fdelm (=1) can be omitted.";

 mRS2mSI::usage =
    "mRS2mSI[mRS, mq, asfct, nuf, nl, loops, fdelm]
computes mMS(mMS) from mRS(nuf). asfct is the head of the 
function which computes \\alpha_s^(nl).
fdelm (=1) can be omitted.";

 mOS2mRSp::usage =
    "As mOS2mRS[] but for RS' mass.";

 mMS2mRSp::usage =
    "As mMS2mRS[] but for RS' mass.";

 mRSp2mMS::usage =
    "As mRS2mMS[] but for RS' mass."; 

 mRSp2mSI::usage =
    "As mRS2mSI[] but for RS' mass.";

 mOS2m1S::usage =
    "mOS2m1S[mOS, mq, asmu, mu, nl, loops] computes the 1S 
mass at the scale mu. {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nf)(\\mu), nf=nl+1 is the number
of active flavours and loops represents the number of loops.";

 mMS2m1S::usage =
    "mMS2m1S[mMS, mq, asmu, mu, nl, loops, fdelm]
computes m1S from mMS(mu). asmu = \\alpha_s^(nl)(mu).
fdelm (=1) can be omitted.";

 m1S2mMS::usage =
    "m1S2mMS[m1S, mq, asnlmu, Mu, nl, loops, fdelm]
compues mMS(Mu) from m1S. asnlmu = \\alpha_s^(nl)(Mu).
fdelm (=1) can be omitted.";

 m1S2mSI::usage =
    "m1S2mSI[m1S, mq, asfct, nl, loops, fdelm]
computes mMS(mMS) from m1S. asfct is the head of the 
function which computes \\alpha_s^(nl).
fdelm (=1) can be omitted.";

 DecAsUpOS::usage = 
    "DecAsUpOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the on-shell
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownOS::usage = 
    "DecAsDownOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpOS[als ,massth ,muth ,nl ,l ]";

 DecAsUpMS::usage = 
    "DecAsUpMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the MS-bar
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownMS::usage = 
    "DecAsDownMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpMS[als ,massth ,muth ,nl ,l ]";

 DecAsUpSI::usage = 
    "DecAsUpSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the
scale invariant
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownSI::usage = 
    "DecAsDownSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpSI[als ,massth ,muth ,nl ,l ]";

 DecMqUpOS::usage = 
    "DecMqUpOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the on-shell mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownOS::usage = 
    "DecMqDownOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpOS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpMS::usage = 
    "DecMqUpMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the MS-bar mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownMS::usage = 
    "DecMqDownMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpMS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpSI::usage = 
    "DecMqUpSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the scale invariant mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownSI::usage = 
    "DecMqDownSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpSI[mq, als ,massth ,muth ,nl ,l ]";

 DecLambdaUp::usage = 
    "DecLambdaUp[lam, massth, nl ,l ] computes \\Lambda^(nl+1) 
from the knowledge of lam=\\Lambda^(nl).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 DecLambdaDown::usage = 
    "DecLambdaDown[lam, massth, nl ,l ] computes \\Lambda^(nl) 
from the knowledge of lam=\\Lambda^(nl+1).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 AlL2AlH::usage = 
    "AlL2Alh[als, mu1, decpar, mu2, nloops] computes alphas(mu2) with h active
flavours from the knowledge of als=alphas(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsUpOS[..].
    ";

 AlH2AlL::usage = 
    "AlH2AlL[als, mu1, decpar, mu2, nloops] computes alphas(mu2) with l active
flavours from the knowledge of als=alphas(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsDownOS[..].
    ";

 mL2mH::usage = 
    "mL2mH[mql, asl, mu1, decpar, mu2, l] computes mq(mu2) with h active
flavours from the knowledge of mql=mq(mu1) with l active flavours.
asl=\\alpha_s(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqUpOS[..] and DecAsUpOS[..].
    ";

 mH2mL::usage = 
    "mH2L[mqh ,ash ,mu1 , decpar, mu2, l] computes mq(mu2) with l active
flavours from the knowledge of mqh=mq(mu1) with h active flavours.
ash=\\alpha_s(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqDownOS[..] and DecAsDownOS[..].
    ";

 AsRunDec::usage = 
    "AsRunDec[als ,mu0 ,mu , l] computes \\alpha_s(mu) from the knowledge of
als=\\alpha_s(mu0) The number of active flavours is determined automatically
from the choices of mu and mu0, respectively. The decoupling is performed at 
the pole mass of the respective heavy quark.";

 AsmMSrunexact::usage = 
    "AsmMSrunexact[mmu0, asmu0, mu0, mu, nf, l] solves simultanesously
the renormlaization group functions for \\alpha_s and mq.
Computes \\alpha_s(mu) and mq(mu) using asmu0=\\alpha_s(mu0) and 
mmu0 = mq(mu0) as starting values.";

 Mc5Mzfrommuc4::usage = 
    "";

 AsMbfromAsMz::usage = 
    "";

 As4McfromAs5Mz::usage = 
    "";

 setbeta::usage = 
    "Coefficients of QCD beta function."

(* ************************************************************ *)

Begin["`Modules`"];

cut[x_, n_Integer] := (x^p_Integer /; p > n -> 0);

(*
 argument used in N[] 
*)

$NumPrec=20;
numprec=$NumPrec;

(*
 minimal allowed ratio of \mu/\Lambda for which no warning is printed
*)
rmulam = 1.5;

(*
 numerical expression of some symbols ... 
*)

num1 = {
    z2 -> Zeta[2],
    z3 -> Zeta[3],
    z4 -> Zeta[4],
    z5 -> Zeta[5],
    log2 -> Log[2],
    B4->16*PolyLog[4,1/2]+2/3*Log[2]^4-2/3*Pi^2*Log[2]^2-13/180*Pi^4
};

(*
 numerical expression for colour factors (QCD)
*)

cf2num = {
    cf -> 4/3,
    ca -> 3,
    tr -> 1/2
};

(*
 coefficients of \beta function 
*)

setbeta = {
    b0 -> 11/4 - nf/6, 
    b1 -> 51/8 - (19*nf)/24,
    b2 -> 2857/128 - (5033*nf)/1152 + (325*nf^2)/3456,
    b3 -> 149753/1536 - (1078361*nf)/41472 + (50065*nf^2)/41472 +
	(1093*nf^3)/186624 + (891*Zeta[3])/64 - (1627*nf*Zeta[3])/1728 +
	    (809*nf^2*Zeta[3])/2592
	    };

(*
 coefficients of \gamma_m function 
*)

setgamma = {
    g0 -> 1,
    g1 -> (101/2 - (5*nf)/3)/12,
    g2 -> (3747/4 - (nf*(1108/9 + (70*nf)/27 + 80*Zeta[3]))/2)/48,
    g3 -> 4603055/41472 - (91723*nf)/6912 + (2621*nf^2)/31104 - 
	(83*nf^3)/15552 +
     (11*nf*Pi^4)/288 - (nf^2*Pi^4)/432 + (530*Zeta[3])/27 -
	 (2137*nf*Zeta[3])/144 + (25*nf^2*Zeta[3])/72 + (nf^3*Zeta[3])/108 -
	     (275*Zeta[5])/8 + (575*nf*Zeta[5])/72,
    g4 -> (* 559.7069 - 143.6864*nf + 7.4824*nf^2 + 0.1083*nf^3 - 0.000085359*nf^4 *)
          (1/4)^5*(99512327/162 + 46402466/243*Zeta[3] + 96800*Zeta[3]^2-698126/9*Zeta[4]
            -231757160/243*Zeta[5] + 242000*Zeta[6] + 412720*Zeta[7]
          +nf*(-150736283/1458 - 12538016/81*Zeta[3] - 75680/9*Zeta[3]^2 + 2038742/27*Zeta[4] 
               +49876180/243*Zeta[5] - 638000/9*Zeta[6] - 1820000/27*Zeta[7])
          +nf^2*(1320742/729 + 2010824/243*Zeta[3] + 46400/27*Zeta[3]^2 - 166300/27*Zeta[4] -264040/81*Zeta[5] + 92000/27*Zeta[6])
          +nf^3*(91865/1458 + 12848/81*Zeta[3] +448/9*Zeta[4] - 5120/27*Zeta[5]) 
          +nf^4*(-260/243 - 320/243*Zeta[3] + 64/27*Zeta[4]))    
    };

(* *** analytic result for g4:
(99512327/162 - (349063*Pi^4)/405 + (48400*Pi^6)/189 + 
  nf^4*(-260/243 + (32*Pi^4)/1215 - (320*Zeta[3])/243) + 
  (46402466*Zeta[3])/243 + 96800*Zeta[3]^2 + 
  nf^2*(1320742/729 - (16630*Pi^4)/243 + (18400*Pi^6)/5103 + 
    (2010824*Zeta[3])/243 + (46400*Zeta[3]^2)/27 - (264040*Zeta[5])/81) + 
  nf^3*(91865/1458 + (224*Pi^4)/405 + (12848*Zeta[3])/81 - 
    (5120*Zeta[5])/27) - (231757160*Zeta[5])/243 + 
  nf*(-150736283/1458 + (1019371*Pi^4)/1215 - (127600*Pi^6)/1701 - 
    (12538016*Zeta[3])/81 - (75680*Zeta[3]^2)/9 + (49876180*Zeta[5])/243 - 
    (1820000*Zeta[7])/27) + 412720*Zeta[7])/1024
*** *)

(*
 \alpha_s as function of \mu and \Lambda expanded in 1/\Lambda\beta_0
*)

setasL = {
    as0 -> 1/L/b0,
    as1 -> 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L],
    as2 -> ( 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L] +
	    1/(b0*L)^3*(b1^2/b0^2*(Log[L]^2-Log[L]-1)+b2/b0)
	    ),
    as3 -> (
	    1/(b0*L) - (b1*Log[L])/(b0^3*L^2) + 
	    (b2/b0^4 + (b1^2*(-1 - Log[L] + Log[L]^2))/b0^5)/L^3 + 
	    (b3/(2*b0^5) - (3*b1*b2*Log[L])/b0^6 + 
	     (b1^3*(-1/2 + 2*Log[L] + (5*Log[L]^2)/2 - Log[L]^3))/b0^7)/L^4
	    )
    };


(*
 \zeta_g, 1/\zeta_g, ...
 as6to5:  api^(nf) = api^(nf-1)*(1 + O(api^(nf-1)) + ...)
*)

sa4 = PolyLog[4,1/2];
sa5 = PolyLog[5,1/2];

setdec = {

    as6to5os ->
  1 + (7*api^2)/24 + (58933*api^3)/124416 + (api*lmm)/6 + (19*api^2*lmm)/24 +
   (8941*api^3*lmm)/1728 + (api^2*lmm^2)/36 + (511*api^3*lmm^2)/576 +
   (api^3*lmm^3)/216 - (2479*api^3*nl)/31104 - (409*api^3*lmm*nl)/1728 +
   (2*api^3*z2)/3 - (api^3*nl*z2)/9 + (2*api^3*z2*Log[2])/9 +
   (80507*api^3*Zeta[3])/27648

+ 
 api^4*(39.754401025244204746144247755796743`19.7659 + (47039*lmm^2)/3456 + 
   (14149*lmm^3)/10368 + lmm^4/1296 - (2965705*sa4)/54432 + (121*sa5)/36 + 
   (697121*z2)/19440 + (49*log2*z2)/54 + (243892631*z4)/8709120 - 
   (605*z2*z4)/2688 - (330575*z5)/41472 - (587*z2*Log[2])/81 + 
   (2057*z4*Log[2])/576 + (2699593*z2*Log[2]^2)/217728 + 
   (121*z2*Log[2]^3)/432 + nl*(-1773073/746496 - (9115*lmm^2)/10368 - 
     (107*lmm^3)/1728 - (173*sa4)/5184 - (557*z2)/162 + (697709*z4)/165888 - 
     (115*z5)/576 - (22*z2*Log[2])/81 + (1709*z2*Log[2]^2)/20736 - 
     (173*Log[2]^4)/124416 + lmm*(-1140191/373248 - (47*z2)/54 - 
       (2*log2*z2)/27 - (132283*Zeta[3])/82944) - (4756441*Zeta[3])/995328) + 
   nl^2*(140825/1492992 + (493*lmm^2)/20736 + lmm*(1679/186624 + z2/27) + 
     (13*z2)/162 + (19*Zeta[3])/1728) - (1439*z2*Zeta[3])/216 + 
   lmm*(21084715/746496 + (35*z2)/9 + (35*log2*z2)/27 + 
     (2922161*Zeta[3])/165888))

,

    as6to5ms ->
  1 + (api*lmm)/6 + api^2*(-11/72 + (11*lmm)/24 + lmm^2/36) + 
   api^3*(-564731/124416 + (2645*lmm)/1728 + (167*lmm^2)/576 + lmm^3/216 + 
   (2633/31104 - (67*lmm)/576 + lmm^2/36)*nl + (82043*Zeta[3])/27648),

    as6to5si ->
  1 + (api*lmmu)/6 + api^2*(-11/72 + (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(-564731/124416 + (2191*lmmu)/576 + (511*lmmu^2)/576 + lmmu^3/216 + 
   (2633/31104 - (281*lmmu)/1728)*nl + (82043*Zeta[3])/27648),

    as5to6os ->
  1 - (api*lmm)/6 + api^2*(-7/24 - (19*lmm)/24 + lmm^2/36) +
   api^3*(-58933/124416 - (8521*lmm)/1728 - (131*lmm^2)/576 - lmm^3/216 +
      nl*(2479/31104 + (409*lmm)/1728 + z2/9) - (2*z2)/3 - (2*z2*Log[2])/9 -
      (80507*Zeta[3])/27648)

       +
 api^4*(-38.60084204117141666826686318045613207152`19.74447411993687 + 
   (3179149*sa4)/54432 + (121*sa5)/36 - (7693*lmm^2)/1152 - 
   (8371*lmm^3)/10368 + lmm^4/1296 - (697121*z2)/19440 - 
   (128623*z2^2)/69120 - (121*z2^3)/1344 - (243892631*z4)/8709120 + 
   (605*z2*z4)/2688 + (49309*z5)/20736 + (1027*z2*Log[2])/162 + 
   (2057*z2^2*Log[2])/720 - (2057*z4*Log[2])/576 - 
   (2913037*z2*Log[2]^2)/217728 + (121*z2*Log[2]^3)/432 + 
   lmm*(-19696909/746496 - (29*z2)/9 - (29*z2*Log[2])/27 - 
     (2439119*Zeta[3])/165888) + nl^2*(-140825/1492992 - (493*lmm^2)/20736 + 
     lmm*(-1679/186624 - z2/27) - (13*z2)/162 - (19*Zeta[3])/1728) + 
   (1439*z2*Zeta[3])/216 + nl*(1773073/746496 + (173*sa4)/5184 + 
     (6661*lmm^2)/10368 + (107*lmm^3)/1728 + (557*z2)/162 - 
     (697709*z4)/165888 + (115*z5)/576 + (22*z2*Log[2])/81 - 
     (1709*z2*Log[2]^2)/20736 + (173*Log[2]^4)/124416 + 
     (4756441*Zeta[3])/995328 + lmm*(1110443/373248 + (41*z2)/54 + 
       (2*z2*Log[2])/27 + (132283*Zeta[3])/82944)))
,

    as5to6ms ->
  1 - (api*lmm)/6 + api^2*(11/72 - (11*lmm)/24 + lmm^2/36) +
   api^3*(564731/124416 - (955*lmm)/576 + (53*lmm^2)/576 - lmm^3/216 +
	(-2633/31104 + (67*lmm)/576 - lmm^2/36)*nl - (82043*Zeta[3])/27648),

    as5to6si ->
  1 - (api*lmmu)/6 + api^2*(11/72 - (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(564731/124416 - (6793*lmmu)/1728 - (131*lmmu^2)/576 - lmmu^3/216 + 
   (-2633/31104 + (281*lmmu)/1728)*nl - (82043*Zeta[3])/27648),

    mq5to6os ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) +
   api^3*(1871/2916 - B4/36 + (121*lmm)/2592 + (319*lmm^2)/432 +
      (29*lmm^3)/216 + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 -
         (2*z3)/27) - (407*z3)/864 - (5*lmm*z3)/6 + (5*z4)/4)
   + api^4*(
   0.22333826337861673 + 2.673908502292928*lmm1OS +
   6.226671960090351*lmm1OS^2 + 2.1651234567901234*lmm1OS^3 +
   0.2647569444444444*lmm1OS^4 +
   (-1.5035995976469407 - 0.6469697615283445*lmm1OS -
     0.9260223765432098*lmm1OS^2 - 0.16319444444444445*lmm1OS^3 -
     0.034722222222222224*lmm1OS^4)*nl +
   (0.05615582122474318 + 0.005015862378999767*lmm1OS +
     0.023919753086419752*lmm1OS^2 + 0.0011574074074074073*lmm1OS^4)*nl^2
   ) /. {lmm1OS->lmm},

    mq5to6ms ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) + 
   api^3*(2951/2916 - B4/36 + (175*lmm^2)/432 + (29*lmm^3)/216 + 
   lmm*(-311/2592 - (5*z3)/6) + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 - 
     (2*z3)/27) - (407*z3)/864 + (5*z4)/4),

    mq5to6si ->
  1 + api^2*(89/432 - (5*lmmu)/36 + lmmu^2/12) + 
   api^3*(2951/2916 - B4/36 - (1031*lmmu)/2592 + (319*lmmu^2)/432 + 
   (29*lmmu^3)/216 + nl*(1327/11664 - (53*lmmu)/432 - lmmu^3/108 - 
     (2*z3)/27) - (407*z3)/864 - (5*lmmu*z3)/6 + (5*z4)/4),

    mq6to5os ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) +
   api^3*(-1871/2916 + B4/36 - (299*lmm)/2592 - (299*lmm^2)/432 -
      (35*lmm^3)/216 - (1327*nl)/11664 + (53*lmm*nl)/432 + (lmm^3*nl)/108 +
      (407*z3)/864 + (5*lmm*z3)/6 + (2*nl*z3)/27 - (5*z4)/4)
   + api^4 * (
   -0.30107210254141137 - 3.714941523611081*lmm1OS -
   5.541401336860276*lmm1OS^2 - 2.677854938271605*lmm1OS^3 -
   0.33188657407407407*lmm1OS^4 +
   (1.5035995976469416 + 0.6346059568442685*lmm1OS +
     0.9873649691358025*lmm1OS^2 + 0.16319444444444445*lmm1OS^3 +
     0.03935185185185185*lmm1OS^4)*nl +
   (-0.05615582122474318 - 0.005015862378999758*lmm1OS -
     0.023919753086419752*lmm1OS^2 - 0.0011574074074074073*lmm1OS^4)*nl^2
   /. {lmm1OS -> lmm} ),

    mq6to5ms ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) + 
   api^3*(-2951/2916 + B4/36 - (155*lmm^2)/432 - (35*lmm^3)/216 + 
   nl*(-1327/11664 + (53*lmm)/432 + lmm^3/108 + (2*z3)/27) + 
   lmm*(133/2592 + (5*z3)/6) + (407*z3)/864 - (5*z4)/4),

    mq6to5si ->
  1 + api^2*(-89/432 + (5*lmmu)/36 - lmmu^2/12) + 
   api^3*(-2951/2916 + B4/36 + (853*lmmu)/2592 - (299*lmmu^2)/432 - 
   (35*lmmu^3)/216 + nl*(-1327/11664 + (53*lmmu)/432 + lmmu^3/108 + 
     (2*z3)/27) + (407*z3)/864 + (5*lmmu*z3)/6 - (5*z4)/4)

   };

(ddd4 = b1b^4*g0b/4 - 3*b1b^2*b2b*g0b/4 + b2b^2*g0b/4 
 + b1b*b3b*g0b/2 - b4b*g0b/4 - b1b^3*g1b/4 + b1b*b2b*g1b/2
 - b3b*g1b/4 + b1b^2*g2b/4 - b2b*g2b/4 - b1b*g3b/4 + g4b/4
 /. {g0b->g0/b0, g1b->g1/b0, g2b->g2/b0, g3b->g3/b0, g4b->g4/b0,
     b1b->b1/b0, b2b->b2/b0, b3b->b3/b0, b4b->Global`b4num/b0});

ddd1 = g1/b0 - b1*g0/b0^2;
ddd2 = 1/2*(+ g2/b0 + b1^2*g0/b0^3 - b1*g1/b0^2 - b2*g0/b0^2);
ddd3 = 1/3*( g3/b0 - b1^3*g0/b0^4 + 2*b1*b2*g0/b0^3
	    - b3*g0/b0^2 + b1^2*g1/b0^3 - b2*g1/b0^2 - b1*g2/b0^2);

setMS2OS = {

  ms2msc ->  1 
      + ddd1 * x
      + ( 1/2 * ddd1^2 + ddd2 ) * x^2
      + ( 1/6*ddd1^3 + ddd1*ddd2 + ddd3) * x^3
      + ( ddd1^4/24 + ddd1^2*ddd2/2 + ddd2^2/2 + ddd1*ddd3 + ddd4) * x^4 ,

    msfromos -> 1 + 
	api*(-cf - (3*cf*lmM)/4) + api^2*((-1111*ca*cf)/384 + (7*cf^2)/128 -
   (185*ca*cf*lmM)/96 + (21*cf^2*lmM)/32 - (11*ca*cf*lmM^2)/32 +
   (9*cf^2*lmM^2)/32 + (143*cf*tr)/96 + (13*cf*lmM*tr)/24 + (cf*lmM^2*tr)/8 +
   (71*cf*nl*tr)/96 + (13*cf*lmM*nl*tr)/24 + (cf*lmM^2*nl*tr)/8 +
   (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 + 3*cf^2*log2*z2 -
   cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 - (3*cf^2*z3)/4) +
 api^3*((lmM^2*(-2341*ca^2*cf + 1962*ca*cf^2 - 243*cf^3 + 1492*ca*cf*tr -
      468*cf^2*tr + 1492*ca*cf*nl*tr - 468*cf^2*nl*tr - 208*cf*tr^2 -
      416*cf*nl*tr^2 - 208*cf*nl^2*tr^2))/1152 +
   (lmM^3*(-242*ca^2*cf + 297*ca*cf^2 - 81*cf^3 + 176*ca*cf*tr -
      108*cf^2*tr + 176*ca*cf*nl*tr - 108*cf^2*nl*tr - 32*cf*tr^2 -
      64*cf*nl*tr^2 - 32*cf*nl^2*tr^2))/1152 +
   (lmM*(-105944*ca^2*cf + 52317*ca*cf^2 - 13203*cf^3 + 74624*ca*cf*tr -
      5436*cf^2*tr + 55616*ca*cf*nl*tr + 2340*cf^2*nl*tr - 12608*cf*tr^2 -
      18304*cf*nl*tr^2 - 5696*cf*nl^2*tr^2 + 12672*ca^2*cf*z2 -
      52704*ca*cf^2*z2 + 19440*cf^3*z2 - 38016*ca^2*cf*log2*z2 +
      91584*ca*cf^2*log2*z2 - 31104*cf^3*log2*z2 - 29952*ca*cf*tr*z2 +
      27648*cf^2*tr*z2 + 13824*ca*cf*log2*tr*z2 - 27648*cf^2*log2*tr*z2 +
      8064*ca*cf*nl*tr*z2 + 12096*cf^2*nl*tr*z2 + 13824*ca*cf*log2*nl*tr*z2 -
      27648*cf^2*log2*nl*tr*z2 + 9216*cf*tr^2*z2 + 4608*cf*nl*tr^2*z2 -
      4608*cf*nl^2*tr^2*z2 + 9504*ca^2*cf*z3 - 22896*ca*cf^2*z3 +
      7776*cf^3*z3 + 6912*ca*cf*tr*z3 - 3456*cf^2*tr*z3 +
      6912*ca*cf*nl*tr*z3 - 3456*cf^2*nl*tr*z3))/13824
	),

    msfromos3set[0] -> -202,
    msfromos3set[1] -> -176,
    msfromos3set[2] -> -150,
    msfromos3set[3] -> -126,
    msfromos3set[4] -> -103,
    msfromos3set[5] -> -82,

    msfromos3 -> -9478333/93312 - (644201*Pi^2)/38880 + (695*Pi^4)/7776 + 
 (587*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + (55*Log[2]^4)/162 + 
 (220*PolyLog[4, 1/2])/27 + nl^2*(-2353/23328 - (13*Pi^2)/324 - 
   (7*Zeta[3])/54) - (61*Zeta[3])/27 + (1439*Pi^2*Zeta[3])/432 + 
 nl*(246643/23328 + (967*Pi^2)/648 - (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - 
   (2*Pi^2*Log[2]^2)/81 - Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + 
   (241*Zeta[3])/72) - (1975*Zeta[5])/216,


    osfromms -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfrommsset -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) - (ca*cf^2*z2)/4 +
   (15*cf^3*z2)/16 + (3*ca*cf^2*log2*z2)/4 - (3*cf^3*log2*z2)/2 +
   (cf^2*tr*z2)/2 - (cf^2*nl*tr*z2)/4 - (3*ca*cf^2*z3)/16 + (3*cf^3*z3)/8 +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfromms3set[0] -> 194,
    osfromms3set[1] -> 168,
    osfromms3set[2] -> 143,
    osfromms3set[3] -> 119,
    osfromms3set[4] -> 96,
    osfromms3set[5] -> 75,

    osfromms3 -> 8481925/93312 + 
 (137*nl)/216 + (652841*Pi^2)/38880 - (nl*Pi^2)/27 - 
 (695*Pi^4)/7776 - (575*Pi^2*Log[2])/162 - (22*Pi^2*Log[2]^2)/81 - 
 (55*Log[2]^4)/162 - (220*PolyLog[4, 1/2])/27 - 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) + (58*Zeta[3])/27 - 
 (1439*Pi^2*Zeta[3])/432 - nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) + 
 (1975*Zeta[5])/216,

    mumfromos -> 1 
   - apiM*cf + apiM^2*((-1111*ca*cf)/384 + (199*cf^2)/128 + (143*cf*tr)/96 +
   (71*cf*nl*tr)/96 + (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 +
   3*cf^2*log2*z2 - cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 -
   (3*cf^2*z3)/4),

    mumfromos3set[0] -> -170,
    mumfromos3set[1] -> -146,
    mumfromos3set[2] -> -123,
    mumfromos3set[3] -> -101,
    mumfromos3set[4] -> -81,
    mumfromos3set[5] -> -62,

    mumfromos3 -> -7172965/93312 - 
	(293*nl)/216 - (618281*Pi^2)/38880 - (nl*Pi^2)/9 + 
 (695*Pi^4)/7776 + (623*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + 
 (55*Log[2]^4)/162 + (220*PolyLog[4, 1/2])/27 + 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) - (70*Zeta[3])/27 + 
 (1439*Pi^2*Zeta[3])/432 + nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) - 
 (1975*Zeta[5])/216,

    msfromri -> 1 - (4*api)/3 + api^2*(-995/72 + (89*nf)/144 + (19*z3)/6) + 
 api^3*(-6663911/41472 + (118325*nf)/7776 - (4459*nf^2)/23328 + 
   (5*nf*z4)/12 - (185*z5)/36 + (408007*z3)/6912 - 
   (617*nf*z3)/216 - (nf^2*z3)/54),

    rifromms -> 1 + (4*api)/3 + api^2*(1123/72 - (89*nf)/144 - (19*z3)/6) + 
 api^3*(6663911/41472 - (118325*nf)/7776 + (4459*nf^2)/23328 + 
   (4*(1123/72 - (89*nf)/144 - (19*z3)/6))/3 - (408007*z3)/6912 + 
   (617*nf*z3)/216 + (nf^2*z3)/54 - (4*(-995/72 + (89*nf)/144 + (19*z3)/6))/
    3 - (5*nf*z4)/12 + (185*z5)/36)

    };

(* ************************************************************ *)

(* Compute \\Lambda^(nf) using 3 methods:
   1. Explicit formulae for \\Lambda: lamexpl[alphas_,mu_,nn_,loops_]
   2. Implicit formulae for \\Lambda: lamimpl[alphas_,mu_,nn_,loops_]
   3. Solving the integral exactly:   lamexact[alphas_,mu_,nn_]
*)

(* Determine \Lamba for given \alpha_s				*)

LamExpl[alphas_,mu_,nn_,loops_] := Module[
    {logmu2lam2,nnn},
    If[(loops<1)||(loops>4), Print["LamExpl: Invalid # loops."];
       Return[];];
    (* \log (\mu^2/\Lambda^2) =	
     *)
    logmu2lam2 = 1/(as*b0) + (as*(-b1^2 + b0*b2))/b0^3 +
	(as^2*(b1^3 - 2*b0*b1*b2 + b0^2*b3))/(2*b0^4) + (b1*Log[as])/b0^2 +
	    (b1*Log[b0])/b0^2;
    (* \Lambda/\mu = Exp[ lamomu[nn_] ]; nn: # active flavours
     *)
    lamomu[nnn_]:=Collect[Expand[ 
	-1/2*logmu2lam2/.setbeta/.nf->nnn
	],as]/.as:>ToExpression["api"<>ToString[nnn]];
    Return[N[mu*Exp[
	Expand[
	    (Expand[ lamomu[nn] * ToExpression["api"<>ToString[nn]]^2 ]
	     /.cut[ToExpression["api"<>ToString[nn]],loops]
	     )
		*1/ToExpression["api"<>ToString[nn]]^2
		]
	/.ToExpression["api"<>ToString[nn]]->alphas/Pi] 
	     ,numprec]];
];

(* ************************************************************ *)

LamImpl[alphas_,mu_,nn_,loops_] := Module[
    {lam,astmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==1, astmp=as0];
    If[loops==2, astmp=as1];
    If[loops==3, astmp=as2];
    If[loops==4, astmp=as3];
    Return[lam/.FindRoot[
	(Expand[
	    astmp/.setasL
		/.setbeta/.nf->nn]/.L->Log[mu^2/lam^2])==alphas/Pi,
	    {lam,LamExpl[alphas,mu,nn,loops]} 
    ]
       ]
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Input: \Lambda, \mu, # active flavours and # loops
*)

AlphasLam[lambda_,mu_,nn_,loops_] := Module[
    {},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[ N[mu/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[mu/lambda]," is very small!"];
    ];
    If[loops==1, Return[ N[Pi*as0/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==2, Return[ N[Pi*as1/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==3, Return[ N[Pi*as2/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==4, Return[ N[Pi*as3/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Solve differential equation numerically
   Input: \alpha_s(\mu0), \mu0, \mu, # active flavours and # loops
*)

AlphasExact[alphasmu0_,mu0_,muend_,nnf_,loops_] := Module[
    {rhs,solx,a,x,lambda,alphasmu0tmp},
    alphasmu0tmp = Rationalize[alphasmu0,10^-numprec];
    If[(loops<1)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    lambda = LamExpl[alphasmu0,mu0,nnf,4];
    If[ N[muend/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[muend/lambda]," is very small!"];
    ];
    If[ mu0 == muend , Return[alphasmu0] ];
    If[ loops==1, rhs = -a[x]^2*(b0) ];
    If[ loops==2, rhs = -a[x]^2*(b0 + b1*a[x]) ];
    If[ loops==3, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2) ];
    If[ loops==4, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2 + b3*a[x]^3) ];
    If[ loops==5, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2 + b3*a[x]^3 
				 + b4*a[x]^4);
	Print["b4 is set to ",Global`b4num /. nf->nnf //N]; ];
    sol = N[(Pi*a[x]/.
	     NDSolve[{a[mu0^2]==alphasmu0tmp/Pi,
		      x*a'[x] == rhs
         /.setbeta/.{b4->Global`b4num}/.nf->nnf}, 
         a[x], 
         {x,muend^2,mu0^2},        WorkingPrecision->26,
                                   AccuracyGoal->16,
                                   PrecisionGoal->16,
                                   MaxSteps->5000]
/.x->muend^2)[[1]],numprec];
Return[sol];
];

(* ************************************************************ *)

(*--#[ AsmMSrunexactOLD: *)

(* Solve simultaneously RGE for \alpha_s and m_q. *)
AsmMSrunexactOLD[mmu0_, asmu0_, mu0_, mu_, nf_, l_] := Module[
    {res,b0,b1,g0,g1,beta,gammam,x,xl,api},
    If[l>4,Print["AsmMSrunexactOLD: l>4 not (yet) implemented."];Abort[];];
    If[l<=0,Return[{mmu0,asmu0}]];
    beta   = Expand[-api[x]^2*(xl*b0+xl^2*api[x]*b1+
			       xl^3*api[x]^2*b2+xl^4*api[x]^3*b3)
		    ]/.cut[xl,l]/.{xl->1};
    gammam = Expand[-api[x]  *(xl*g0+xl^2*api[x]*g1+
			       xl^3*api[x]^2*g2+xl^4*api[x]^3*g3)
		    ]/.cut[xl,l]/.{xl->1};
    b0 = 11/4 - nf/6;
    b1 = 51/8 - (19*nf)/24;
    b2 = 2857/128 - (5033*nf)/1152 + (325*nf^2)/3456;
    (b3 = 149753/1536 - (1078361*nf)/41472 + (50065*nf^2)/41472 +
     (1093*nf^3)/186624 + (891*Zeta[3])/64 - (1627*nf*Zeta[3])/1728 +
     (809*nf^2*Zeta[3])/2592);
    g0 = 1;
    g1 = (101/2 - (5*nf)/3)/12;
    g2 = (3747/4 - (nf*(1108/9 + (70*nf)/27 + 80*Zeta[3]))/2)/48;
    (g3 = 4603055/41472 - (91723*nf)/6912 + (2621*nf^2)/31104 - 
     (83*nf^3)/15552 +
     (11*nf*Pi^4)/288 - (nf^2*Pi^4)/432 + (530*Zeta[3])/27 -
     (2137*nf*Zeta[3])/144 + (25*nf^2*Zeta[3])/72 + (nf^3*Zeta[3])/108 -
     (275*Zeta[5])/8 + (575*nf*Zeta[5])/72);
    
    res = NDSolve[ {x*api'[x]      == beta,
		    x/mq[x]*mq'[x] == gammam,
		    api[mu0^2]     == asmu0/Pi, 
		    mq[mu0^2]      == mmu0},    {api[x],mq[x]}, {x,mu^2,mu0^2}
(*,
			WorkingPrecision->26,
			AccuracyGoal->16,
			PrecisionGoal->16,
			MaxSteps->5000
*)
		    ];
    Return[{mq[x],Pi*api[x]}/.res[[1]]/.{x->mu^2}];
];

(*--#] *)

(* ************************************************************ *)

(*
 Light-mass corrections at order \alpha_s^2
*)

delta[mOS_,mq_] := Module[
    {i,del,delta,r},
    del=0;
    delta[r_] := If[(r<=1) && (r>=0),
		    Pi^2/8*r-0.597*r^2+0.230*r^3,
		    Print["\\Delta(",N[r],
                          ") IS CALLED; THE FUNCTION IS NOT IMPLEMENTED 
FOR ARGUMENTS OUTSIDE THE INTERVAL [0,1]."];
		    Abort[]
		];
    For[i=1,i<=Length[mq],i++,
	del=del+delta[mq[[i]]/mOS];
    ];
    Return[del];
];

(* ************************************************************ *)

(*
 m_MS = M_OS * ( 1 + ... )
*)

mOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_] := 
    mOS2mMS[mOS,mq,asmu,mu,nnf,loops,1];

mOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_,fdelm_] := Module[
    {extmp,ex4ltmp,Mu},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromos 
	    + api^2*(-4/3 * delta[mOS,mq])
	    (* + api^3*msfromos3set[nnf-1] /.setMS2OS *)
	    + api^3*msfromos3 /.setMS2OS
	    )/.cut[api,loops] /.{lmM->Log[mu^2/mOS^2],nl->nnf-1
				 } /.num1 /.cf2num;
   ];

    If[ loops == 4, 
	(ex4ltmp = 
 - 3654.15040757339*fdelm - 1524.2292266911543*Log[Mu^2/mOS^2] - 
 288.778291935394*Log[Mu^2/mOS^2]^2 - 32.54735725308642*Log[Mu^2/mOS^2]^3 - 
 1.85546875*Log[Mu^2/mOS^2]^4 + (-1 + nnf)^3*(0. + 0.678141025604516*fdelm + 
   0.3205521521864135*Log[Mu^2/mOS^2] + 0.0800290327210927*
    Log[Mu^2/mOS^2]^2 + 0.010030864197530864*Log[Mu^2/mOS^2]^3) + 
 (-1 + nnf)^2*(0. - 43.48241924867489*fdelm - 
   19.82672048099557*Log[Mu^2/mOS^2] - 4.482957520194182*Log[Mu^2/mOS^2]^2 - 
   0.5270061728395061*Log[Mu^2/mOS^2]^3 - 0.04108796296296297*
    Log[Mu^2/mOS^2]^4) + (-1 + nnf)*(0. + 756.9421565599532*fdelm + 
   330.1770776731065*Log[Mu^2/mOS^2] + 67.99849534415492*Log[Mu^2/mOS^2]^2 + 
   7.595293209876542*Log[Mu^2/mOS^2]^3 + 0.48119212962962954*
    Log[Mu^2/mOS^2]^4)
	 );
	extmp = extmp + api^4 * ex4ltmp /. {nl->nnf-1, Mu->mu} // Expand;
      ];

    Return[ N[ mOS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_OS = M_MS * ( 1 + ... )
*)

mMS2mOS[mMS_,mq_,asmu_,mu_,nnf_,loops_] := mMS2mOS[mMS,mq,asmu,mu,nnf,loops,1];
mMS2mOS[mMS_,mq_,asmu_,mu_,nnf_,loops_,fdelm_] := Module[
    {extmp,ex4ltmp,Mu},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( osfromms 
	    + api^2*(4/3 * delta[mMS,mq])
	    (* + api^3*osfromms3set[nnf-1] /.setMS2OS *)
	    + api^3*osfromms3 /.setMS2OS
          )/.cut[api,loops] /.{lmm->Log[mu^2/mMS^2],nl->nnf-1} /.num1 /.cf2num;
   ];
    If[ loops == 4, 
	(ex4ltmp = 
3567.602784989066*fdelm + 1727.2260148986106*fdelm*Log[(1.*mu^2)/mMS^2] + 
 409.2429990574718*fdelm*Log[(1.*mu^2)/mMS^2]^2 + 
 66.93663194444443*fdelm*Log[(1.*mu^2)/mMS^2]^3 + 
 8.056278935185185*fdelm*Log[(1.*mu^2)/mMS^2]^4 + 
 (-1 + nnf)^3*(-0.678141025604516*fdelm - 0.3205521521864134*fdelm*
    Log[(1.*mu^2)/mMS^2] - 0.0800290327210927*fdelm*Log[(1.*mu^2)/mMS^2]^2 - 
   0.010030864197530864*fdelm*Log[(1.*mu^2)/mMS^2]^3) + 
 (-1 + nnf)*(-745.7207145811878*fdelm - 358.29765085086774*fdelm*
    Log[(1.*mu^2)/mMS^2] - 87.39262571554698*fdelm*Log[(1.*mu^2)/mMS^2]^2 - 
   11.883873456790122*fdelm*Log[(1.*mu^2)/mMS^2]^3 - 
   1.2705439814814814*fdelm*Log[(1.*mu^2)/mMS^2]^4) + 
 (-1 + nnf)^2*(43.396250117985666*fdelm + 20.528466368867228*fdelm*
    Log[(1.*mu^2)/mMS^2] + 4.971905254812516*fdelm*Log[(1.*mu^2)/mMS^2]^2 + 
   0.6304012345679011*fdelm*Log[(1.*mu^2)/mMS^2]^3 + 
   0.06655092592592593*fdelm*Log[(1.*mu^2)/mMS^2]^4)
	 );
	extmp = extmp + api^4 * ex4ltmp /. {nl->nnf-1, Mu->mu} // Expand;
      ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS) and afterwards m_MS(mu).
*)

mOS2mMSrun[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {mum,asM,asmum},
    asM = AlphasExact[asmu,mu,mOS,nnf,loops];
    mum = mOS2mSI[mOS,mq,asM,nnf,loops];
    asmum = AlphasExact[asmu,mu,mum,nnf,loops];
    Return[ N[ mMS2mMS[mum,asmum,asmu,nnf,loops] , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mOSrun[mMS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];
    (* print mMS2mSI2 and \alpha_s(mMS2mSI2)
     Print[ mMS2mSI2/.res1 ];
     Print[ asmMS/.res1 ];
     *)
    Return[ N[ mMS2mOS[mMS2mSI2/.res1,mq,asmMS/.res1,mMS2mSI2/.res1,nnf,loops] 
	       , numprec ] ];
];

(* ************************************************************ *)

mOS2mMSit[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( osfromms 
	    + api^2*(4/3 * delta[mOS,mq])
	    (* + api^3*osfromms3set[nnf-1] /.setMS2OS *)
	    + api^3*osfromms3 /.setMS2OS
          )/.cut[api,loops] /.{lmm->Log[mu^2/mMS^2],nl->nnf-1} /.num1 /.cf2num;
   ];
    Return[ N[ mMS/.FindRoot[ mOS == mMS * (extmp/.api->asmu/Pi) , {mMS,mOS}]
	       , numprec ] ];
];

(* ************************************************************ *)

(*
 \mu_m = M_OS * ( 1 + ... )
*)

mOS2mSI[mOS_,mq_,asM_,nnf_,loops_] := mOS2mSI[mOS,mq,asM,nnf,loops,1];
mOS2mSI[mOS_,mq_,asM_,nnf_,loops_,fdelm_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( mumfromos 
	    + apiM^2*(-4/3 * delta[mOS,mq])
	    (* + apiM^3*mumfromos3set[nnf-1] /.setMS2OS *)
	    + apiM^3*mumfromos3 /.setMS2OS
          )/.cut[apiM,loops] /.{nl->nnf-1} /.num1 /.cf2num;
   ];
    If[ loops == 4, 
        (ex4ltmp =
fdelm*(-3214.227044839041 + 692.4809215366435*(-1 + nnf) - 
  41.95978562498058*(-1 + nnf)^2 + 0.678141025604516*(-1 + nnf)^3)
	 );
	extmp = extmp + apiM^4 * ex4ltmp /. {nl->nnf-1} // Expand;
      ];
    Return[ N[ mOS * (extmp/.apiM->asM/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Running MS-bar mass
*)

mMS2mMS[mmu0_,asmu0_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    If[loops==5, Print["b4 is set to ",Global`b4num /. nf->nnf //N] ];

    Return[ N[ mmu0 * (ccc/.x->asmu/Pi)/(ccc/.x->asmu0/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mSI[mMS_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];

    Return[ N[ mMS2mSI2/.res1 , numprec ] ];
];


(mMS2mSInew[mqmu_,asmu_,mu_,nfnum_,nolo_] := Module[
    {mcmc,i,imax,asmc},
    imax = 30;
    mcmc = mqmu;
    For[i=1,i<=imax,i++,
	asmc = AlphasExact[asmu ,mu ,mcmc ,nfnum ,nolo];
        mcmc = mMS2mMS[mqmu, asmu, asmc, nfnum, nolo];
(*         Print[mcmc];						  *)
    ];
    Return[mcmc];
]);

(* ************************************************************ *)

(*
 m_RI = m_MS * ( 1 + ... )
*)

mMS2mRI[mMS_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( rifromms 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_MS = m_RI * ( 1 + ... )
*)

mRI2mMS[mRI_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromri 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mRI * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)
(*
 m_RGI = m_MS * ( 1 + ... )
*)
(*
mMS2mRGI[1,0.1,10,5,1]
*)
mMS2mRGI[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mMS, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x->asmu/Pi) , numprec ] ];
];

mMS2mRGImod[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc,bet0},
    bet0 = 11/4 - nnf/6;
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x-> 2*bet0*asmu/Pi ) , numprec ] ];
(*                             ^^^^^^ *)
];


(* ************************************************************ *)

(*
 m_MS = m_RGI * ( 1 + ... )
*)

mRGI2mMS[mRGI_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mRGI * (ccc/.x->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf+1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsUpOS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5os/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsUpMS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5ms/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsUpSI[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5si/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf-1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsDownOS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6os/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsDownMS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6ms/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsDownSI[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6si/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
 the decoupling is performed at "loops-1" order
*)

DecMqUpOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decouplingof quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];


(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl+1)
 input: \Lambda^(nl), ...
*)

DecLambdaUp[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0m,ex1,ex2,ex3,ex4},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[(L*(1/2 - b0[-1 + nf]/b0times2[nf]) +
 ((-b1p[-1 + nf] + b1p[nf])*Log[L] +
   (-b1p[-1 + nf]^2 + b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] + c2[nf, -1 + nf] +
     (-b1p[-1 + nf]^2 + b1p[-1 + nf]*b1p[nf])*Log[L])/Lb0m +
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 + b3p[-1 + nf]/2 - b3p[nf]/2 +
     b1p[nf]*(b1p[-1 + nf]^2 - b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]) +
     c3[nf, -1 + nf] + (b1p[-1 + nf]^2*b1p[nf] - b1p[-1 + nf]*b1p[nf]^2 +
       b1p[-1 + nf]*(-b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]))*Log[L] +
     (b1p[-1 + nf]^3/2 - (b1p[-1 + nf]^2*b1p[nf])/2)*Log[L]^2)/Lb0m^2 +
  b1p[nf]*Log[b0[-1 + nf]/b0[nf]])/b0times2[nf]
		   )/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0m->(L*b0[nf-1])}
	      ];
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];

    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];

    exra2 = Expand[ exra  
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]->11/72}
		    /.{c3[mm_,nn_]->-564731/124416+82043/27648*Zeta[3]
		       +2633/31104*nn}
		][[1]];
    Return[ lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl)
 input: \Lambda^(nl+1), ...
*)

DecLambdaDown[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[(L*(1/2 - b0[nf]/b0times2[-1 + nf]) +
 ((b1p[-1 + nf] - b1p[nf])*Log[L] + (b1p[-1 + nf]^2 - b1p[nf]^2 -
     b2p[-1 + nf] + b2p[nf] + c2[-1 + nf, nf] +
     (b1p[-1 + nf]*b1p[nf] - b1p[nf]^2)*Log[L])/Lb0 +
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 - b3p[-1 + nf]/2 + b3p[nf]/2 +
     b1p[-1 + nf]*(b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]) +
     c3[-1 + nf, nf] + (-(b1p[-1 + nf]^2*b1p[nf]) + b1p[-1 + nf]*b1p[nf]^2 +
       b1p[nf]*(b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]))*Log[L] +
     (-(b1p[-1 + nf]*b1p[nf]^2)/2 + b1p[nf]^3/2)*Log[L]^2)/Lb0^2 +
   b1p[-1 + nf]*Log[b0[nf]/b0[-1 + nf]])/b0times2[-1 + nf]
		   )/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0->(L*b0[nf])}
	      ];
    
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];

    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];
    
    exra2 = Expand[ exra  /.{iLb0->1/(L*b0[nf])}
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]->11/72}
		    /.{c3[mm_,nn_]->-564731/124416+82043/27648*Zeta[3]
		       +2633/31104*nn}
		][[1]];
    
    Return[lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(l)(mu1), ...
 output: \alpha_s^(h)(mu2)
*)

AlL2AlH[asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=asl;
    muini=mu1;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(h)(mu1), ...
 output: \alpha_s^(l)(mu2)
*)

AlH2AlL[ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=ash;
    muini=mu1;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(l)(mu1), ...
 output: m_q^(h)(mu2)
*)

mL2mH[mql_,asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=asl;
    muini=mu1;
    mqini=mql;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	res4 = DecMqUpOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]],loops];
    Return[res6];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(h)(mu1), ...
 output: m_q^(l)(mu2)
*)

mH2mL[mqh_,ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=ash;
    muini=mu1;
    mqini=mqh;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res4 = DecMqDownOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]]-1,loops];
    Return[res6];
];

(* ************************************************************ *)

(*
 AsRunDec[als, mu0, mu, l] computes \alpha_s^(m)(mu) from the knowledge
 of \alpha_s^(n)(mu0)  
*)

AsRunDec[als_,mu0_,mu_,loops_] := Module[
    {m,n,kk,decpar},
    If[mu>=Global`Mt/.Global`NumDef,m=6,
       If[mu>=Global`Mb/.Global`NumDef,m=5,
	  If[mu>=Global`Mc/.Global`NumDef,m=4,
	     If[mu<Global`Mc/.Global`NumDef,m=3
		]]]];
    If[mu0>=Global`Mt/.Global`NumDef,n=6,
       If[mu0>=Global`Mb/.Global`NumDef,n=5,
	  If[mu0>=Global`Mc/.Global`NumDef,n=4,
	     If[mu0<Global`Mc/.Global`NumDef,n=3
		]]]];
    If[m==n,
       (* Return[{AlphasExact[als,mu0,mu,n,loops],m,n}]; *)
       Return[AlphasExact[als,mu0,mu,n,loops]];
   ];
    decpar = {};
    If[m>n,
       For[kk=n+1,kk<=m,kk++,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Return[{AlL2AlH[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlL2AlH[als,mu0,decpar,mu,loops]];
   ];
    If[m<n,
       For[kk=n,kk>=m+1,kk--,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Print[decpar]; *)
       (* Return[{AlH2AlL[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlH2AlL[als,mu0,decpar,mu,loops]];
   ];
];

(* ************************************************************ *)

(* 
AsmMSrunexact: solve simultaneously RGE for \alpha_s and m_q.
*)

AsmMSrunexact[mmu0_, asmu0_, mu0_, mu_, nnf_, l_] := Module[
    {res,beta,gammam,x,xl,api},
    If[l>5,Print["AsmMSrunexact: l>5 not implemented."];Abort[];];

    beta   = Expand[-api[x]^2*(xl*b0+xl^2*api[x]*b1+
			       xl^3*api[x]^2*b2+xl^4*api[x]^3*b3+xl^5*api[x]^4*b4)
		    ]/.cut[xl,l]/.{xl->1};
    gammam = Expand[-api[x]  *(xl*g0+xl^2*api[x]*g1+
			       xl^3*api[x]^2*g2+xl^4*api[x]^3*g3+xl^5*api[x]^4*g4)
		    ]/.cut[xl,l]/.{xl->1};
    (*Print[beta];
    Print[gammam];*)

    res = NDSolve[ {x*api'[x]      == beta /. setbeta  /.{b4->Global`b4num} /. nf->nnf ,
		    x/mq[x]*mq'[x] == gammam /. setgamma /. nf->nnf,
		    api[mu0^2]     == asmu0/Pi, 
		    mq[mu0^2]      == mmu0},    {api[x],mq[x]}, {x,mu^2,mu0^2}
(*			,
		   WorkingPrecision->26,
		   AccuracyGoal->16,
		   PrecisionGoal->16,
		   MaxSteps->5000 *)
		   ];
    Return[{mq[x],Pi*api[x]}/.res[[1]]/.{x->mu^2}];
];

(* ************************************************************ *)

(*
 mOS2mMS[mOS, nf, l] computes MS-bar mass corresponding to mOS
*)

mOS2mMS[mOS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mOS,loops];
    Return[mOS2mMS[mOS,{},as,mOS,nnf,loops]];
];

(* ************************************************************ *)

(*
 mMS2mOS[mOS, nf, l] computes the on-shell mass corresponding to mMS
*)

mMS2mOS[mMS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mMS,loops];
    Return[mMS2mOS[mMS,{},as,mMS,nnf,loops]];
];

(* ************************************************************ *)

(* 
 input: \mu_c^(4)
 output: m_c(MZ)^(5)
*)

Mc5Mzfrommuc4[asMz_,muc4_,Mb_,mub_,Mz_,loops_] := Module[
    {alsmuth,alsmuthp,alsmuc,mcthp,mcth,mcMZ},

    alsmuth = AlphasExact[asMz,Mz,mub,5,loops];
    alsmuthp = DecAsDownOS[alsmuth,Mb,mub,4,loops];
    alsmuc = AlphasExact[alsmuthp,mub,muc4,4,loops];

    mcthp = mMS2mMS[muc4,alsmuc,alsmuthp,4,loops];
    mcth = DecMqUpOS[mcthp,alsmuthp,Mb,mub,4,loops];
    mcMZ = mMS2mMS[mcth,alsmuth,asMz,5,loops];

(*    Print["alsmuc, alsmuthp, alsmuth: ", alsmuc," ",alsmuthp," ",alsmuth]; *)
(*    Print["mcth, mcthp, mcMZ: ", mcth," ",mcthp," ",mcMZ]; *)
    
    Return[mcMZ];
];

(* ************************************************************ *)

(* 
 compare running of "exact" vs "\lambda" 
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(5)(Mb)
 rem: if no value is given for Mb the default one is chosen
*)

AsMbfromAsMz[asMz_,Mb_,loops_] := Module[
    {res1,res2},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,Mb,5,loops];
    res2 = AlphasLam[LamImpl[asMz,Global`Mz/.Global`NumDef,5,loops],Mb,5,loops];
    Return[ {res1, res2}];
];

AsMbfromAsMz[asMz_,loops_] := AsMbfromAsMz[asMz,Global`Mb/.Global`NumDef,
					   loops];

(* ************************************************************ *)

(* 
 running-decoupling-running
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(4)(Mc)
*)

As4McfromAs5Mz[asMz_,Mb_,mub_,Mc_,loops_] := Module[
    {res1,res2,res3},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,mub,5,loops];
    res2 = DecAsDownOS[res1,Mb,mub,4,loops];
    res3 = AlphasExact[res2,mub,Mc,4,loops];
    Return[res3];
];

(* ************************************************************ *)

(* Additions August 2015: *)

(* ************************************************************ *)

(*
 mOS2mPS[mOS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] computes the PS mass from the OS mass
*)

mOS2mPS[mOS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] := 
    Module[{extmp, api}, If[loops < 0 || loops > 4, 
       Print["PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; If[loops == 0, extmp = 0, 
       extmp = Expand[deltamf[api, Muf, Mu, nl]] /. cut[api, loops]; ]; 
      Return[N[mOS - (extmp /. api -> asmu/Pi), numprec]]; ];
  
deltamf[api_, Muf_, Mu_, nl_] := 
    Muf*((4*api)/3 + api^2*(97/9 + nl*(-22/27 - (2*Log[Mu^2/Muf^2])/9) + 
        (11*Log[Mu^2/Muf^2])/3) + api^3*(33623/216 + 3*Pi^2 - (3*Pi^4)/16 + 
        (610*Log[Mu^2/Muf^2])/9 + (121*Log[Mu^2/Muf^2]^2)/12 + 
        nl^2*(157/243 + (22*Log[Mu^2/Muf^2])/81 + Log[Mu^2/Muf^2]^2/27) + 
        nl*(-7145/324 - (493*Log[Mu^2/Muf^2])/54 - (11*Log[Mu^2/Muf^2]^2)/9 - 
          (13*Zeta[3])/9) + (11*Zeta[3])/2) + 
      api^4*(3567.723056629293 + (7271*Log[Mu^2/Muf^2]^2)/24 + 
        (1331*Log[Mu^2/Muf^2]^3)/48 + nl^3*(-2951/4374 - 
          (157*Log[Mu^2/Muf^2])/486 - (11*Log[Mu^2/Muf^2]^2)/162 - 
          Log[Mu^2/Muf^2]^3/162) + nl*(-701.2303148875468 - 
          (8485*Log[Mu^2/Muf^2]^2)/144 - (121*Log[Mu^2/Muf^2]^3)/24 + 
          Log[Mu^2/Muf^2]*(-253189/864 - (3*Pi^2)/2 + (3*Pi^4)/32 - 
            (44*Zeta[3])/3)) + nl^2*(1751971/46656 + Pi^4/135 + 
          (773*Log[Mu^2/Muf^2]^2)/216 + (11*Log[Mu^2/Muf^2]^3)/36 + 
          Log[Mu^2/Muf^2]*(15355/864 + (13*Zeta[3])/18) + 
          (259*Zeta[3])/108) + Log[Mu^2/Muf^2]*(26125/18 + (99*Pi^2)/4 - 
          (99*Pi^4)/64 + (363*Zeta[3])/8)));

(* ************************************************************ *)

(*
 mMS2mPS[mMS, mq, asmu, mu, muf, nl, loops] computes the PS mass from the MSbar mass
*)

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, 1, "no"];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, marker_String] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, 1, marker];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, fdelm_/;NumericQ[fdelm]] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, fdelm, "no"];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, fdelm_, marker_String] := 
    Module[{extmp, api, delmuf, exmOS, exmOSprime, apiprime, exals, mark}, 
           If[ marker == "no", mark = 1, mark = "RunDecXMS2OS" ];
	   If[loops < 0 || loops > 4, 
	      Print["PROCEDURE IS NOT IMPLEMENTED FOR ",  loops, " LOOPS."]; 
	      Return[]]; 
	   If[loops == 0, extmp = mMS, 
	      delmuf = deltamf[api, Muf, Mu, nl]; 
	      exmOS = mMS2mOS[mMS, {}, api*Pi, Mu, nl + 1, loops, fdelm];
	      exals = (1*DecAsUpMS[apiprime*Pi, mMS, Mu, nl, 4])/(apiprime*Pi); 
	      exmOSprime = Normal[Series[
		  exmOS /. {api -> apiprime*exals}, {apiprime, 0, loops}]]; 
	      extmp = Expand[mark * (exmOSprime /. apiprime -> api) - delmuf] /. 
		  cut[api, loops]; ]; Return[N[extmp /. api -> asmu/Pi, numprec]]; 
       ];

(* ************************************************************ *)

(*
 mPS2mMS[] computes the MSbar mass from the PS mass
*)

mPS2mMS[mPS_, mq_, asnlmu_, Mu_, Muf_, nl_, loops_] := 
    mPS2mMS[mPS, mq, asnlmu, Mu, Muf, nl, loops, 1];
 
mPS2mMS[mPS_, mq_, asnlmu_, Mu_, Muf_, nl_, loops_, fdelm_] :=
    Module[{mtmp, mMS}, mtmp = mMS2mPS[mMS, mq, asnlmu, Mu, Muf, nl, loops, fdelm];
    Return[mMS /. FindRoot[mPS == mtmp, {mMS, mPS}]]; ];

(* ************************************************************ *)

(*
 mPS2mSI[] computes the SI mass from the PS mass
*)

mPS2mSI[mPS_, mq_, asfct_, Muf_, nl_, loops_] := 
    mPS2mSI[mPS, mq, asfct, Muf, nl, loops, 1];
 
mPS2mSI[mPS_, mq_, asfct_, Muf_, nl_, loops_, fdelm_] := 
    Module[{mMS1, mMS, acc}, mMS1 = 0; mMS = mPS; acc = 1.*^-6; 
      While[Abs[mMS1 - mMS] > acc, mMS1 = mMS; 
        mMS = mPS2mMS[mPS, mq, asfct[mMS1], mMS1, Muf, nl, loops, fdelm]; ]; 
      Return[mMS]; ];

(* ************************************************************ *)

(*
 mOS2mPS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME] computes the RS or RS' mass from the OS mass
*)

mOS2mRS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mOS2mRS[mOS, mq, asmu, Mu, nuf, nl, loops, "no"]

mOS2mRSp[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mOS2mRS[mOS, mq, asmu, Mu, nuf, nl, loops, "yes"]
 
mOS2mRS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_] := 
    Module[{extmp, api}, If[loops < 0 || loops > 4, 
       Print["PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; If[loops == 0, extmp = 0, 
       If[PRIME === "yes", extmp = exOS2RSp[api, Mu, nuf, nl]; , 
         extmp = exOS2RS[api, Mu, nuf, nl]; ]; 
        extmp = extmp /. cut[api, loops]; ]; 
      Return[N[mOS - (extmp /. api -> asmu/Pi), numprec]]; ];
 
exOS2RSp[aapi_, mmu_, nnuf_, nnl_] = 
    (aapi^2*(11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
         Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
         Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
         Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/2 + 
     aapi^3*(((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        4 + ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((11*Log[mmu^2/nnuf^2])/2 - (nnl*Log[mmu^2/nnuf^2])/3)*Nm[nnl])/2) + 
     aapi^4*(((11 - (2*nnl)/3)^3*nnuf*Pi*((ctil[3, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[3 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[4 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        8 + ((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((33*Log[mmu^2/nnuf^2])/4 - (nnl*Log[mmu^2/nnuf^2])/2)*Nm[nnl])/4 + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         (612*Log[mmu^2/nnuf^2] - 76*nnl*Log[mmu^2/nnuf^2] + 
          1089*Log[mmu^2/nnuf^2]^2 - 132*nnl*Log[mmu^2/nnuf^2]^2 + 
          4*nnl^2*Log[mmu^2/nnuf^2]^2)*Nm[nnl])/96);

ctil[0,nf_] = 1;
ctil[1,3] = -0.1638; ctil[2,3] = 0.2372; ctil[3,3] = -0.1205; nu[3] = 0.3951;
ctil[1,4] = -0.1054; ctil[2,4] = 0.2736; ctil[3,4] = -0.1610; nu[4] = 0.3696;
ctil[1,5] =  0.0238; ctil[2,5] = 0.3265; ctil[3,5] = -0.2681; nu[5] = 0.3289;
Nm[3] = 0.563; (*\pm 26*)
Nm[4] = 0.547; (*\pm 33*)
Nm[5] = 0.527; (*\pm 51*)
 
exOS2RS[aapi_, mmu_, nnuf_, nnl_] = 
    aapi*nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*Nm[nnl] + 
     aapi^3*(((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        4 + nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
        Log[mmu^2/nnuf^2]*((102 - (38*nnl)/3)/16 + 
         ((11 - (2*nnl)/3)^2*Log[mmu^2/nnuf^2])/16)*Nm[nnl] + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((11*Log[mmu^2/nnuf^2])/2 - (nnl*Log[mmu^2/nnuf^2])/3)*Nm[nnl])/2) + 
     aapi^2*(((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        2 - (nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
         (-33*Log[mmu^2/nnuf^2] + 2*nnl*Log[mmu^2/nnuf^2])*Nm[nnl])/12) + 
     aapi^4*(((11 - (2*nnl)/3)^3*nnuf*Pi*((ctil[3, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[3 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[4 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        8 + ((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((33*Log[mmu^2/nnuf^2])/4 - (nnl*Log[mmu^2/nnuf^2])/2)*Nm[nnl])/4 + 
       nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
        Log[mmu^2/nnuf^2]*((2857/2 - (5033*nnl)/18 + (325*nnl^2)/54)/64 + 
         (5*(102 - (38*nnl)/3)*(11 - (2*nnl)/3)*Log[mmu^2/nnuf^2])/128 + 
         ((11 - (2*nnl)/3)^3*Log[mmu^2/nnuf^2]^2)/64)*Nm[nnl] + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         (612*Log[mmu^2/nnuf^2] - 76*nnl*Log[mmu^2/nnuf^2] + 
          1089*Log[mmu^2/nnuf^2]^2 - 132*nnl*Log[mmu^2/nnuf^2]^2 + 
          4*nnl^2*Log[mmu^2/nnuf^2]^2)*Nm[nnl])/96);

(* ************************************************************ *)
(*
mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the RS mass form the MS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", 1];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", 1, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", fdelm, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_String, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, PRIME, fdelm, "no"];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", 1];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", 1, "no"];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", fdelm, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_String, fdelm_, marker_String] := 
    Module[{extmp, api, delmuf, exmOS, exmOSprime, apiprime, exals}, 
     If[ marker == "no", mark = 1, mark = "RunDecXMS2OS" ];
     If[loops < 0 || loops > 4, Print["PROCEDURE IS NOT IMPLEMENTED FOR ", 
         loops, " LOOPS."]; Return[]]; If[loops == 0, extmp = mMS, 
       If[PRIME === "yes", delmuf = exOS2RSp[api, Mu, nuf, nl]; , 
         delmuf = exOS2RS[api, Mu, nuf, nl]; ]; 
        exmOS = mMS2mOS[mMS, {}, api*Pi, Mu, nl + 1, loops, fdelm]; 
	exals = (DecAsUpMS[apiprime*Pi, mMS, Mu, nl, 4])/(apiprime*Pi); 
        exmOSprime = Normal[Series[exmOS /. {api -> apiprime*exals}, 
           {apiprime, 0, loops}]]; extmp = 
         Expand[mark * (exmOSprime /. apiprime -> api) - delmuf] /. 
          cut[api, loops]; ]; Return[N[extmp /. api -> asmu/Pi, numprec]]; ];


(* ************************************************************ *)
(*
mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the MS mass form the RS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "no", 1];

mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "no", fdelm];

mRSp2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "yes", 1];

mRSp2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "yes", fdelm];
 
mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_] :=  Module[
    {mtmp, mMS}, 
    mtmp = mMS2mRS[mMS, mq, asnlmu, Mu, nuf, nl, loops, PRIME, fdelm]; 
    Return[mMS /. FindRoot[mRS == mtmp, {mMS, mRS}]]; ];


(* ************************************************************ *)
(*
mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the SI mass form the RS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "no", 1];

mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_ , fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "no", fdelm];

mRSp2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "yes", 1];

mRSp2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_ , fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "yes", fdelm];
 
mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_, PRIME_, fdelm_] := 
    Module[{mMS1, mMS, acc}, 
	   mMS1 = 0; mMS = mRS; 
	   acc = 1.*^-6; 
	   While[Abs[mMS1 - mMS] > acc, mMS1 = mMS; 
		 mMS = mRS2mMS[mRS, mq, asfct[mMS1], 
			       mMS1, nuf, nl, loops, PRIME, fdelm]; ]; 
	   Return[mMS]; ];


(* ************************************************************ *)
(*
 mOS2m1S[mOS_, mq_, asmu_, mu_, nl_, loops_] computes the 1S mass from the OS mass
*)

mOS2m1S[mOS_, mq_, asmu_, mu_, nl_, loops_] := 
    Module[{xx, EnC, E1p, nn, eps}, If[loops < 0 || loops > 4, 
       Print["mOS2m1S: PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; E1p = 
       EnC*(1 + (asmu*xx*(97/6 + nl*(-11/9 - (2*Log[(3*mu)/(4*asmu*mOS)])/
                 3) + 11*Log[(3*mu)/(4*asmu*mOS)]))/Pi + 
           (asmu^2*xx^2*(1793/12 + (2917/216 - (11*nl)/18 + nl^2/54)*Pi^2 - 
              (9*Pi^4)/32 + (927*Log[(3*mu)/(4*asmu*mOS)])/4 + 
              (363*Log[(3*mu)/(4*asmu*mOS)]^2)/4 + nl*(-1693/72 - 
                (193*Log[(3*mu)/(4*asmu*mOS)])/6 - 11*Log[(3*mu)/(4*asmu*
                     mOS)]^2 - (19*Zeta[3])/2) + nl^2*(77/108 + 
                Log[(3*mu)/(4*asmu*mOS)] + Log[(3*mu)/(4*asmu*mOS)]^2/3 + 
                (2*Zeta[3])/9) + (275*Zeta[3])/4))/Pi^2 + 
           (asmu^3*xx^3*(1267919/1728 + Pi^4*(-723119/51840 + (11*nl^2)/
                 1080 - nl^3/4860 + nl*(59677/77760 + (3*Log[(3*mu)/(4*asmu*
                       mOS)])/8) - (99*Log[(3*mu)/(4*asmu*mOS)])/16) + 
              (4521*Log[(3*mu)/(4*asmu*mOS)]^2)/2 + (1331*
                Log[(3*mu)/(4*asmu*mOS)]^3)/2 + (114917*Zeta[3])/48 + 
              Log[(3*mu)/(4*asmu*mOS)]*(247675/96 + (3025*Zeta[3])/2) + 
              Pi^2*(265.389067842508 + (865*Log[mu/mOS])/18 + 
                (26897*Log[(3*mu)/(4*asmu*mOS)])/108 + nl^2*(905/432 + 
                  (11*Log[(3*mu)/(4*asmu*mOS)])/9 - (11*Zeta[3])/6) + 
                nl^3*(-19/486 - (2*Log[(3*mu)/(4*asmu*mOS)])/81 + 
                  Zeta[3]/27) + nl*(-397591/7776 - (5095*Log[(3*mu)/(4*asmu*
                       mOS)])/162 + (121*Zeta[3])/4)) + (13432.614375 - 
                3289.906669391583*nl - (1000*nl^3)/729 + nl^2*
                 ((14002/81 - (416*Zeta[3])/3)/3 + (3*(12541/243 + (64*Pi^4)/
                      135 + (368*Zeta[3])/3))/4))/32 + nl*(-52033/288 - 
                (10955*Log[(3*mu)/(4*asmu*mOS)]^2)/24 - 
                121*Log[(3*mu)/(4*asmu*mOS)]^3 + Log[(3*mu)/(4*asmu*mOS)]*
                 (-166309/288 - (902*Zeta[3])/3) - (8797*Zeta[3])/18 - 
                363*Zeta[5]) + nl^3*(-98/729 - (5*Log[(3*mu)/(4*asmu*mOS)]^2)/
                 9 - (4*Log[(3*mu)/(4*asmu*mOS)]^3)/27 + 
                Log[(3*mu)/(4*asmu*mOS)]*(-50/81 - (8*Zeta[3])/27) - 
                (44*Zeta[3])/81 - (4*Zeta[5])/9) + (3993*Zeta[5])/2 + 
              nl^2*(3073/288 + (1027*Log[(3*mu)/(4*asmu*mOS)]^2)/36 + 
                (22*Log[(3*mu)/(4*asmu*mOS)]^3)/3 + (3239*Zeta[3])/108 + 
                Log[(3*mu)/(4*asmu*mOS)]*(10351/288 + (158*Zeta[3])/9) + 
                22*Zeta[5])))/Pi^3) /. {xx -> 1, 
          EnC -> -((asmu^2*cf^2*mOS)/(4*1^2))} /. {cf -> 4/3}; 
      E1p = Expand[(1*(E1p /. asmu -> eps*asmu))/eps] /. 
        {eps^(nn_.) /; nn > loops :> 0}; Return[mOS + (1*E1p)/2 /. 
        eps -> 1]; ];
 
(* ************************************************************ *)
(*
 mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_] 
 computes the 1S mass from the MS mass
*)

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, 1, "no"];

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, marker_String] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, 1, marker];

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_/;NumericQ[fdelm]] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, fdelm, "no"];
 
mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_, marker_String] := 
    Module[{mOS, m1S, res, epstmp, asmutmp, asmutmp2, keeplab, mark, markeps}, 
	   keeplab = "no";
     If[ marker == "no", mark = 1; markeps = 1, mark = "RunDecXMS2OS"; markeps = "RunDecXeps" ];
     If[loops == 0, Return[mMS]]; mOS = mMS2mOS[mMS, {}, asmutmp, mu, nl + 1, 
        loops, fdelm]; mOS = mOS /. {asmutmp -> DecAsUpMS[asmutmp2, mMS, mu, nl, 
           loops]}; mOS = Normal[Series[mOS, {asmutmp2, 0, loops}]] /. 
        {asmutmp2 -> asmu*eps*xmMS2mOS}; 
      m1S = Normal[Series[mOS2m1S[mOS, {}, asmu*epstmp, mu, nl, loops], 
         {epstmp, 0, loops + 3}]]; 
      m1S = Normal[Series[m1S /. {Log[epstmp] -> 0} /. 
          {epstmp^(nn_) :> eps^(nn - 1)*xmOS2m1S}, {eps, 0, loops}]]; 
      If[keeplab === "yes", 
       m1S = Expand[PowerExpand[m1S] /. Log[eps] -> 0 /. Log[xmOS2m1S] -> 
                0 /. Log[xmMS2mOS] -> 0 /. Log[xx_] :> Log[xx /. 
                {eps -> 1, xmMS2mOS -> 1}]] //. {xmOS2m1S^(nn_) :> xmOS2m1S, 
            xmMS2mOS^(nn_) :> xmMS2mOS} /. {xmMS2mOS*xmOS2m1S -> xmOS2m1S}; 
        res = Collect[{mOS, m1S}, {eps}, Expand]; Print["mOS = ", 
         InputForm[res[[1]]]]; m1S = res[[2]]; Print["m1S = ", 
         InputForm[m1S]]; 
         , 
         m1S = m1S /. {eps -> markeps, xmOS2m1S -> 1, xmMS2mOS -> mark}; ]; 
         Return[m1S]; ];


(* ************************************************************ *)
(*
 m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_, fdelm_]
 computes the MS mass from the 1S mass
*)

m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_] := m1S2mMS[m1S, mq, asnlmu, Mu, 
     nl, loops, 1];
 
m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_, fdelm_] := 
    Module[{mtmp, mMS}, mtmp = mMS2m1S[mMS, mq, asnlmu, Mu, nl, loops, fdelm];
      Return[mMS /. FindRoot[m1S == mtmp, {mMS, m1S}]]; ];


(* ************************************************************ *)
(*
 m1S2mSI[m1S_, mq_, asfct_, nl_, loops_,fdelm_]
 computes the SI mass from the 1S mass
*)

m1S2mSI[m1S_, mq_, asfct_, nl_, loops_] := m1S2mSI[m1S, mq, asfct, nl, loops, 1];
 
m1S2mSI[m1S_, mq_, asfct_, nl_, loops_, fdelm_] := 
    Module[{mMS1, mMS, acc}, 
	   mMS1 = 0; mMS = m1S; 
	   acc = 1.*^-6; 
	   While[Abs[mMS1 - mMS] > acc, 
		 mMS1 = mMS; mMS = m1S2mMS[m1S, mq, asfct[mMS1], mMS1, nl, loops, 
					   fdelm]; ]; 
	   Return[mMS]; ];

(* ************************************************************ *)

End[];

(* ************************************************************ *)

EndPackage[];

(* ************************************************************ *)
(* ************************************************************ *)

