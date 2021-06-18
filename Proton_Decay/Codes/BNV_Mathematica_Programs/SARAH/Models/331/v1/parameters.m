ParameterDefinitions = { 

{g1, {   Description -> "X-Coupling", 
	 LaTeX -> "g_1", 
	 LesHouches -> {gauge,1},
	 OutputName-> g1 }},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    { Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}}, 

{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{mS,   {LaTeX -> "m_S",
	LesHouches -> MS,
        OutputName->mS}},
{yl,   {LaTeX -> "y_\\ell",
	LesHouches -> YL,
	OutputName->yl }},
{ys,   {LaTeX -> "y_s",
	LesHouches -> YS,
	OutputName->ys }},
{ya,   {LaTeX -> "y_a",
	LesHouches -> YA,
	OutputName->ya }},
{yd1,  {LaTeX -> "y^d_1",
	LesHouches -> YD1,
	OutputName->yd1 }},
{yd2,  {LaTeX -> "y^d_2",
	LesHouches -> YD2,
	OutputName->yd2 }},
{yd3,  {LaTeX -> "\\tilde{y}^d",
	LesHouches -> YD3,
	OutputName->yd3 }},
{yd1t, {LaTeX -> "\\bar{y}^d_1",
	LesHouches -> YD1T,
	OutputName->yd1t }},
{yd2t, {LaTeX -> "\\bar{y}^d_2",
	LesHouches -> YD2T,
	OutputName->yd2t }},
{yd1Xt, {LaTeX -> "\\bar{y}^d_{X 1}",
	 LesHouches -> YD1XT,
	 OutputName->yd1Xt }},
{yd2Xt, {LaTeX -> "\\bar{y}^d_{X 1}",
	 LesHouches -> YD2XT,
	 OutputName->yd2Xt }},
{yu1,   {LaTeX -> "y^u_1",
	 LesHouches -> YU1,
	 OutputName->yu1 }},
{yu2,   {LaTeX -> "y^u_2",
	 LesHouches -> YU2,
	 OutputName->yu2 }},
{yu3,   {LaTeX -> "\\tilde{y}^u",
	 LesHouches -> YU3,
	 OutputName->yu3 }},
{yu3t,  {LaTeX -> "\\bar{y}^u",
	 LesHouches -> YU3T,
	 OutputName->yu3t }},
{yu3Xt, {LaTeX -> "\\bar{y}^u_X",
	 LesHouches -> YU3XT,
	 OutputName->yu3Xt }},
{mu12,  {LaTeX -> "\\mu_1^2",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,1},
	 OutputName->mu12 }},
{mu22,  {LaTeX -> "\\mu_2^2",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,2},
	 OutputName->mu22 }},
{mu32,  {LaTeX -> "\\mu_3^2",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,3},
	 OutputName->mu32 }},
{muX2,  {LaTeX -> "\\mu_X^2",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,4},
	 OutputName->muX2 }},
{f,     {LaTeX -> "f",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,5},
	 OutputName->fcoup }},
{kap,   {LaTeX -> "\\kappa",
	 Real -> True,
	 LesHouches -> {POTENTIAL331,6},
	 OutputName->kap }},
{l1,    {LaTeX -> "\\lambda_1",
	 LesHouches -> {POTENTIAL331,7},
	 OutputName-> l1 }},
{l2,    {LaTeX -> "\\lambda_2",
	 LesHouches -> {POTENTIAL331,8},
	 OutputName-> l2 }},
{l3,    {LaTeX -> "\\lambda_3",
	 LesHouches -> {POTENTIAL331,9},
	 OutputName-> l3 }},
{lX,    {LaTeX -> "\\lambda_X",
	 LesHouches -> {POTENTIAL331,10},
	 OutputName-> lX }},
{l12,   {LaTeX -> "\\lambda_{12}",
	 LesHouches -> {POTENTIAL331,11},
	 OutputName-> l12 }},
{l13,   {LaTeX -> "\\lambda_{13}",
	 LesHouches -> {POTENTIAL331,12},
	 OutputName-> l13 }},
{l23,   {LaTeX -> "\\lambda_{23}",
	 LesHouches -> {POTENTIAL331,13},
	 OutputName-> l23 }},
{l1X,   {LaTeX -> "\\lambda_{1X}",
	 LesHouches -> {POTENTIAL331,14},
	 OutputName-> l1X }},
{l2X,   {LaTeX -> "\\lambda_{2X}",
	 LesHouches -> {POTENTIAL331,15},
	 OutputName-> l2X }},
{l3X,   {LaTeX -> "\\lambda_{3X}",
	 LesHouches -> {POTENTIAL331,16},
	 OutputName-> l3X }},
{vevk1, {LaTeX -> "k_1",
	 LesHouches -> {POTENTIAL331,17},
	 OutputName-> vevk1 }},
{vevk3, {LaTeX -> "k_3",
	 LesHouches -> {POTENTIAL331,18},
	 OutputName-> vevk3 }},  
{vevn,  {LaTeX -> "n",
	 LesHouches -> {POTENTIAL331,19},
	 OutputName-> vevn }},                          
                                                                           

{mH2,       { Description -> "SM Higgs Mass Parameter"}},

{ThetaW,    { Description -> "Weinberg-Angle",
              DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]}},

{ZW,{ 
     Description -> "Charged gauge boson Mixing Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZWMIX, 
     LaTeX -> "Z^{W}", 
     OutputName -> ZW}}, 
{ZZ,{ 
     Description -> "Neutral gauge boson Mixing Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZZMIX, 
     LaTeX -> "Z^{Z}", 
     OutputName -> ZZ}}, 
{ZH,{ 
     Description -> "Scalar-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZHMIX, 
     LaTeX -> "Z^{H}", 
     OutputName -> ZH}},
{ZHd,{ 
     Description -> "Scalar-Mixing-Matrix 2", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZHdMIX, 
     LaTeX -> "Z^{H}_2", 
     OutputName -> ZHd}},
{ZA,{ 
     Description -> "Pseudo-Scalar-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZAMIX, 
     LaTeX -> "Z^{A}", 
     OutputName -> ZA}},
{ZAd,{ 
     Description -> "Pseudo-Scalar-Mixing-Matrix 2", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZAdMIX, 
     LaTeX -> "Z^{A}_2", 
     OutputName -> ZAd}},
{ZP,{ 
     Description -> "Charged-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZPMIX, 
     LaTeX -> "Z^{\\pm}", 
     OutputName -> ZP}},
{ZPd,{ 
     Description -> "Charged-Mixing-Matrix 2", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> ZPdMIX, 
     LaTeX -> "Z^{\\pm}_2", 
     OutputName -> ZPd}},
{UV, {Description -> "Neutral fermion Mixing Matrix",
      LaTeX -> "U_\\nu",
      LesHouches -> UVMIX,
      OutputName -> UV }},

{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}}

 }; 
 

