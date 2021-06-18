#include <TCanvas.h>
#include <TGraph.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <iostream>
#include <stdio.h>
#include <random>
#include <math.h>
#include <complex>
#include </Users/noahsteinberg/Physics/eigen/Eigen/SVD>
#include </Users/noahsteinberg/Physics/eigen/Eigen/Dense>

typedef std::complex<double> cp;
typedef std::mt19937 MyRNG;
std::random_device rd;

MyRNG rng;

using namespace Eigen;
using std::cout;
using std::abs;

//Initialize some useful parameters
double M_pl = 2.435e18;
double vev = 246;
//double lambda = 162.9;

class RSFitter
{
	int num_rolls;
	double lambda;
	int index;

	//All values given in GeV
	double l_dat[3] ={.0004853,.01025, 1.742};
	double l_dat_err[3] = {7.3e-6,.00015,.026};//arbitrarily gave 1.5% errors (perscribed in Iyer & Vempati)
	double nu_dat[2] ={7.53e-23,2.51e-21};//Dmsq_21, Dmsq_31
	double nu_dat_err[2] = {.18e-23,.05e-21};
	double nu1_mass = 1e-11;
	double s12 = .307, s12_err = .013;
	double s23 = .417, s23_err = .028;
	double s13 = .0212, s13_err = .0008;
	double jarl = -.03, jarl_err = .01;

	std::vector<double> masses;
	std::vector<double> masses_err;
	double chi_squared;
	double cl1,cl2,cl3;
	double ce1,ce2,ce3;
	double cn1,cn2,cn3;
	std::complex<double> Yl11,Yl12,Yl13,Yl21,Yl22,Yl23,Yl31,Yl32,Yl33;
	std::complex<double> Ynu11,Ynu12,Ynu13,Ynu22,Ynu23,Ynu33;
	std::complex<double> Ynu21,Ynu31,Ynu32;
	TTree *tree;
	bool penalty_on;

	public:
		RSFitter(int idx, int rolls, double lam, std::vector<double> m, std::vector<double> m_err);

	std::complex<double> div(cp a, cp b);
	double penalty(std::complex<double> x, double low_bound,double up_bound);
	void turn_on_penalty(bool b) { penalty_on = b; };
	double m2y(double m);
	double y2m(double m);
	double kr(double lam);
	double lmass(double cl, double cr);
	double numass(double ci,double cj);
	double Chisq(const double *par);
	TTree *fit();

};

RSFitter::RSFitter(int idx, int rolls, double lam, std::vector<double> m, std::vector<double> m_err) {
	index = idx;
	num_rolls = rolls;
	lambda = lam;
	masses = m;
	masses_err = m_err;
	std::string lam_string = std::to_string(int(index));
	penalty_on = false;

	tree = new TTree(lam_string.c_str(),lam_string.c_str());
	tree->Branch("Lambda",&lambda);
	tree->Branch("chi_squared",&chi_squared);
	tree->Branch("cl1",&cl1);
	tree->Branch("cl2",&cl2);
	tree->Branch("cl3",&cl3);
	tree->Branch("ce1",&ce1);
	tree->Branch("ce2",&ce2);
	tree->Branch("ce3",&ce3);
	tree->Branch("cn1",&cn1);
	tree->Branch("cn2",&cn2);
	tree->Branch("cn3",&cn3);

	tree->Branch("Yl11", &Yl11); tree->Branch("Yl12", &Yl12); tree->Branch("Yl13", &Yl13);
	tree->Branch("Yl21", &Yl21); tree->Branch("Yl22", &Yl22); tree->Branch("Yl23", &Yl23);
	tree->Branch("Yl31", &Yl31); tree->Branch("Yl32", &Yl32); tree->Branch("Yl33", &Yl33);
	tree->Branch("Ynu11", &Ynu11); tree->Branch("Ynu12", &Ynu12); tree->Branch("Ynu13", &Ynu13);
	tree->Branch("Ynu21", &Ynu21); tree->Branch("Ynu22", &Ynu22); tree->Branch("Ynu23", &Ynu23);
	tree->Branch("Ynu31", &Ynu31); tree->Branch("Ynu32", &Ynu32); tree->Branch("Ynu33", &Ynu33);
}

std::complex<double> RSFitter::div(cp a, cp b) {
	double real = (a.real()*b.real() + a.imag()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));
	double imaginary = (a.imag()*b.real() - a.real()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));

	return cp(real, imaginary);
}

double RSFitter::penalty(std::complex<double> x, double low_bound,double up_bound) {
	double mag = abs(x);
	if (mag >= low_bound && mag <= up_bound) return 0.0;
	else {
		return 1.0e20;
	}
}

//Function that converts masses to yukawa couplings
double RSFitter::m2y(double m) {
	return sqrt(2)*m/vev;
}

double RSFitter::y2m(double y) {
	return y*vev/sqrt(2);
}

//KR as a function of the IR scale
double RSFitter::kr(double lam) {
	return log(M_pl/lam)/M_PI;
}

double RSFitter::lmass(double cl, double cr) {
	return (vev/sqrt(2))*exp((1-cl-cr)*M_PI*kr(lambda))*sqrt((.5 - cl)/(exp((1 - 2*cl)*M_PI*kr(lambda)) - 1 ))*sqrt((.5 - cr)/(exp((1 - 2*cr)*M_PI*kr(lambda)) - 1 ));
}

double RSFitter::numass(double ci, double cj) {
	return (vev*vev/M_pl)*exp(M_PI*kr(lambda)*(2 - ci - cj))*sqrt((.5 - ci)/(exp((1 - 2*ci)*M_PI*kr(lambda)) - 1 ))*sqrt((.5 - cj)/(exp((1 - 2*cj)*M_PI*kr(lambda)) - 1 ));
}

double RSFitter::Chisq(const double *par) {

	//5d c parameters
	Vector3d cl(par[0],par[1],par[2]);
	Vector3d ce(par[3],par[4],par[5]);
	//Added in for dirac neutrinos
	Vector3d cn(par[6],par[7],par[8]);

	/* LLHH operator
	//5d yukawa matrices
	Matrix3cd yl;
	yl << cp(par[6],par[21]), cp(par[7],par[22]), cp(par[8],par[23]),
		   cp(par[9],par[24]), cp(par[10],par[25]), cp(par[11],par[26]),
		   cp(par[12],par[27]), cp(par[13],par[28]), cp(par[14],par[29]);
	Matrix3cd ynu;
	ynu << cp(par[15],par[30]), cp(par[16],par[31]), cp(par[17],par[32]),
		   cp(par[16],par[31]), cp(par[18],par[33]), cp(par[19],par[34]),
		   cp(par[17],par[32]), cp(par[19],par[34]), cp(par[20],par[35]);
	*/


	// Dirac neutrinos
	Matrix3cd yl;
	yl << cp(par[9],par[27]), cp(par[10],par[28]), cp(par[11],par[29]),
		   cp(par[12],par[30]), cp(par[13],par[31]), cp(par[14],par[32]),
		   cp(par[15],par[33]), cp(par[16],par[34]), cp(par[17],par[35]);
	Matrix3cd ynu;
	ynu << cp(par[18],par[36]), cp(par[19],par[37]), cp(par[20],par[38]),
		   cp(par[21],par[39]), cp(par[22],par[40]), cp(par[23],par[41]),
		   cp(par[24],par[42]), cp(par[25],par[43]), cp(par[26],par[44]);

	Matrix3cd cl5;
	Matrix3cd cnu5;
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			cl5(i,j) = lmass(cl(i),ce(j));
			//cnu5(i,j) = numass(cl(i),cl(j));
			cnu5(i,j) = lmass(cl(i),cn(j));
		}
	}

	yl = yl.cwiseProduct(cl5);
	ynu = ynu.cwiseProduct(cnu5);

	JacobiSVD<Matrix3cd> l_svd(yl,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3cd> nu_svd(ynu,ComputeFullU | ComputeFullV);

	Matrix3cd Vp = nu_svd.matrixU().adjoint()*l_svd.matrixU();
	double yl_pred[3] = {l_svd.singularValues()[2],l_svd.singularValues()[1],l_svd.singularValues()[0]};
	double ynu_pred[3] = {nu_svd.singularValues()[2],nu_svd.singularValues()[1],nu_svd.singularValues()[0]};

	double m21_pred = ynu_pred[1]*ynu_pred[1] - ynu_pred[2]*ynu_pred[2];
	double m31_pred = ynu_pred[0]*ynu_pred[0] - ynu_pred[2]*ynu_pred[2];

	double chisq = 0;
	for(int i = 0; i < 3; i++) {
		chisq += pow((yl_pred[i] - l_dat[i])/l_dat_err[i],2);
	}

	chisq += pow((m21_pred - nu_dat[0])/nu_dat_err[0],2) + pow((m31_pred - nu_dat[1])/nu_dat_err[1],2) + pow((ynu_pred[2] - nu1_mass),2);

	double s12_pred = pow(abs(Vp(0,1)),2)/(1 - pow(abs(Vp(0,2)),2));
	double s13_pred = pow(abs(Vp(0,2)),2);
	double s23_pred = pow(abs(Vp(1,2)),2)/(1 - pow(abs(Vp(0,2)),2));
	double jarl_pred = (Vp(1,2)*conj(Vp(0,2))*Vp(0,1)*conj(Vp(1,1))).imag();

	//Leaving out CP violation right now

	chisq += pow((s12_pred - s12)/s12_err,2) + pow((s13_pred - s13)/s13_err,2) + pow((s23_pred - s23)/s23_err,2) /*+ pow((jarl_pred - jarl)/jarl_err,2)*/;

	//Penalize .08 < Abs(Y) < 4.0
	if (penalty_on == true) {
		//Dirac
		for(int i = 9; i < 27; i++) {
		//LLHH
		//for(int i = 6; i < 21; i++) { //LLHH
			//chisq += penalty(cp(par[i],par[i+15]),.01,4); // LLHH
			chisq += penalty(cp(par[i],par[i+18]),.01,4); //Dirac
		}
	}

	return chisq;

}

TTree *RSFitter::fit() {

	ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
	minimum->SetMaxIterations(100000);  // for GSL
	minimum->SetTolerance(0.01);
	minimum->SetPrintLevel(0);
	//LLHH
	//ROOT::Math::Functor func(this,&RSFitter::Chisq,36);
	//Dirac
	ROOT::Math::Functor func(this,&RSFitter::Chisq,45);
	minimum->SetFunction(func);

	int i = 0;
	int good_sols = 0;
	while (i < num_rolls) {
		minimum->Clear();

		rng.seed(rd());
		//std::uniform_real_distribution<> cls(-1,1);
		//std::uniform_real_distribution<> ces(-1e8,0);
		std::uniform_real_distribution<> cs(-1.0,2.0);
	  std::uniform_real_distribution<> lams(-4.0,4.0);

	  //LLHH case
	  /*
		minimum->SetVariable(0,"cl1",cls(rng),.01);
		minimum->SetVariable(1,"cl2",cls(rng),.01);
		minimum->SetVariable(2,"cl3",cls(rng),.01);
		minimum->SetVariable(3,"ce1",cls(rng),.01);
		minimum->SetVariable(4,"ce2",cls(rng),.01);
		minimum->SetVariable(5,"ce3",cls(rng),.01);
		minimum->SetLimitedVariable(6,"yl1",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(7,"yl2",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(8,"yl3",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(9,"yl4",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(10,"yl5",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(11,"yl6",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(12,"yl7",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(13,"yl8",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(14,"yl9",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(15,"ynu1",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(16,"ynu2",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(17,"ynu3",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(18,"ynu5",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(19,"ynu6",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(20,"ynu9",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(21,"yl1i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(22,"yl2i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(23,"yl3i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(24,"yl4i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(25,"yl5i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(26,"yl6i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(27,"yl7i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(28,"yl8i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(29,"yl9i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(30,"ynu1i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(31,"ynu2i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(32,"ynu3i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(33,"ynu5i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(34,"ynu6i",lams(rng),.1,-4.0,4.0);
		minimum->SetLimitedVariable(35,"ynu9i",lams(rng),.1,-4.0,4.0);
		*/

	  //Dirac mass case
	  minimum->SetVariable(0,"cl1",cs(rng),.01);
		minimum->SetVariable(1,"cl2",cs(rng),.01);
		minimum->SetVariable(2,"cl3",cs(rng),.01);
		minimum->SetVariable(3,"ce1",cs(rng),.01);
		minimum->SetVariable(4,"ce2",cs(rng),.01);
		minimum->SetVariable(5,"ce3",cs(rng),.01);
		minimum->SetVariable(6,"cn1",cs(rng),.01);
		minimum->SetVariable(7,"cn2",cs(rng),.01);
		minimum->SetVariable(8,"cn3",cs(rng),.01);
		minimum->SetLimitedVariable(9,"yl1",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(10,"yl2",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(11,"yl3",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(12,"yl4",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(13,"yl5",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(14,"yl6",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(15,"yl7",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(16,"yl8",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(17,"yl9",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(18,"ynu1",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(19,"ynu2",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(20,"ynu3",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(21,"ynu4",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(22,"ynu5",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(23,"ynu6",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(24,"ynu7",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(25,"ynu8",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(26,"ynu9",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(27,"yl1i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(28,"yl2i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(29,"yl3i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(30,"yl4i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(31,"yl5i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(32,"yl6i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(33,"yl7i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(34,"yl8i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(35,"yl9i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(36,"ynu1i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(37,"ynu2i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(38,"ynu3i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(39,"ynu4i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(40,"ynu5i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(41,"ynu6i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(42,"ynu7i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(43,"ynu8i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(44,"ynu9i",lams(rng),.01,-4.0,4.0);


		minimum->Minimize();
		const double *x = minimum->X();

		chi_squared = minimum->MinValue();

		/*
		cl1 = x[0]; cl2 = x[1]; cl3 = x[2];
		ce1 = x[3]; ce2 = x[4]; ce3 = x[5];

		Yl11 = cp(x[6],x[21]); Yl12 = cp(x[7],x[22]); Yl13 = cp(x[8],x[23]);
		Yl21 = cp(x[9],x[24]); Yl22 = cp(x[10],x[25]); Yl23 = cp(x[11],x[26]);
		Yl31 = cp(x[12],x[27]); Yl32 = cp(x[13],x[28]); Yl33 = cp(x[14],x[29]);

		Ynu11 = cp(x[15],x[30]); Ynu12 = cp(x[16],x[31]); Ynu13 = cp(x[17],x[32]);
		Ynu22 = cp(x[18],x[33]); Ynu23 = cp(x[19],x[34]);
		Ynu33 = cp(x[20],x[35]);
		*/

		cl1 = x[0]; cl2 = x[1]; cl3 = x[2];
		ce1 = x[3]; ce2 = x[4]; ce3 = x[5];
		cn1 = x[6]; cn2 = x[7]; cn3 = x[8];

		Yl11 = cp(x[9],x[27]); Yl12 = cp(x[10],x[28]); Yl13 = cp(x[11],x[29]);
		Yl21 = cp(x[12],x[30]); Yl22 = cp(x[13],x[31]); Yl23 = cp(x[14],x[32]);
		Yl31 = cp(x[15],x[33]); Yl32 = cp(x[16],x[34]); Yl33 = cp(x[17],x[35]);

		Ynu11 = cp(x[18],x[36]); Ynu12 = cp(x[19],x[37]); Ynu13 = cp(x[20],x[38]);
		Ynu21 = cp(x[21],x[39]); Ynu22 = cp(x[22],x[40]); Ynu23 = cp(x[23],x[41]);
		Ynu31 = cp(x[24],x[42]); Ynu32 = cp(x[25],x[43]); Ynu33 = cp(x[26],x[44]);

		cout << "Chisquared = " << chi_squared << "\n";

		//Accept fit if chi_squared < 10
		if (chi_squared < 10.0) {
			tree->Fill();
			good_sols++;
		}
		i++;
	}

	cout << "Number of good solutions = " << good_sols << "\n";
	return tree;
}












