#include <TCanvas.h>
#include <TGraph.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TVectorD.h>
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
double lam_ckm = .2257, lam_ckm_err = .0010;
double A_ckm = .814, A_ckm_err = .022;
double rho_ckm = .135, rho_ckm_err = .031;
double eta_ckm = .349, eta_ckm_err = .017;

class RSFitter
{
	int num_rolls;
	std::string id_string;
	double lambda;
	int index;
	std::vector<double> masses;
	std::vector<double> masses_err;
	double chi_squared;
	double cq1,cq2,cq3;
	double cu1,cu2,cu3;
	double cd1,cd2,cd3;
	std::complex<double> Yu11,Yu12,Yu13,Yu21,Yu22,Yu23,Yu31,Yu32,Yu33;
	std::complex<double> Yd11,Yd12,Yd13,Yd21,Yd22,Yd23,Yd31,Yd32,Yd33;
	bool penalty_on;

	public:
		RSFitter(int idx, int rolls, double lam, std::vector<double> m, std::vector<double> m_err);


	std::complex<double> div(cp a, cp b);
	double penalty(std::complex<double> x, double low_bound,double up_bound);
	void turn_on_penalty(bool b) { penalty_on = b; };
	double m2y(double m);
	double y2m(double m);
	double kr(double lam);
	double yuk(double cl, double cr);
	double Chisq(const double *par);
	TTree *fit(int idx, double lam, std::vector<double> m, std::vector<double> m_err);

};

RSFitter::RSFitter(int idx, int rolls, double lam, std::vector<double> m, std::vector<double> m_err) {
	num_rolls = rolls;
	lambda = lam;
	index = idx;
	masses = m;
	masses_err = m_err;
	id_string = std::to_string(int(index));
	penalty_on = false;
	chi_squared = 99999999;
}

std::complex<double> RSFitter::div(cp a, cp b) {
	double real = (a.real()*b.real() + a.imag()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));
	double imaginary = (a.imag()*b.real() - a.real()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));

	return cp(real, imaginary);
}

double RSFitter::penalty(std::complex<double> x, double low_bound, double up_bound) {
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

//Quark yukawa couplings as a function of 5d parameters c
double RSFitter::yuk(double cl, double cr) {
	return exp((1-cl-cr)*M_PI*kr(lambda))*sqrt((.5 - cl)/(exp((1 - 2*cl)*M_PI*kr(lambda)) - 1 ))*sqrt((.5 - cr)/(exp((1 - 2*cr)*M_PI*kr(lambda)) - 1 ));
}

double RSFitter::Chisq(const double *par) {

	//5d c parameters
	Vector3d cq(par[0],par[1],par[2]);
	Vector3d cu(par[3],par[4],par[5]);
	Vector3d cd(par[6],par[7],par[8]);

	//5d yukawa matrices
	Matrix3cd yu;
	yu << cp(par[9],par[27]), cp(par[10],par[28]), cp(par[11],par[29]),
		   cp(par[12],par[30]), cp(par[13],par[31]), cp(par[14],par[32]),
		   cp(par[15],par[33]), cp(par[16],par[34]), cp(par[17],par[35]);
	Matrix3cd yd;
	yd << cp(par[18],par[36]), cp(par[19],par[37]), cp(par[20],par[38]),
		   cp(par[21],par[39]), cp(par[22],par[40]), cp(par[23],par[41]),
		   cp(par[24],par[42]), cp(par[25],par[43]), cp(par[26],par[44]);

	Matrix3cd cu5;
	Matrix3cd cd5;
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			cu5(i,j) = yuk(cq(i),cu(j));
			cd5(i,j) = yuk(cq(i),cd(j));
		}
	}

	yu = yu.cwiseProduct(cu5);
	yd = yd.cwiseProduct(cd5);

	JacobiSVD<Matrix3cd> u_svd(yu,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3cd> d_svd(yd,ComputeFullU | ComputeFullV);

	Matrix3cd Vp = u_svd.matrixU().adjoint()*d_svd.matrixU();
	double y_pred[6] = {u_svd.singularValues()[2], d_svd.singularValues()[2], d_svd.singularValues()[1], u_svd.singularValues()[1], d_svd.singularValues()[0], u_svd.singularValues()[0]};

	double chisq = 0;
	for(int i = 0; i < 6; i++) {
		chisq += pow((y_pred[i] - m2y(masses[i]))/m2y(masses_err[i]), 2);
	}

	double lam_ckm_pred = abs(Vp(0,1))/sqrt(pow(abs(Vp(0,0)),2) + pow(abs(Vp(0,1)),2));
	double A_ckm_pred = abs(div(Vp(1,2),Vp(0,1)))/lam_ckm_pred;
	double rho_ckm_pred = Vp(0,2).real()/(A_ckm_pred*pow(lam_ckm_pred,3));
	double eta_ckm_pred = -Vp(0,2).imag()/(A_ckm_pred*pow(lam_ckm_pred,3));

	chisq += pow((lam_ckm_pred - lam_ckm)/lam_ckm_err,2) + pow((A_ckm_pred - A_ckm)/A_ckm_err,2) + pow((rho_ckm_pred - rho_ckm)/rho_ckm_err,2) + pow((eta_ckm_pred - eta_ckm)/eta_ckm_err,2);

	//Penalize .08 < Abs(Y) < 4.0
	if (penalty_on == true) {
		for(int i = 9; i < 27; i++) {
			chisq += penalty(cp(par[i],par[i+18]),.01,4);
		}
	}

	return chisq;

}

TTree *RSFitter::fit(int idx, double lam, std::vector<double> m, std::vector<double> m_err) {

	id_string = std::to_string(int(idx)); lambda = lam;
	masses = m; masses_err = m_err;

	cout << "Index = " << id_string << "\n";
	cout << "Lambda = " << lambda << "\n";
	cout << "Masses[0] = " << masses[0] << "\n";
	cout << "Masses_Err[0] = " << masses_err[0] << "\n";

	TTree *tree = new TTree(id_string.c_str(),id_string.c_str());

	tree->Branch("chi_squared",&chi_squared);
	tree->Branch("cq1",&cq1);
	tree->Branch("cq2",&cq2);
	tree->Branch("cq3",&cq3);
	tree->Branch("cu1",&cu1);
	tree->Branch("cu2",&cu2);
	tree->Branch("cu3",&cu3);
	tree->Branch("cd1",&cd1);
	tree->Branch("cd2",&cd2);
	tree->Branch("cd3",&cd3);
	tree->Branch("Yu11", &Yu11); tree->Branch("Yu12", &Yu12); tree->Branch("Yu13", &Yu13);
	tree->Branch("Yu21", &Yu21); tree->Branch("Yu22", &Yu22); tree->Branch("Yu23", &Yu23);
	tree->Branch("Yu31", &Yu31); tree->Branch("Yu32", &Yu32); tree->Branch("Yu33", &Yu33);
	tree->Branch("Yd11", &Yd11); tree->Branch("Yd12", &Yd12); tree->Branch("Yd13", &Yd13);
	tree->Branch("Yd21", &Yd21); tree->Branch("Yd22", &Yd22); tree->Branch("Yd23", &Yd23);
	tree->Branch("Yd31", &Yd31); tree->Branch("Yd32", &Yd32); tree->Branch("Yd33", &Yd33);

	ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
	minimum->SetMaxIterations(10000);  // for GSL
	minimum->SetTolerance(0.01);
	minimum->SetPrintLevel(0);
	ROOT::Math::Functor func(this,&RSFitter::Chisq,45);
	minimum->SetFunction(func);

	int i = 0;
	while (i < num_rolls) {
		minimum->Clear();

		rng.seed(rd());
		std::uniform_real_distribution<> cs(-1.0,2.0);
	  std::uniform_real_distribution<> lams(-4.0,4.0);

		minimum->SetLimitedVariable(0,"cq1",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(1,"cq2",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(2,"cq3",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(3,"cu1",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(4,"cu2",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(5,"cu3",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(6,"cd1",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(7,"cd2",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(8,"cd3",cs(rng),.01,-1.0,2.0);
		minimum->SetLimitedVariable(9,"yu1",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(10,"yu2",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(11,"yu3",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(12,"yu4",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(13,"yu5",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(14,"yu6",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(15,"yu7",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(16,"yu8",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(17,"yu9",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(18,"yd1",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(19,"yd2",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(20,"yd3",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(21,"yd4",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(22,"yd5",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(23,"yd6",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(24,"yd7",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(25,"yd8",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(26,"yd9",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(27,"yu1i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(28,"yu2i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(29,"yu3i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(30,"yu4i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(31,"yu5i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(32,"yu6i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(33,"yu7i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(34,"yu8i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(35,"yu9i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(36,"yd1i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(37,"yd2i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(38,"yd3i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(39,"yd4i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(40,"yd5i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(41,"yd6i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(42,"yd7i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(43,"yd8i",lams(rng),.01,-4.0,4.0);
		minimum->SetLimitedVariable(44,"yd9i",lams(rng),.01,-4.0,4.0);

		minimum->Minimize();
		const double *x = minimum->X();

		chi_squared = minimum->MinValue();

		cq1 = x[0]; cq2 = x[1]; cq3 = x[2];
		cu1 = x[3]; cu2 = x[4]; cu3 = x[5];
		cd1 = x[6]; cd2 = x[7]; cd3 = x[8];

		Yu11 = cp(x[9],x[27]); Yu12 = cp(x[10],x[28]); Yu13 = cp(x[11],x[29]);
		Yu21 = cp(x[12],x[30]); Yu22 = cp(x[13],x[31]); Yu23 = cp(x[14],x[32]);
		Yu31 = cp(x[15],x[33]); Yu32 = cp(x[16],x[34]); Yu33 = cp(x[17],x[35]);

		Yd11 = cp(x[18],x[36]); Yd12 = cp(x[19],x[37]); Yd13 = cp(x[20],x[38]);
		Yd21 = cp(x[21],x[39]); Yd22 = cp(x[22],x[40]); Yd23 = cp(x[23],x[41]);
		Yd31 = cp(x[24],x[42]); Yd32 = cp(x[25],x[43]); Yd33 = cp(x[26],x[44]);

		//Accept fit if chi_squared < 10
		if (chi_squared < 10.0) {
			tree->Fill();
		}
		i++;
	}
	return tree;
}












