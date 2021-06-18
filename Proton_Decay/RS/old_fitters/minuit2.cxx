#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TString.h"
#include <iostream>
#include <stdio.h>
#include <random>
#include <math.h>
#include </Users/noahsteinberg/Physics/eigen/Eigen/SVD>
#include </Users/noahsteinberg/Physics/eigen/Eigen/Dense>


typedef std::mt19937 MyRNG;
std::random_device rd;

MyRNG rng;

using namespace Eigen;
using std::cout;

//Initialize some useful parameters
double M_pl = 2.435e18;
double vev = 246;
double lambda = 162.9;

double lam_ckm = .2257, lam_ckm_err = .0010;
double A_ckm = .814, A_ckm_err = .022;
double rhoeta_ckm = .349, rhoeta_ckm_err = .017;

//Function that converts masses to yukawa couplings
double m2y(double m) {
	return sqrt(2)*m/vev;
}

double y2m(double y) {
	return y*vev/sqrt(2);
}

//Initialize quark masses and their errors at mtop(mtop)
double yu_dat[3] ={m2y(.00122),m2y(.59), m2y(162.9)};
double yu_dat_err[3] = {m2y(.00045),m2y(.080),m2y(2.8)};
//double yu_dat_err[3] = {1,1,1};
double yd_dat[3] ={m2y(.00276),m2y(.052),m2y(2.75)};
double yd_dat_err[3] = {m2y(.00117),m2y(.015),m2y(.09)};
//double yd_dat_err[3] = {1,1,1};

//CKM matrix
double v00 = .97427, v01 = .22536, v02 = .00355;
double v10 = .22522, v11 = .97343, v12 = .0414;
double v20 = .00886, v21 = .0405, v22 = .99914;
double Vckm[3][3] = {{v00, v01, v02}, {v10,v11,v12}, {v20,v21,v22}};
double Vckm_err[3][3] = {{.00014,.00061,.00015},{.00061,.00015,.0012},{.00033,.0011,.00005}};

//KR as a function of the IR scale
double kr(double lam) {
	return log(M_pl/lam)/M_PI;
}

//Quark yukawa couplings as a function of 5d parameters c
double yuk(double cl, double cr, double lam) {
	return exp((1-cl-cr)*M_PI*kr(lam))*sqrt((.5 - cl)/(exp((1 - 2*cl)*M_PI*kr(lam)) - 1 ))*sqrt((.5 - cr)/(exp((1 - 2*cr)*M_PI*kr(lam)) - 1 ));
}

//neutrino mass parameters (weinberg operators)
double vyuk(double ci, double cj, double lam) {
	return (1/M_pl)*exp(M_PI*kr(lam)*(2 - ci - cj))*sqrt((.5 - ci)/(exp((1 - 2*ci)*M_PI*kr(lam)) - 1 ))*sqrt((.5 - cj)/(exp((1 - 2*cj)*M_PI*kr(lam)) - 1 ));
}

double Chisq(const double *par) {

	//5d c parameters
	Vector3d cq(par[0],par[1],par[2]);
	Vector3d cu(par[3],par[4],par[5]);
	Vector3d cd(par[6],par[7],par[8]);

	//5d yukawa matrices
	Matrix3d yu;
	yu << par[9], par[10], par[11],
		   par[12], par[13], par[14],
		   par[15], par[16], par[17];
	Matrix3d yd;
	yd << par[18], par[19], par[20],
		   par[21], par[22], par[23],
		   par[24], par[25], par[26];

	Matrix3d cu5;
	Matrix3d cd5;
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			cu5(i,j) = yuk(cq(i),cu(j),lambda);
			cd5(i,j) = yuk(cq(i),cd(j),lambda);
		}
	}

	yu = yu.cwiseProduct(cu5);
	yd = yd.cwiseProduct(cd5);	
	
	JacobiSVD<Matrix3d> u_svd(yu,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3d> d_svd(yd,ComputeFullU | ComputeFullV);

	Matrix3d Vp = u_svd.matrixU().adjoint()*d_svd.matrixU();
	double yu_pred[3] = {u_svd.singularValues()[2],u_svd.singularValues()[1],u_svd.singularValues()[0]};
	double yd_pred[3] = {d_svd.singularValues()[2],d_svd.singularValues()[1],d_svd.singularValues()[0]};


	double chisq = 0;
	for(int i = 0; i < 3; i++) {
		chisq += pow((yu_pred[i] - yu_dat[i])/yu_dat_err[i],2);
		chisq += pow((yd_pred[i] - yd_dat[i])/yd_dat_err[i],2);
	}

	double lam_ckm_pred = Vp(0,1)/sqrt(Vp(0,0)*Vp(0,0) + Vp(0,1)*Vp(0,1));
	double rhoeta_ckm_pred = Vp(0,2)*(Vp(0,0)*Vp(0,0) + Vp(0,1)*Vp(0,1))/(Vp(1,2)*Vp(0,1));
	double A_ckm_pred = Vp(1,2)*sqrt(Vp(0,0)*Vp(0,0) + Vp(0,1)*Vp(0,1))/(Vp(0,1)*Vp(0,1));

	//chisq += pow((lam_ckm_pred - lam_ckm)/lam_ckm_err,2) + pow((A_ckm_pred - A_ckm)/A_ckm_err,2) + pow((rhoeta_ckm_pred - rhoeta_ckm)/rhoeta_ckm_err,2);

	
	for(int i = 0; i < 3; i ++) {
		for(int j = 0; j < 3; j++) {
			chisq += pow((Vp(i,j) - Vckm[i][j])/Vckm_err[0][0],2);
		}
	}
	

	//chisq += pow((Vp(0,1) - Vckm[0][1])/Vckm_err[0][1],2) + pow((Vp(0,2) - Vckm[0][2])/Vckm_err[0][2],2) + pow((Vp(1,2) - Vckm[1][2])/Vckm_err[1][2],2);

	return chisq;

}


std::vector<double> Minimize(const unsigned int input_seed = 0) {

	std::vector<double> output;
	unsigned int current_seed;

	if(input_seed != 0) current_seed = input_seed;
	else current_seed = rd();
	rng.seed(current_seed);
    std::uniform_real_distribution<> lams(-3.0,3.0);

	ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
   	minimum->SetMaxIterations(10000);  // for GSL
   	minimum->SetTolerance(0.01);
   	minimum->SetPrintLevel(0);

   	ROOT::Math::Functor f(&Chisq,27);
   	minimum->SetFunction(f);

   	minimum->SetVariable(0,"cq1",.7,.01);
	minimum->SetVariable(1,"cq2",.7,.01);
	minimum->SetVariable(2,"cq3",.4,.01);
	minimum->SetVariable(3,"cu1",.7,.01);
	minimum->SetVariable(4,"cu2",.7,.01);
	minimum->SetVariable(5,"cu3",-1,.01);
	minimum->SetVariable(6,"cd1",.7,.01);
	minimum->SetVariable(7,"cd2",.7,.01);
	minimum->SetVariable(8,"cd3",.7,.01);
	minimum->SetVariable(9,"yu1",lams(rng),.1);
	minimum->SetVariable(10,"yu2",lams(rng),.1);
	minimum->SetVariable(11,"yu3",lams(rng),.1);
	minimum->SetVariable(12,"yu4",lams(rng),.1);
	minimum->SetVariable(13,"yu5",lams(rng),.1);
	minimum->SetVariable(14,"yu6",lams(rng),.1);
	minimum->SetVariable(15,"yu7",lams(rng),.1);
	minimum->SetVariable(16,"yu8",lams(rng),.1);
	minimum->SetVariable(17,"yu9",lams(rng),.1);
	minimum->SetVariable(18,"yd1",lams(rng),.1);
	minimum->SetVariable(19,"yd2",lams(rng),.1);
	minimum->SetVariable(20,"yd3",lams(rng),.1);
	minimum->SetVariable(21,"yd4",lams(rng),.1);
	minimum->SetVariable(22,"yd5",lams(rng),.1);
	minimum->SetVariable(23,"yd6",lams(rng),.1);
	minimum->SetVariable(24,"yd7",lams(rng),.1);
	minimum->SetVariable(25,"yd8",lams(rng),.1);
	minimum->SetVariable(26,"yd9",lams(rng),.1);

	minimum->Minimize();

	if (input_seed != 0) {

		const double *x = minimum->X();

		//5d yukawa matrices
		Matrix3d yu;
		yu << x[9], x[10], x[11],
			   x[12], x[13], x[14],
			   x[15], x[16], x[17];
		Matrix3d yd;
		yd << x[18], x[19], x[20],
			   x[21], x[22], x[23],
			   x[24], x[25], x[26];

		Matrix3d cu5;
		cu5 << yuk(x[0],x[3],lambda), yuk(x[0],x[4],lambda), yuk(x[0],x[5],lambda),
				yuk(x[1],x[3],lambda), yuk(x[1],x[4],lambda), yuk(x[1],x[5],lambda),
				yuk(x[2],x[3],lambda), yuk(x[2],x[4],lambda), yuk(x[2],x[5],lambda);
		Matrix3d cd5;
		cd5 << yuk(x[0],x[6],lambda), yuk(x[0],x[7],lambda), yuk(x[0],x[8],lambda),
				yuk(x[1],x[6],lambda), yuk(x[1],x[7],lambda), yuk(x[1],x[8],lambda),
				yuk(x[2],x[6],lambda), yuk(x[2],x[7],lambda), yuk(x[2],x[8],lambda);

		yu = yu.cwiseProduct(cu5);
		yd = yd.cwiseProduct(cd5);	
		
		JacobiSVD<Matrix3d> u_svd(yu,ComputeFullU | ComputeFullV);
		JacobiSVD<Matrix3d> d_svd(yd,ComputeFullU | ComputeFullV);

		Matrix3d V_pred = u_svd.matrixU().adjoint()*d_svd.matrixU();
		double mu_pred[3] = {y2m(u_svd.singularValues()[2]),y2m(u_svd.singularValues()[1]),y2m(u_svd.singularValues()[0])};
		double md_pred[3] = {y2m(d_svd.singularValues()[2]),y2m(d_svd.singularValues()[1]),y2m(d_svd.singularValues()[0])};

		//Test prediction for quark masses
		cout << "Predicted up quark mass " << mu_pred[0]*1000 << " MeV." << "\n";
		cout << "Observed up quark mass " << y2m(yu_dat[0])*1000 << " MeV." << "\n";
		cout << "\n";
		cout << "Predicted down quark mass " << md_pred[0]*1000 << " MeV." << "\n";
		cout << "Observed down quark mass " << y2m(yd_dat[0])*1000 << " MeV." << "\n";
		cout << "\n";
		cout << "Predicted strange quark mass " << md_pred[1]*1000 << " MeV." << "\n";
		cout << "Observed strange quark mass " << y2m(yd_dat[1])*1000 << " MeV." << "\n";
		cout << "\n";
		cout << "Predicted charm quark mass " << mu_pred[1] << " GeV." << "\n";
		cout << "Observed charm quark mass " << y2m(yu_dat[1]) << " GeV." << "\n";
		cout << "\n";
		cout << "Predicted bottom quark mass " << md_pred[2] << " GeV." << "\n";
		cout << "Observed bottom quark mass " << y2m(yd_dat[2]) << " GeV." << "\n";
		cout << "\n";
		cout << "Predicted top quark mass " << mu_pred[2] << " GeV." << "\n";
		cout << "Observed top quark mass " << y2m(yu_dat[2]) << " GeV." << "\n";
		cout << "\n";
		cout << "Predicted CKM Matrix = " << V_pred; 
		cout << ".\n";
	}

	output.push_back(minimum->MinValue());
	output.push_back(current_seed);

	return output;

}

int main(int num_rolls) {

	double min_chisq = 1000000000;
	unsigned int min_seed;
	std::vector<double> output;
	int i = 0;
	while (i < num_rolls) {
		output = Minimize();
		if (output[0] < min_chisq) {
			min_chisq = output[0];
			min_seed = output[1];
		}
		i++;
	}

	Minimize(min_seed);
	return 0;
}