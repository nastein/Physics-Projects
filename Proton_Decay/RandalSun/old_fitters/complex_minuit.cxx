#include "TRandom2.h"
#include "TError.h"
#include "TMinuit.h"
#include "TString.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include </Users/noahsteinberg/Physics/eigen/Eigen/SVD>
#include </Users/noahsteinberg/Physics/eigen/Eigen/Dense>
#include <cstdlib>
#include <random>
#include <complex>

#define cp std::complex<double> 

typedef std::mt19937 MyRNG;
std::random_device rd;

MyRNG rng;

using namespace Eigen;
using std::cout;
using std::abs;

//Initialize some useful parameters
double M_pl = 2.435e18;
double vev = 246;
double lambda = 162.9;

double lam_ckm = .2257, lam_ckm_err = .0010;
double A_ckm = .814, A_ckm_err = .022;
double rho_ckm = .135, rho_ckm_err = .031;
double eta_ckm = .349, eta_ckm_err = .017;
double rhoeta_ckm = .349, rhoeta_ckm_err = .017;

std::complex<double> div(cp a, cp b) {
	double real = (a.real()*b.real() + a.imag()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));
	double imaginary = (a.imag()*b.real() - a.real()*b.imag())/(pow(b.real(),2) + pow(b.imag(),2));

	return cp(real, imaginary);
}

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
double yd_dat[3] ={m2y(.00276),m2y(.052),m2y(2.75)};
double yd_dat_err[3] = {m2y(.00117),m2y(.015),m2y(.09)};

//CKM matrix
double Vckm[3][3] = {{.97427, .22536, .00355}, {.22522,.97343,.0414}, {.00886,.0405,.99914}};
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

double Chisq(Double_t *par) {

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
			cu5(i,j) = yuk(cq(i),cu(j),lambda);
			cd5(i,j) = yuk(cq(i),cd(j),lambda);
		}
	}

	yu = yu.cwiseProduct(cu5);
	yd = yd.cwiseProduct(cd5);	
	
	JacobiSVD<Matrix3cd> u_svd(yu,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3cd> d_svd(yd,ComputeFullU | ComputeFullV);

	Matrix3cd Vp = u_svd.matrixU().adjoint()*d_svd.matrixU();
	double yu_pred[3] = {u_svd.singularValues()[2],u_svd.singularValues()[1],u_svd.singularValues()[0]};
	double yd_pred[3] = {d_svd.singularValues()[2],d_svd.singularValues()[1],d_svd.singularValues()[0]};


	double chisq = 0;
	for(int i = 0; i < 3; i++) {
		chisq += pow((yu_pred[i] - yu_dat[i])/yu_dat_err[i],2);
		chisq += pow((yd_pred[i] - yd_dat[i])/yd_dat_err[i],2);
	}

	double lam_ckm_pred = abs(Vp(0,1))/sqrt(pow(abs(Vp(0,0)),2) + pow(abs(Vp(0,1)),2));
	double A_ckm_pred = abs(div(Vp(1,2),Vp(0,1)))/lam_ckm_pred;
	double rho_ckm_pred = Vp(0,2).real()/(A_ckm_pred*pow(lam_ckm_pred,3));
	double eta_ckm_pred = -Vp(0,2).imag()/(A_ckm_pred*pow(lam_ckm_pred,3));

	chisq += pow((lam_ckm_pred - lam_ckm)/lam_ckm_err,2) + pow((A_ckm_pred - A_ckm)/A_ckm_err,2) + pow((rho_ckm_pred - rho_ckm)/rho_ckm_err,2) + pow((eta_ckm_pred - eta_ckm)/eta_ckm_err,2);

	return chisq;

}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
	f = Chisq(par);
}

std::vector<double> Minimize(const unsigned int input_seed = 0) {
	
	std::vector<double> output;
	unsigned int current_seed;

	TMinuit *minuit = new TMinuit(45);
	minuit->SetPrintLevel(-1);
	minuit->SetFCN(fcn);

	Double_t arglist[30];
	Int_t ierflg = 0;

	if(input_seed != 0) current_seed = input_seed;
	else current_seed = rd();
	
	rng.seed(current_seed);

    std::uniform_real_distribution<> cs(-1.0,1.0);
    std::uniform_real_distribution<> lams(-3.0,3.0);

    /*
	minuit->mnparm(0,"cq1",.7,.01,0,0,ierflg);
	minuit->mnparm(1,"cq2",.7,.01,0,0,ierflg);
	minuit->mnparm(2,"cq3",.4,.01,0,0,ierflg);
	minuit->mnparm(3,"cu1",.7,.01,0,0,ierflg);
	minuit->mnparm(4,"cu2",.7,.01,0,0,ierflg);
	minuit->mnparm(5,"cu3",-1,.01,0,0,ierflg);
	minuit->mnparm(6,"cd1",.7,.01,0,0,ierflg);
	minuit->mnparm(7,"cd2",.7,.01,0,0,ierflg);
	minuit->mnparm(8,"cd3",.7,.01,0,0,ierflg);*/
	minuit->mnparm(0,"cq1",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(1,"cq2",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(2,"cq3",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(3,"cu1",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(4,"cu2",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(5,"cu3",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(6,"cd1",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(7,"cd2",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(8,"cd3",cs(rng),.01,0,0,ierflg);
	minuit->mnparm(9,"yu1",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(10,"yu2",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(11,"yu3",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(12,"yu4",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(13,"yu5",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(14,"yu6",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(15,"yu7",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(16,"yu8",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(17,"yu9",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(18,"yd1",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(19,"yd2",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(20,"yd3",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(21,"yd4",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(22,"yd5",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(23,"yd6",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(24,"yd7",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(25,"yd8",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(26,"yd9",lams(rng),.1,-3,3,ierflg);

	minuit->mnparm(27,"yu1i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(28,"yu2i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(29,"yu3i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(30,"yu4i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(31,"yu5i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(32,"yu6i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(33,"yu7i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(34,"yu8i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(35,"yu9i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(36,"yd1i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(37,"yd2i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(38,"yd3i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(39,"yd4i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(40,"yd5i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(41,"yd6i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(42,"yd7i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(43,"yd8i",lams(rng),.1,-3,3,ierflg);
	minuit->mnparm(44,"yd9i",lams(rng),.1,-3,3,ierflg);

	//arglist[0] = 500;
	arglist[0] = 1;

	minuit->mnexcm("MIGRAD", arglist, 0, ierflg);

	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	if (input_seed != 0) {

		double idc;
		double cq1,cq2,cq3,cu1,cu2,cu3,cd1,cd2,cd3;
		double yu1,yu2,yu3,yu4,yu5,yu6,yu7,yu8,yu9;
		double yd1,yd2,yd3,yd4,yd5,yd6,yd7,yd8,yd9;
		double yu1i,yu2i,yu3i,yu4i,yu5i,yu6i,yu7i,yu8i,yu9i;
		double yd1i,yd2i,yd3i,yd4i,yd5i,yd6i,yd7i,yd8i,yd9i;

		minuit->GetParameter(0,cq1,idc);
		minuit->GetParameter(1,cq2,idc);
		minuit->GetParameter(2,cq3,idc);
		minuit->GetParameter(3,cu1,idc);
		minuit->GetParameter(4,cu2,idc);
		minuit->GetParameter(5,cu3,idc);
		minuit->GetParameter(6,cd1,idc);
		minuit->GetParameter(7,cd2,idc);
		minuit->GetParameter(8,cd3,idc);
		minuit->GetParameter(9,yu1,idc);
		minuit->GetParameter(10,yu2,idc);
		minuit->GetParameter(11,yu3,idc);
		minuit->GetParameter(12,yu4,idc);
		minuit->GetParameter(13,yu5,idc);
		minuit->GetParameter(14,yu6,idc);
		minuit->GetParameter(15,yu7,idc);
		minuit->GetParameter(16,yu8,idc);
		minuit->GetParameter(17,yu9,idc);
		minuit->GetParameter(18,yd1,idc);
		minuit->GetParameter(19,yd2,idc);
		minuit->GetParameter(20,yd3,idc);
		minuit->GetParameter(21,yd4,idc);
		minuit->GetParameter(22,yd5,idc);
		minuit->GetParameter(23,yd6,idc);
		minuit->GetParameter(24,yd7,idc);
		minuit->GetParameter(25,yd8,idc);
		minuit->GetParameter(26,yd9,idc);
		minuit->GetParameter(27,yu1i,idc);
		minuit->GetParameter(28,yu2i,idc);
		minuit->GetParameter(29,yu3i,idc);
		minuit->GetParameter(30,yu4i,idc);
		minuit->GetParameter(31,yu5i,idc);
		minuit->GetParameter(32,yu6i,idc);
		minuit->GetParameter(33,yu7i,idc);
		minuit->GetParameter(34,yu8i,idc);
		minuit->GetParameter(35,yu9i,idc);
		minuit->GetParameter(36,yd1i,idc);
		minuit->GetParameter(37,yd2i,idc);
		minuit->GetParameter(38,yd3i,idc);
		minuit->GetParameter(39,yd4i,idc);
		minuit->GetParameter(40,yd5i,idc);
		minuit->GetParameter(41,yd6i,idc);
		minuit->GetParameter(42,yd7i,idc);
		minuit->GetParameter(43,yd8i,idc);
		minuit->GetParameter(44,yd9i,idc);
		
		//5d yukawa matrices
		Matrix3cd yu;
		yu << cp(yu1,yu1i), cp(yu2,yu2i), cp(yu3,yu3i),
			   cp(yu4,yu4i), cp(yu5,yu5i), cp(yu6,yu6i),
			   cp(yu7,yu7i), cp(yu8,yu8i), cp(yu9,yu9i);
		Matrix3cd yd;
		yd << cp(yd1,yd1i), cp(yd2,yd2i), cp(yd3,yd3i),
			   cp(yd4,yd4i), cp(yd5,yd5i), cp(yd6,yd6i),
			   cp(yd7,yd7i), cp(yd8,yd8i), cp(yd9,yd9i);

		Matrix3cd cu5;
		cu5 << yuk(cq1,cu1,lambda), yuk(cq1,cu2,lambda), yuk(cq1,cu3,lambda),
				yuk(cq2,cu1,lambda), yuk(cq2,cu2,lambda), yuk(cq2,cu3,lambda),
				yuk(cq3,cu1,lambda), yuk(cq3,cu2,lambda), yuk(cq3,cu3,lambda);
		Matrix3cd cd5;
		cd5 << yuk(cq1,cd1,lambda), yuk(cq1,cd2,lambda), yuk(cq1,cd3,lambda),
				yuk(cq2,cd1,lambda), yuk(cq2,cd2,lambda), yuk(cq2,cd3,lambda),
				yuk(cq3,cd1,lambda), yuk(cq3,cd2,lambda), yuk(cq3,cd3,lambda);


		yu = yu.cwiseProduct(cu5);
		yd = yd.cwiseProduct(cd5);	
		
		JacobiSVD<Matrix3cd> u_svd(yu,ComputeFullU | ComputeFullV);
		JacobiSVD<Matrix3cd> d_svd(yd,ComputeFullU | ComputeFullV);

		Matrix3cd Vp = u_svd.matrixU().adjoint()*d_svd.matrixU();
		double mu_pred[3] = {y2m(u_svd.singularValues()[2]),y2m(u_svd.singularValues()[1]),y2m(u_svd.singularValues()[0])};
		double md_pred[3] = {y2m(d_svd.singularValues()[2]),y2m(d_svd.singularValues()[1]),y2m(d_svd.singularValues()[0])};

		Matrix3d V_pred(3,3);
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				V_pred(i,j) = abs(Vp(i,j));
			}
		}


		//Test prediction for quark masses
		cout << "\n";
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
		cout << "Predicted CKM Matrix = \n" << V_pred; 
		cout << ".\n";
		cout << "The final Chisq was " << amin << "\n";
	
	}
	
	output.push_back(amin);
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





