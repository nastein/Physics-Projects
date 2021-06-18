#include "TRandom2.h"
#include "TError.h"
#include "TMinuit.h"
#include "TString.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include </Users/noahsteinberg/Physics/eigen/Eigen/SVD>
#include </Users/noahsteinberg/Physics/eigen/Eigen/Dense>

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

double yuk_logs(double cl,double cr, double lam) {
	return sqrt(((cl*cr/exp(1)),2*M_PI*kr(lambda))*(.5 - log(cl))*(.5 - log(cr))/(((1-Pow((cl/exp(1)),M_PI*kar(lam)))*(1-(cr/exp(1))**(M_PI*kr(lam))))**2));
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
	Matrix3f yu;
	yu << par[9], par[10], par[11],
		   par[12], par[13], par[14],
		   par[15], par[16], par[17];
	Matrix3f yd;
	yd << par[18], par[19], par[20],
		   par[21], par[22], par[23],
		   par[24], par[25], par[26];

	Matrix3f cu5;
	Matrix3f cd5;
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			cu5(i,j) = yuk_logs(cq(i),cu(j),lambda);
			cd5(i,j) = yuk_logs(cq(i),cd(j),lambda);
		}
	}

	yu = yu.cwiseProduct(cu5);
	yd = yd.cwiseProduct(cd5);	
	
	JacobiSVD<Matrix3f> u_svd(yu,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3f> d_svd(yd,ComputeFullU | ComputeFullV);

	Matrix3f Vp = u_svd.matrixU().adjoint()*d_svd.matrixU();
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

	/*
	for(int i = 0; i < 3; i ++) {
		for(int j = 0; j < 3; j++) {
			chisq += pow((Vp(i,j) - Vckm[i][j])/Vckm_err[0][0],2);
		}
	}*/
	

	return chisq;

}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
	f = Chisq(par);
}


int Minimize() {
	
	TMinuit *minuit = new TMinuit(27);
	minuit->SetPrintLevel();
	minuit->SetFCN(fcn);

	Double_t arglist[30];
	Int_t ierflg = 0;

	minuit->mnparm(0,"cq1",.7,.01,0,0,ierflg);
	minuit->mnparm(1,"cq2",.7,.01,0,0,ierflg);
	minuit->mnparm(2,"cq3",.4,.01,0,0,ierflg);
	minuit->mnparm(3,"cu1",.7,.01,0,0,ierflg);
	minuit->mnparm(4,"cu2",.7,.01,0,0,ierflg);
	minuit->mnparm(5,"cu3",-1,.01,0,0,ierflg);
	minuit->mnparm(6,"cd1",.7,.01,0,0,ierflg);
	minuit->mnparm(7,"cd2",.7,.01,0,0,ierflg);
	minuit->mnparm(8,"cd3",.7,.01,0,0,ierflg);
	minuit->mnparm(9,"yu1",1.5,.2,-2,2,ierflg);
	minuit->mnparm(10,"yu2",1.5,.2,-2,2,ierflg);
	minuit->mnparm(11,"yu3",1.5,.2,-2,2,ierflg);
	minuit->mnparm(12,"yu4",1.5,.2,-2,2,ierflg);
	minuit->mnparm(13,"yu5",1.5,.2,-2,2,ierflg);
	minuit->mnparm(14,"yu6",1.5,.2,-2,2,ierflg);
	minuit->mnparm(15,"yu7",1.5,.2,-2,2,ierflg);
	minuit->mnparm(16,"yu8",1.5,.2,-2,2,ierflg);
	minuit->mnparm(17,"yu9",1.5,.2,-2,2,ierflg);
	minuit->mnparm(18,"yd1",1.5,.2,-2,2,ierflg);
	minuit->mnparm(19,"yd2",1.5,.2,-2,2,ierflg);
	minuit->mnparm(20,"yd3",1.5,.2,-2,2,ierflg);
	minuit->mnparm(21,"yd4",1.5,.2,-2,2,ierflg);
	minuit->mnparm(22,"yd5",1.5,.2,-2,2,ierflg);
	minuit->mnparm(23,"yd6",1.5,.2,-2,2,ierflg);
	minuit->mnparm(24,"yd7",1.5,.2,-2,2,ierflg);
	minuit->mnparm(25,"yd8",1.5,.2,-2,2,ierflg);
	minuit->mnparm(26,"yd9",1.5,.2,-2,2,ierflg);

	arglist[0] = 1;

	minuit->mnexcm("MIGRAD", arglist, 0, ierflg);

	double idc;
	double cq1,cq2,cq3,cu1,cu2,cu3,cd1,cd2,cd3;
	double yu1,yu2,yu3,yu4,yu5,yu6,yu7,yu8,yu9;
	double yd1,yd2,yd3,yd4,yd5,yd6,yd7,yd8,yd9;
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
	
	//5d yukawa matrices
	Matrix3f yu;
	yu << yu1, yu2, yu3,
		   yu4, yu5, yu6,
		   yu7, yu8, yu9;
	Matrix3f yd;
	yd << yd1, yd2, yd3,
		   yd4, yd5, yd6,
		   yd7, yd8, yd9;

	Matrix3f cu5;
	cu5 << yuk_logs(cq1,cu1,lambda), yuk_logs(cq1,cu2,lambda), yuk_logs(cq1,cu3,lambda),
			yuk_logs(cq2,cu1,lambda), yuk_logs(cq2,cu2,lambda), yuk_logs(cq2,cu3,lambda),
			yuk_logs(cq3,cu1,lambda), yuk_logs(cq3,cu2,lambda), yuk_logs(cq3,cu3,lambda);
	Matrix3f cd5;
	cd5 << yuk_logs(cq1,cd1,lambda), yuk_logs(cq1,cd2,lambda), yuk_logs(cq1,cd3,lambda),
			yuk_logs(cq2,cd1,lambda), yuk_logs(cq2,cd2,lambda), yuk_logs(cq2,cd3,lambda),
			yuk_logs(cq3,cd1,lambda), yuk_logs(cq3,cd2,lambda), yuk_logs(cq3,cd3,lambda);


	yu = yu.cwiseProduct(cu5);
	yd = yd.cwiseProduct(cd5);	
	
	JacobiSVD<Matrix3f> u_svd(yu,ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3f> d_svd(yd,ComputeFullU | ComputeFullV);

	Matrix3f V_pred = u_svd.matrixU().adjoint()*d_svd.matrixU();
	double mu_pred[3] = {y2m(u_svd.singularValues()[2]),y2m(u_svd.singularValues()[1]),y2m(u_svd.singularValues()[0])};
	double md_pred[3] = {y2m(d_svd.singularValues()[2]),y2m(d_svd.singularValues()[1]),y2m(d_svd.singularValues()[0])};

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
	cout << "Predicted CKM Matrix = " << V_pred; 
	cout << ".\n";
	return 0;

}

int main(int argc, char *argv[]) {

	Minimize();
	return 0;
}