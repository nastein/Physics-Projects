#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include "reader.h"
#include "RS_quark_fitter.h"

using std::cout;

const int num_lambdas = 20;
double lambdas[num_lambdas] = {200.0, 1052.03, 5533.84, 29108.9, 153117.0, 805419.0, 4.23663e6, 2.22853e7, 1.17224e8, 6.16617e8, 3.2435e9, 1.70613e10, 8.97452e10, 4.72073e11, 2.48318e12, 1.30619e13, 6.87076e13, 3.61412e14, 1.90109e15, 1.0e16};

CSVReader reader("fermion_mass_data/quark_masses.csv");
CSVReader reader_err("fermion_mass_data/quark_masses_err.csv");

std::vector<double> mu, mu_err;
std::vector<double> md, md_err;
std::vector<double> ms, ms_err;
std::vector<double> mc, mc_err;
std::vector<double> mb, mb_err;
std::vector<double> mt, mt_err;
std::vector<std::vector<double>> masses;
std::vector<std::vector<double>> masses_err;

std::vector<std::vector<std::string>> masses_str = reader.getData();
std::vector<std::vector<std::string>> masses_str_err = reader_err.getData();


int main(int argc, char const *argv[]) {

	for (int l = 0; l < 6; l++) {
		for (int k = 0; k < num_lambdas; k++) {
			if (l == 0) {
				mu.push_back(atof(masses_str[l][k].c_str()));
				mu_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
			if (l == 1) {
				md.push_back(atof(masses_str[l][k].c_str()));
				md_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
			if (l == 2) {
				ms.push_back(atof(masses_str[l][k].c_str()));
				ms_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
			if (l == 3) {
				mc.push_back(atof(masses_str[l][k].c_str()));
				mc_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
			if (l == 4) {
				mb.push_back(atof(masses_str[l][k].c_str()));
				mb_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
			if (l == 5) {
				mt.push_back(atof(masses_str[l][k].c_str()));
				mt_err.push_back(atof(masses_str_err[l][k].c_str()));
			}
		}
	}

	masses.push_back(mu); masses.push_back(md); masses.push_back(ms);
	masses.push_back(mc); masses.push_back(mb); masses.push_back(mt);
	masses_err.push_back(mu_err); masses_err.push_back(md_err); masses_err.push_back(ms_err);
	masses_err.push_back(mc_err); masses_err.push_back(mb_err); masses_err.push_back(mt_err);

	int num_rolls = atof(argv[1]);
	std::string output_file = argv[2];
	bool penalty = bool(argv[3]);

	if (penalty == true) cout << "Penalizing yukawa couplings outside of acceptable range!" << "\n";
	else cout << "Not penalizing yukawa couplings at all" << "\n";

	TString output = (output_file).c_str();
	TFile *outfile = new TFile(output,"RECREATE");

	std::vector<double> ms_organized;
	std::vector<double> ms_err_organized;
	int l = 0;

	//Create the fitter
	RSFitter fitter(l, num_rolls, lambdas[l], ms_organized, ms_err_organized);
	fitter.turn_on_penalty(penalty);

	for(l = 0; l < num_lambdas; l++) {

		for(int j = 0; j < 6; j++) {
			ms_organized.push_back(masses[j][l]);
			ms_err_organized.push_back(masses_err[j][l]);
		}

		TTree *out_tree = fitter.fit(l, lambdas[l], ms_organized, ms_err_organized);
		out_tree->Write();

		ms_organized.clear(); ms_err_organized.clear();
	}

	outfile->Close();

	return 0;
}












