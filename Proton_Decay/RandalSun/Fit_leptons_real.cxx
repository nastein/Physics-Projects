#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include "reader.h"
#include "RS_lepton_real_fitter.h"

using std::cout;

const int num_lambdas = 1;
//const int num_lambdas = 20;
//double lambdas[num_lambdas] = {200.0, 1052.03, 5533.84, 29108.9, 153117.0, 805419.0, 4.23663e6, 2.22853e7, 1.17224e8, 6.16617e8, 3.2435e9, 1.70613e10, 8.97452e10, 4.72073e11, 2.48318e12, 1.30619e13, 6.87076e13, 3.61412e14, 1.90109e15, 1.0e16};

double lambdas[num_lambdas] = {200.0};


int main(int argc, char const *argv[]) {

  int num_rolls = atof(argv[1]);
  std::string output_file = argv[2];
  bool penalty;

  if (strcmp(argv[3],"true") == 0) {penalty = true; cout << "Penalizing yukawa couplings outside of acceptable range!" << "\n";}
  else {penalty = false; cout << "Not penalizing yukawa couplings at all" << "\n";}

  TString output = (output_file).c_str();
  TFile *outfile = new TFile(output,"RECREATE");

  std::vector<double> ms_organized;
  std::vector<double> ms_err_organized;
  int l = 0;

  RSFitter fitter(l, num_rolls, lambdas[l], ms_organized, ms_err_organized);
  fitter.turn_on_penalty(penalty);

  for(l = 0; l < num_lambdas; l++) {

    //for(int j = 0; j < 6; j++) {
      //ms_organized.push_back(masses[j][l]);
      //ms_err_organized.push_back(masses_err[j][l]);
    //}

    TTree *out_tree = fitter.fit(l, lambdas[l], ms_organized, ms_err_organized);
    out_tree->Write();

    ms_organized.clear(); ms_err_organized.clear();
  }

  outfile->Close();

  return 0;
}
