#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>

using std::cout;
double M_pl = 2.435e18;
double vev = 246;

const int num_lambdas = 20;
double lambdas[num_lambdas] = {200.0, 1052.03, 5533.84, 29108.9, 153117.0, 805419.0, 4.23663e6, 2.22853e7, 1.17224e8, 6.16617e8, 3.2435e9, 1.70613e10, 8.97452e10, 4.72073e11, 2.48318e12, 1.30619e13, 6.87076e13, 3.61412e14, 1.90109e15, 1.0e16};

double kr(double lam) {
  return log(M_pl/lam)/M_PI;
}

double Ninv(double c, double lambda) {
  return sqrt((.5 - c)/(exp((1 - 2*c)*M_PI*kr(lambda)) - 1));
}

double Proton_Decay_Sup(double c1,double c2,double c3,double c4, double lambda) {

  double Ninvall = Ninv(c1,lambda)*Ninv(c2,lambda)*Ninv(c3,lambda)*Ninv(c4,lambda);
  return sqrt(
    ((2-2*exp(-2*M_PI*kr(lambda)))/pow(M_pl,2))*Ninvall*(exp((4 - c1 - c2 - c3 - c4)*M_PI*kr(lambda)) - 1)/(4 - c1 - c2 -c3 - c4)
  );
}

int main(int argc, char* argv[]) {

  std::string file = argv[1];
  auto infile = TFile::Open(file.c_str(),"Update");

  TTree* tree[num_lambdas];
  double cq1,cq2,cq3;
  double cu1,cu2,cu3;
  double cd1,cd2,cd3;
  double chisq;

  double qqql_proton_decay_suppression;
  double duql_proton_decay_suppression;

  for (int i = 0; i < num_lambdas; i++) {
    tree[i] = (TTree*)infile->Get(std::to_string(i).c_str());
    tree[i]->SetBranchAddress("cq1",&cq1);tree[i]->SetBranchAddress("cq2",&cq2);tree[i]->SetBranchAddress("cq3",&cq3);
    tree[i]->SetBranchAddress("cu1",&cu1);tree[i]->SetBranchAddress("cu2",&cu2);tree[i]->SetBranchAddress("cu3",&cu3);
    tree[i]->SetBranchAddress("cd1",&cd1);tree[i]->SetBranchAddress("cd2",&cd2);tree[i]->SetBranchAddress("cd3",&cd3);
    tree[i]->SetBranchAddress("chi_squared", &chisq);
    TBranch *qqqls = tree[i]->Branch("qqql_suppression", &qqql_proton_decay_suppression);
    TBranch *duqls = tree[i]->Branch("duql_suppression", &duql_proton_decay_suppression);

    //Loop over all the entries in each tree and calculate the proton decay suppression scale
    Long64_t nentries = tree[i]->GetEntries();
    for (Long64_t j = 0; j < nentries; j++) {
      tree[i]->GetEntry(j);
      qqql_proton_decay_suppression = Proton_Decay_Sup(cq1,cq1,cq1,1,lambdas[i]);
      duql_proton_decay_suppression = Proton_Decay_Sup(cd1,cu1,cq1,1,lambdas[i]);
      qqqls->Fill();
      duqls->Fill();
    }

    tree[i]->Write();
  }

  infile->Close();

  return 0;

}