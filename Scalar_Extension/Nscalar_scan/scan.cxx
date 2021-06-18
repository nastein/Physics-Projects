#include <iostream>
#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <TVectorD.h>
#include <TVector.h>
#include <TImage.h>
#include <TFile.h>
#include <TTree.h>
#include "range.hpp"
#include </Users/noahsteinberg/Physics/eigen/Eigen/Dense>
#include </Users/noahsteinberg/Physics/eigen/Eigen/Eigenvalues>

int numrolls = 1e8;

typedef std::mt19937 MyRNG;
std::random_device rd;
MyRNG rng;

double vev = 246;
double eps[5] = {1e-2,5e-3,1e-3,5e-4,1e-4};
double higgs_mass_exp = 125.1;
double higgs_unc = 3*.14;

using std::cout;
using util::lang::range;
using util::lang::indices;
using namespace Eigen;

int main(int argc, char const *argv[]) {

  int N = atoi(argv[1]);//Number of extra scalars

  std::string output_file = argv[2];
  TString output = (output_file).c_str();
  TFile *outfile = new TFile(output,"RECREATE");

  TVector num(1); num[0] = N;
  num.Write("Number_of_Scalars");

  cout << "Scalar extension with " << N << " singlet scalars. \n";
  cout << "Creating output file " << output_file << "\n";

  cout.precision(8);

  double xi[N], kappa[N], epsilon[N], lambda;
  rng.seed(rd());
  std::uniform_real_distribution<> lambdas(-M_PI,M_PI);
  std::uniform_real_distribution<> epsilons(-1,1);
  std::uniform_real_distribution<> xis(1000,10000);

  int higgs_count;
  int higgs_num;
  int sols;
  double suppression = 0;

  bool PSD;

  //LDLT< Matrix<long double,N,N> > ldlt_solver;
  LDLT<MatrixXd> ldlt_solver;

  for(auto idx : indices(eps)) {
    TVectorD v(1); v[0] = eps[idx];
    std::string id_string = std::to_string(int(idx));
    TTree *outtree = new TTree(id_string.c_str(),id_string.c_str());
    outtree->GetUserInfo()->Add(&v);

    outtree->Branch("suppression", &suppression);

    sols = 0;
    cout << "Using epsilon_0 = " << eps[idx] << "\n";

    for(auto k : range (0,numrolls)) {
      suppression = 0;
      PSD = false;

      //Fill array of random values for our model parameters
      lambda = lambdas(rng);
      for(auto i : range (0,N)) {
        xi[i] = xis(rng);
        kappa[i] = lambdas(rng);
        epsilon[i] = epsilons(rng)*eps[idx];
      }
      //Fill matrix
      MatrixXd M(N,N);
      for(auto i : range(0,N)) {
        for (auto j : range(0,N)) {

          if(i == j) {
            if(i == 0) M(i,j) = 2*vev*vev*lambda;
            else M(i,j) = 2*xi[i]*xi[i]*kappa[i];
          }

          else if(i == 0 && j != 0) M(i,j) = vev*epsilon[j]*xi[j];
          else if((i != 0 && j == 0)) M(i,j) = M(j,i);
          else M(i,j) = 0;
        }
      }

      //We want to see if our parameters have given us a physical spectrum (i.e. all non negative masses). This is equivalent to the matrix being positive semi-definite (https://math.stackexchange.com/questions/1873559/which-non-negative-matrices-have-negative-eigenvalues).
      PSD = M.ldlt().isPositive();
      if (PSD == false) continue;

      //Diagonalize matrix
      SelfAdjointEigenSolver<MatrixXd> eigensolver(M);
      MatrixXd Hvec(N,1);

      higgs_num = 0;//eigenvector number
      higgs_count = 0;//number of eigenvectors with mass which agrees with observed higgs

      //Check if we have a physical mass spectrum and if any of the mass eigenstates match the mass of the observed higgs boson M_H = 125.09 +/- .35 GeV
      for(auto i : range (0,N)) {
        double mass_sq = eigensolver.eigenvalues()[i];
        if(sqrt(mass_sq) >= higgs_mass_exp - higgs_unc && sqrt(mass_sq) <= higgs_mass_exp + higgs_unc) {
          Hvec(i) = 1.0;
          higgs_num = i;
          higgs_count++;
        }
        else Hvec(i) = 0;
      }

      if(higgs_count == 0) continue;
      if(higgs_count > 1) cout << "WARNING!!!! FOUND " << higgs_count << " CANDIDATE HIGGSES \n";
      if(higgs_count == 1) {
        cout << "Found a higgs candidate in this parameter scan" << "\n";
        cout << "PSD = " << PSD << "\n";
        cout << "Mass matrix: " << "\n";
        cout << M << "\n";
        cout << "Hvec = " << "\n";
        cout << Hvec << "\n";

        //Matrix<long double, N, N> U = eigensolver.eigenvectors();
        //MatrixXd U(N,N) = ;
        cout << "U: " << "\n";
        cout << eigensolver.eigenvectors() << "\n";
        suppression = (eigensolver.eigenvectors()*Hvec)(0);
        cout << "Suppression factor = " << std::fixed << suppression <<  "\n";
        outtree->Fill();
        sols++;
      }
    }
    cout << sols << " solutions for epsilon_0 = " << eps[idx] << "\n";
    outtree->Write();
  }

  outfile->Close();

  return 0;
}