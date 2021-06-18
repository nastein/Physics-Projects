#include "PlotUtils.C"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TFile.h>
#include <TApplication.h>
#include <stdlib.h>

using namespace std;

void fcc_plotter(string *files, int num_files, std::vector<double> norms)
{

  gStyle->SetOptStat(0);


  Plot plot4("costhetastar");
  plot4.Title("Reconstructed leptonic Cos(#theta *)");

  Plot plot5("deltaphi");
  plot5.Title("Reconstructed #Delta#phi");

  int colors[8] = {1,2,3,4,6,8,40,47};

  for (int i = 0; i < num_files; i++) {

    TFile *Xsm = new TFile((files[i] + "_hists.root").c_str());

    plot4.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("Costheta_Xsm_" + files[i]).c_str()), true, norms[i]*1E3*3000);
    plot5.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("Delta_phi_Xsm_" + files[i]).c_str()), true, norms[i]*1E3*3000);

  }

  plot4.Draw("Cos(#theta *)", "Events / bin", 0.60, 1.2);
  plot5.Draw("#Delta#phi", "Events / bin", 0.60, 1.2);


}