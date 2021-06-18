#include "PlotUtils.C"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TFile.h>
#include <TApplication.h>
#include <stdlib.h>

using namespace std;

void fcc_plotter(string *files, int num_files)
{

  gStyle->SetOptStat(0);

  /*
  Plot plot("TagMjj");
  plot.Title("Mjj of Tagging Jet Candidates");

  Plot plot2("Vhadmjj-Mw");
  plot2.Title("Difference in mass between W and Vhad candidate jets");

  Plot plot3("MwwT");
  plot3.Title("Reconstructed transverse mass of W boson pair");
  */

  Plot plot4("costhetastar");
  plot4.Title("Reconstructed leptonic Cos(#theta *)");

  int colors[8] = {1,2,3,4,6,8,40,47};

  for (int i = 0; i < num_files; i++) {

    TFile *Xsm = new TFile((files[i] + "_hists.root").c_str());

    //plot.Add((files[i]).c_str(),  colors[i], 1, (TH1D*)Xsm->Get(("tag_mjj_" + files[i]).c_str()));
    //plot2.Add((files[i]).c_str(),  colors[i], 1, (TH1D*)Xsm->Get(("VHad_mjj_diff_" + files[i]).c_str()));
    //plot3.Add((files[i]).c_str(),  colors[i], 1,  (TH1D*)Xsm->Get(("MwwT_Xsm_" + files[i]).c_str()));
    plot4.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("Costheta_Xsm_" + files[i]).c_str()));


  }


  //plot.Draw("M_{jj} (GeV)", "Events / bin", 0.60, 1.2);
  //plot2.Draw("M^{Vhad}_{jj} - M_{W} (GeV)", "Events / bin", 0.60, 1.2);
  //plot3.SetLog(true);plot3.Draw("M_{T}(WW) (GeV)", "Events / bin", 0.60, 1.2);
  plot4.Draw("Cos(#theta *)", "Events / bin", 0.60, 1.2);


}