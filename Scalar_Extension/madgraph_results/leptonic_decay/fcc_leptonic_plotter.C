#include "PlotUtils.C"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TFile.h>
#include <TApplication.h>
#include <stdlib.h>

using namespace std;

void fcc_leptonic_plotter(string *files, int num_files)
{
  //Just hard coding cross sections for now
  //Goes SM/2TeVpt1/2TeVpt2/5TeVpt1/5TeVpt2
  std::vector<double> norms = {1.0/112500.0, 1.0/90000.0, 1.0/110000.0, 1.0/110000.0, 1.0/107500.0};

  gStyle->SetOptStat(0);


  Plot plot5("deltaphi_jet");
  plot5.Title("Reconstructed #Delta#phi_{jj}");

  Plot plot6("deltaphi_lep");
  plot6.Title("Reconstructed #Delta#phi_{ll}");

  Plot plot7("M_{1}");
  plot7.Title("Mass variable");
  plot7.SetLog(true);

  Plot plot10("min lep pt phi ");
  plot10.Title("Number of events vs minimum leading lepton PT in Phi longitudinal region");
  plot10.SetLog(true);

  Plot plot11("min lep pt mass");
  plot11.Title("Number of events vs minimum leading lepton PT in M_{T} longitudinal region");
  plot11.SetLog(true);

  Plot plot12("MTww");
  plot12.Title("Transverse mass of WW system");
  plot12.SetLog(true);

  int colors[8] = {1,2,3,4,6,8,40,47};

  for (int i = 0; i < num_files; i++) {
    TFile *Xsm = new TFile((files[i] + "_hists.root").c_str());

    plot5.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("Delta_phi_jet_" + files[i]).c_str()), true, norms[i]);
    plot6.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("Delta_phi_lep_" + files[i]).c_str()), true, norms[i]);
    plot7.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("m_trans_" + files[i]).c_str()), true, norms[i]);
    plot10.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("min_PT_phi_" + files[i]).c_str()), true, norms[i]);
    plot11.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("min_PT_mass_" + files[i]).c_str()), true, norms[i]);
    plot12.Add((files[i]).c_str(), colors[i], 1, (TH1D*)Xsm->Get(("MTww_" + files[i]).c_str()), true, norms[i]);
  }

  plot5.Draw("#Delta#phi_{jj}", "Events / bin", 0.60, 1.2);
  plot6.Draw("#Delta#phi_{ll}", "Events / bin", 0.60, 1.2);
  plot7.Draw("M_{trans}", "Events / bin", 0.60, 1.2);
  plot10.DrawRatio("min PT_{l} (GeV)", "Events / bin", 0.60, 1.2);
  plot10.AddPlotLabel("|#Delta#phi| > 2.0", .7, .5, 0.05);
  plot11.DrawRatio("min PT_{l} (GeV)", "Events / bin", 0.60, 1.2);
  plot11.AddPlotLabel("M_{T} < 250 GeV", .7, .5, 0.05);
  plot12.DrawRatio("MTww (GeV)", "Events / bin", 0.60, 1.2);



}