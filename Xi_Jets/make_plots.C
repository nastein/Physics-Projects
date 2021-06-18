#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TVectorD.h>
#include <TGraph.h>
#include "/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/madgraph_results/PlotUtils.C"

using namespace std;

void make_plots() {

  Plot plot("const_delta_R");
  plot.Title("DeltaR between constituents of reconstructed photons");

  Plot plot2("true_reco_R");
  plot2.Title("DeltaR generated jet and reconstructed xi jet");

  Plot plot3("matched reco");

  Plot plot4("true");

  Plot plot5("reco");

  Plot plot6("efficiency");
  plot6.Title("Efficiency of XiJet reconstruction");

  TGraph *true_jets = new TGraph();
  TGraph *reco_jets = new TGraph();
  TGraph *matched_reco_jets = new TGraph();
  TGraph *eff = new TGraph();
  TVectorD *info;
  double mass = 0;


  //std::ifstream file("hists1/files.txt");
  std::ifstream file("hists3/files3.txt");
  std::string str;
  int n = 0;
  while (std::getline(file,str)) {
    TFile *f = new TFile(("hists3/" + str).c_str());
    info = (TVectorD*)f->Get("Info");
    //plot.Add(str.substr(14,4), n, 1, (TH1D*)f->Get("genjet_const_R"));
    //plot2.Add(str.substr(14,4), n, 1, (TH1D*)f->Get("True_reco_R"));
    true_jets->SetPoint(n, atof(str.substr(14,4).c_str()), (*info)[1]);
    reco_jets->SetPoint(n, atof(str.substr(14,4).c_str()), (*info)[2]);
    matched_reco_jets->SetPoint(n, atof(str.substr(14,4).c_str()), (*info)[3]);
    eff->SetPoint(n, atof(str.substr(14,4).c_str()), (*info)[3]/(*info)[1]);
    n++;

  }

  //plot.SetLog(true);
  //plot.Draw("DeltaR", "Number of Jets", 0.60, 1.2);

  //plot2.Draw("DeltaR", "Number of Jets", 0.60, 1.2);

  plot3.Add("matched reco xi jets", 4, matched_reco_jets);
  plot3.DrawGraph("Scalar Mass (GeV)", "number of matched reco xi jets", 0.60);

  plot4.Add("true xi jets", 4, true_jets);
  plot4.DrawGraph("Scalar Mass (GeV)", "number of true xi jets", 0.60);

  plot5.Add("reco xi jets", 4, reco_jets);
  plot5.DrawGraph("Scalar Mass (GeV)", "number of reco xi jets", 0.60);

  plot6.Add("Efficiency", 4, eff);
  plot6.DrawGraph("Scalar Mass (GeV)", "Efficiency", 0.60);



}
