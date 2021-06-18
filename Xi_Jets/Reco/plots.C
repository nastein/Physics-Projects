#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TVectorD.h>
#include <TGraph.h>
#include "/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/madgraph_results/PlotUtils.C"

using namespace std;

void plots() {

  Plot plot2("number of matches");
  plot2.Title("");


  TGraph *num_matches = new TGraph(); //all matches
  TGraph *num_matches_withonlyphotons = new TGraph(); //photons only
  TGraph *num_matches_withxijetsandphotons = new TGraph(); //photons + xi jets
  TGraph *num_matches_withxijets = new TGraph(); //xi jets inclusive
  TGraph *num_matches_onlyxijets = new TGraph(); //xi jets only
  TGraph *num_matches_other_combs = new TGraph(); //anything else
  TVectorD *info;
  //double mass = 0;

  int colors[19] = {1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20};

  std::ifstream file("3GeV_matches/files.txt");
  std::string str;
  int n = 1;
  while (std::getline(file,str)) {
    TFile *f = new TFile(("3GeV_matches/" + str).c_str());
    info = (TVectorD*)f->Get("Info");
    num_matches->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[0]/25000);
    num_matches_withxijets->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[1]/25000);
    num_matches_withxijetsandphotons->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[2]/25000);
    num_matches_onlyxijets->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[3]/25000);
    num_matches_withonlyphotons->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[4]/25000);
    num_matches_other_combs->SetPoint(n-1, atof(str.substr(14,4).c_str()), (*info)[5]/25000);
    n++;
}


  plot2.Add("photons only", 3, num_matches_withonlyphotons);
  plot2.Add("photons + xi jets", 2, num_matches_withxijetsandphotons);
  plot2.Add("xi jets only", 7, num_matches_onlyxijets);
  plot2.Add("xi jets inclusive", 6, num_matches_withxijets);
  plot2.Add("other combs", 9, num_matches_other_combs);
  plot2.Add("All", 4, num_matches);
  plot2.DrawGraphRatio("Scalar Mass (GeV)", "Number of matches", 0.60);

}
