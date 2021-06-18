#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TVectorD.h>
#include <TGraph.h>
#include "/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/madgraph_results/PlotUtils.C"

using namespace std;

void xsec() {

  //Plot plot("number of matches");
  //plot.Title("");

  TGraph *x = new TGraph();

  Double_t masses[12] = {2, 2.5, 3, 4, 4.5, 5, 5.5, 6, 7, 10, 15, 20};
  Double_t xsecs[12] = {.00611, .0119, .0206, .0487, .069, .0947, .1258, .1631, .2578, .7374, 2.3773, 5.2628};

  int n = 0;
  for (n; n < 12; n++) {
    x->SetPoint(n, masses[n], xsecs[n]/16);
  }

  x->SetMarkerStyle(21);
  x->SetMarkerSize(0.7);
  x->SetMarkerColor(1);
  x->SetLineColor(1);
  x->SetLineWidth(2);



  TLegend *legend=new TLegend(0.5,0.65,0.80,0.80);
  legend->SetTextFont(60);
  legend->SetTextSize(0.05);
  legend->SetFillColor(0);

  legend->AddEntry(x, "g = 1e-3", "p");
  x->GetYaxis()->SetTitle("#sigma(e+ e- #rightarrow Z #rightarrow #phi #gamma #rightarrow 3#gamma) (pb)");
  x->GetYaxis()->CenterTitle(-1);
  x->GetXaxis()->SetTitle("Mass (GeV)");
  x->GetXaxis()->CenterTitle(-1);
  x->Draw("AP");
  //x->GetXaxis()->SetLimits(1.5,21);
  //x->GetYaxis()->SetLimits(0.1, 200000);
  //bkgd->Draw();
  legend->Draw();

  gStyle->SetOptStat(0);


}