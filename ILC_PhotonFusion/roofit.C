#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "RooPlot.h"
#include "TAxis.h"

#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooAbsArg.h"
#include "RooWorkspace.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include <string>

#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"

#include "RooMinimizer.h"

using namespace RooFit;
using namespace RooStats;

 
double MakeModel(RooWorkspace *wks, TNtupleD* results, double alp_mass)
{

   //Files containing simulated sig and simualted bkgd
   TFile *s500 = new TFile("/Users/noah/Physics/James_Research/ILC_PhotonFusion/500GeV_signal.root");
   TFile *LBL500 = new TFile("/Users/noah/Physics/James_Research/ILC_PhotonFusion/500GeV_LBL.root");

   RooRealVar lumi_nom("lumi_nom", "lumi_nom", 3000000, 0, 10000000);
   lumi_nom.setConstant(kTRUE);

   //POI (signal cross section)
   RooRealVar xsec("xsec", "xsec", 0.1, 0.0, 1.0);

   //Background cross section
   RooRealVar bkgd_xsec("bkgd_xsec", "bkgd_xsec", 0.093765120, 0.05, 0.2);

   RooRealVar sigeff("sigeff", "sigeff", 0.85, 0.01, 1);
   sigeff.setConstant(kTRUE);

   RooProduct sigyield("sigyield", "sigyield", RooArgSet(sigeff,  lumi_nom, xsec));
   RooProduct bkgdyield("bkgdyield", "bkgdyield", RooArgSet(bkgd_xsec, lumi_nom));

   //Observable
   RooRealVar M_inv("M_inv", "Invariant mass of gamma gamma pair", 0, 200, "GeV");
   M_inv.setBins(200);

   //Signal PDF
   RooRealVar gauss_mean("gauss_mean", "mean", alp_mass, alp_mass - 2.0, alp_mass + 2.0, "GeV");
   gauss_mean.setConstant(kTRUE);
   RooFormulaVar gauss_sigma("gauss_sigma", "sigma", "sqrt(pow(@0 * .01,2) + @0 * pow(.17,2))", RooArgList(gauss_mean));
   RooGaussian peak_gaussian("peak_gaussian", "gaussian", M_inv, gauss_mean, gauss_sigma);

   //Bkgd pdf
   RooRealVar cb_mean("cb_mean", "mean", 5.2477, 1, 25);
   RooRealVar cb_sigma("cb_sigma", "sigma", .7889, 0, 10);
   RooRealVar cb_alpha("cb_alpha", "alpha", -.3579, -15, -0.01);
   RooRealVar cb_n("cb_n", "n", 3.69195, 1, 4);
   RooCBShape CBall("CBall", "Crystal ball shape", M_inv, cb_mean, cb_sigma, cb_alpha, cb_n);

   //Get Histograms of signal and background distributions
   TH1D *s500_hist = (TH1D*) s500->Get("M_{#gamma#gamma}35.000000");
   TH1D *LBL500_hist = (TH1D*) LBL500->Get("M_{#gamma#gamma}bkgd35.000000");

   RooDataHist s500_dh("minv", "minv", M_inv, Import(*s500_hist));
   RooDataHist LBL500_dh("minv", "minv", M_inv, Import(*LBL500_hist));

   RooCategory sample("sample","");
   sample.defineType("signal");
   sample.defineType("bkgd");

   //two datahists, one containing sig + bkgd, one containing bkgd only 
   RooDataHist combData("combData", "combined data", M_inv, Index(sample), Import("signal",s500_dh), Import("bkgd",LBL500_dh));
   RooDataHist combData2("combData2", "combined data2", M_inv, Index(sample), Import("bkgd",LBL500_dh));

   
   //Construct combined sig + bkgd pdf
   RooArgList shapes;
   RooArgList yields;

   double nbg = combData2.sum(kTRUE);
   RooRealVar bkgd_yield("bkgd_yield", "bkgd yield", nbg, 0, 5000000);

   shapes.add(peak_gaussian); shapes.add(CBall);
   yields.add(sigyield); yields.add(bkgd_yield);
   RooAddPdf model("model","total", shapes, yields);

   //Import into workspace
   wks->import(model);
   wks->import(combData2,Rename("data"));

   //Do I even need this?
   RooArgSet global_obs("global_obs");

   RooArgSet nuis("nuis");
   nuis.add(*wks->var("cb_n"));
   nuis.add(*wks->var("cb_alpha"));
   nuis.add(*wks->var("cb_mean"));
   nuis.add(*wks->var("cb_sigma"));
   nuis.add(*wks->var("bkgd_yield"));

   RooArgSet poi("poi");
   poi.add(*wks->var("xsec"));

   RooArgSet obs("obs");
   obs.add(*wks->var("M_inv"));

   //Create bgkd model config
   ModelConfig bmodel("bmodel", wks);
   //bmodel.SetWS(*wks);
   bmodel.SetPdf(*wks->pdf("model"));
   bmodel.SetObservables(obs);
   bmodel.SetParametersOfInterest(poi);
   wks->var("xsec")->setVal(0.0);
   //bmodel.SetGlobalObservables(global_obs);
   //bmodel.SetNuisanceParameters(nuis);
   bmodel.SetSnapshot(poi);

   //Create sig + bkgd model config
   ModelConfig sbmodel("sbmodel", wks);
   sbmodel.SetWS(*wks);
   sbmodel.SetPdf(*wks->pdf("model"));
   sbmodel.SetObservables(obs);
   sbmodel.SetParametersOfInterest(poi);
   wks->var("xsec")->setVal(.2);// Does this value matter?
   //bmodel.SetGlobalObservables(global_obs);
   //bmodel.SetNuisanceParameters(nuis);
   sbmodel.SetSnapshot(poi);

   wks->import(bmodel);
   wks->import(sbmodel);
   wks->writeToFile("ALPworkspace.root");

   return 0;

}

void roofit(double starting_mass = 30)
{


   RooRandom::randomGenerator()->SetSeed(1);

   TNtupleD* results = new TNtupleD("results","results","mu_obs:mu_exp:m");

   int num_masses = 1;
   Double_t masses[num_masses];
   //Double_t upperlimit[50];
   for (int i = 0; i < num_masses; i++) {
      // Create a workspace to manage the project.
      RooWorkspace *wspace = new RooWorkspace("myWS");
      masses[i] = starting_mass + i*5;
      MakeModel(wspace, results, masses[i]);
      // cleanup
      delete wspace;
   }

   /*
   TCanvas* scanCan = new TCanvas("scanCan","scan",800,600);
   Int_t n1 = results->Draw("m:mu_exp:mu_obs","","goff");

    //Plot the expected and observed limits
    TGraph* muLimit_exp = new TGraph(n1,results->GetVal(0),results->GetVal(1));
    muLimit_exp->SetLineWidth(2);
    muLimit_exp->SetLineStyle(2);
    muLimit_exp->GetYaxis()->SetTitle("upper limits on cross section");
    muLimit_exp->GetXaxis()->SetTitle("ALP Mass (GeV)");
    TGraph* muLimit_obs = new TGraph(n1,results->GetVal(0),results->GetVal(2));    
    muLimit_obs->SetLineWidth(3);
    muLimit_obs->SetLineColor(kRed);
    muLimit_exp->Draw("APL");
    muLimit_obs->Draw("same");
    gPad->SetLogy();
    */



}
