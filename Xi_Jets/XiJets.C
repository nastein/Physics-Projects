
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <stdlib.h>
#include <TMath.h>
#include "/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/madgraph_results/PlotUtils.C"

#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct Output{
  TH1D *const_deltaR;
  TH1D *true_recoR;
  TH1D *jettiness;
  TH1D *hadfrac;
  TH1D *numtracks;
  double reco_xi_jets;
  double true_xi_jets;
  double matched_reco_xi_jets;
};

//------------------------------------------------------------------------------

Output *AnalyseEvents(const char *inputFile, bool signal = false)
{

  Output *out = new Output;

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet1");
  TClonesArray *branchXiJet  = treeReader->UseBranch("XiJet2");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  TH1D* true_reco_R = new TH1D("True_reco_R", "True_reco_R", 20, 0.0, 0.05);
  TH1D* genjet_const_R = new TH1D("genjet_const_R", "genjet_const_R", 40, 0.0, 0.4);
  TH1D* jettiness = new TH1D("tau_{2/1}", "tau_{2/1}", 40, 0.0, 1.0);
  TH1D* hadfrac = new TH1D("hadfrac", "hadfrac", 20, -3.0, 0.0);
  TH1D* numtracks = new TH1D("numtracks", "numtracks", 3, 0, 3.0);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Jet *jet;
  TObject *object;
  GenParticle *GenJet_Const1;
  GenParticle *GenJet_Const2;
  GenParticle *GenJet_Const;
  double true_xi_jets = 0;
  double reco_xi_jets = 0;
  double matched_reco_xi_jets = 0;
  Long64_t entry;
  double truth_deltaR;

  bool good_genjet;
  int num_genphotons;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(signal == true) {
    // Loop over all genjets in event
      for(int i = 0; i < branchGenJet->GetEntriesFast(); ++i)
      {
        truth_deltaR = 0.0;
        good_genjet = true;
        num_genphotons = 0;

        jet = (Jet*) branchGenJet->At(i);

        if(jet->Constituents.GetEntriesFast() < 2) continue;
        //Defintion of our signal is a Jet with only two photons with PT > 0.5 GeV, and no non photons with PT > 0.5 GeV

        for(int c = 0; c < jet->Constituents.GetEntriesFast(); c++) {
          GenJet_Const = (GenParticle*) jet->Constituents.At(c);
          if (GenJet_Const->PID != 22 && GenJet_Const->PT > 0.5) {
            good_genjet = false; break;
          }
          if (GenJet_Const->PID == 22 && GenJet_Const->PT > 0.5) {
            if (num_genphotons == 2){
              good_genjet = false; break;
            }
            num_genphotons++;
            if (num_genphotons == 1 ) GenJet_Const1 = (GenParticle*) jet->Constituents.At(c);
            if (num_genphotons == 2 ) GenJet_Const2 = (GenParticle*) jet->Constituents.At(c);
          }
        }

        if (num_genphotons != 2) good_genjet = false;
        if (good_genjet == false) continue;

        truth_deltaR = sqrt(pow(GenJet_Const1->Eta - GenJet_Const2->Eta,2) + pow(GenJet_Const1->Phi - GenJet_Const2->Phi,2));
        genjet_const_R->Fill(truth_deltaR);
        if(truth_deltaR < 0.025) continue; //Won't be able to see xi jets smaller than this

        double genjet_eta = jet->Eta;
        double genjet_phi = jet->Phi;

        if(branchXiJet->GetEntriesFast() == 0) continue;
        Jet *closest_jet;
        double closest_jet_R = 100000;

        for(int k = 0; k < branchXiJet->GetEntriesFast(); ++k){

          Jet *xi = (Jet*) branchXiJet->At(k);
          double xi_eta = xi->Eta;
          double xi_phi = xi->Phi;

          double R = sqrt(pow(genjet_eta - xi_eta,2) + pow(genjet_phi - xi_phi,2));
          if(R < closest_jet_R) {
            closest_jet_R = R;
            closest_jet = (Jet*) xi;
          }
        }
        if(closest_jet_R > 0.05) continue;

        double number_of_tracks = 0;

        //Let's look at the constituents of the matched (closest) jet
        for(int j = 0; j < closest_jet->Constituents.GetEntriesFast(); j++) {
            object = closest_jet->Constituents.At(j);
            if(object == 0) continue;
            if(object->IsA() == Track::Class()) {
            number_of_tracks++;
            }
          }

        //Now we want to check the hadronic activity in the jet
        double Ehadfrac = closest_jet->Ehad/(closest_jet->Eem + closest_jet->Ehad);
        hadfrac->Fill(Ehadfrac);
        numtracks->Fill(number_of_tracks);

        if(closest_jet->Tau[0] != 0){
          double tau21 = closest_jet->Tau[1]/closest_jet->Tau[0];
          if(tau21 != 0)
            jettiness->Fill(closest_jet->Tau[1]/closest_jet->Tau[0]);
        }

      }//end of loop over gen jets

    }//end of if signal

    else {

      for(int i = 0; i < branchGenJet->GetEntriesFast(); ++i)
      {
        truth_deltaR = 0.0;
        good_genjet = true;
        num_genphotons = 0;

        jet = (Jet*) branchGenJet->At(i);
        //std::cout << "New jet" << std::endl;

        //if(jet->Constituents.GetEntriesFast() < 2) continue;
        //Defintion of our signal is a Jet with only two photons with PT > 0.5 GeV, and no non photons with PT > 0.5 GeV

        for(int c = 0; c < jet->Constituents.GetEntriesFast(); c++) {

          GenJet_Const = (GenParticle*) jet->Constituents.At(c);
          //std::cout << "PID = " << GenJet_Const->PID << std::endl;
          if (GenJet_Const->PID != 22 && GenJet_Const->PT > 0.5) {
            good_genjet = false; break;
          }
          if (GenJet_Const->PID == 22 && GenJet_Const->PT > 0.5) {
            if (num_genphotons == 2){
              good_genjet = false; break;
            }
            num_genphotons++;
            if (num_genphotons == 1 ) GenJet_Const1 = (GenParticle*) jet->Constituents.At(c);
            if (num_genphotons == 2 ) GenJet_Const2 = (GenParticle*) jet->Constituents.At(c);
          }
        }
      }


      for(int k = 0; k < branchXiJet->GetEntriesFast(); ++k){
        Jet *closest_jet = (Jet*) branchXiJet->At(k);

        double number_of_tracks = 0;

        //Let's look at the constituents of the matched (closest) jet
        for(int j = 0; j < closest_jet->Constituents.GetEntriesFast(); j++) {
            object = closest_jet->Constituents.At(j);
            if(object == 0) continue;
            if(object->IsA() == Track::Class()) {
            number_of_tracks++;
            }
          }

        //Now we want to check the hadronic activity in the jet
        double Ehadfrac = closest_jet->Ehad/(closest_jet->Eem + closest_jet->Ehad);
        hadfrac->Fill(Log10(Ehadfrac));
        numtracks->Fill(number_of_tracks);

        if(closest_jet->Tau[0] != 0){
          double tau21 = closest_jet->Tau[1]/closest_jet->Tau[0];
          if(tau21 != 0)
            jettiness->Fill(closest_jet->Tau[1]/closest_jet->Tau[0]);
        }
      }
    }

    matched_reco_xi_jets+=1;

    //true_reco_R->Fill(closest_jet_R);

  }//loop over events

  out->jettiness = jettiness;
  out->hadfrac = hadfrac;
  out->numtracks = numtracks;

  return out;
}//end of code

//------------------------------------------------------------------------------

void XiJets(const string filename)
{
  gSystem->Load("libDelphes");

  // Show resulting histograms
  int colors[8] = {1,2,3,4,6,8,40,47};

  Output *signal = AnalyseEvents("/Users/noahsteinberg/Physics/James_Research/Xi_Jets/runMass_1.19181765377272094220.lhe.gz.hepmc.root", true);
  Output *bkgd = AnalyseEvents("/Users/noahsteinberg/Physics/James_Research/Xi_Jets/gamma_j_bkgd.root");

  Plot plot("");
  plot.Add("signal", 3, 1, signal->hadfrac, true);
  plot.Add("jet + photon", 4, 1, bkgd->hadfrac, true);
  plot.SetLog();
  plot.Draw("Hadronic energy fraction", "Number of Xi Jets", 0.60, 1.2);

  Plot plot2("");
  plot2.Add("signal", 3, 1, signal->numtracks, true);
  plot2.Add("jet + photon", 4, 1, bkgd->numtracks, true);
  plot2.Draw("Number of tracks", "Number of Xi Jets", 0.60, 1.2);

  Plot plot3("");
  plot3.Add("signal", 3, 1, signal->jettiness, true);
  plot3.Add("jet + photon", 4, 1, bkgd->jettiness, true);
  plot3.Draw("tau2/tau1", "Number of Xi Jets", 0.60, 1.2);

}
