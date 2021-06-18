
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <stdlib.h>
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct Output{
  double reco_xi_jets;
  double true_xi_jets;
  double matched_reco_xi_jets;
};

//------------------------------------------------------------------------------

Output *AnalyseEvents(const char *inputFile)
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
  TClonesArray *branchXiJet1  = treeReader->UseBranch("XiJet1");
  TClonesArray *branchXiJet2  = treeReader->UseBranch("XiJet2");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  //TClonesArray *branchXiJet1Size = treeReader->UseBranch("");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Jet *jet;
  TObject *object;
  GenParticle *GenJet_Const1;
  GenParticle *GenJet_Const2;
  double true_xi_jets = 0;
  double reco_xi_jets = 0;
  double matched_reco_xi_jets = 0;
  Long64_t entry;
  double truth_deltaR;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    //if(branchXiJet1->GetEntriesFast() == 0) continue;
    for(int k = 0; k < branchXiJet1->GetEntriesFast(); ++k){
       reco_xi_jets+=1;
     }

    // Loop over all genjets in event
    for(int i = 0; i < branchGenJet->GetEntriesFast(); ++i)
    {
      truth_deltaR = 0.0;

      jet = (Jet*) branchGenJet->At(i);
      if(jet->Constituents.GetEntriesFast() != 2) continue;

      GenJet_Const1 = (GenParticle*) jet->Constituents.At(0);
      GenJet_Const2 = (GenParticle*) jet->Constituents.At(1);

      if(jet->Constituents.At(0)->PID != 22 || jet->Constituents.At(1)->PID != 22) continue;
      if(jet->Constituents.At(0)->PT < 0.5 || jet->Constituents.At(1)->PT < 0.5 ) continue;

      truth_deltaR = sqrt(pow(GenJet_Const1->Eta - GenJet_Const2->Eta,2) + pow(GenJet_Const1->Phi - GenJet_Const2->Phi,2));

      genjet_const_R->Fill(truth_deltaR);

      //if(truth_deltaR < 0.025) continue; //Won't be able to see xi jets smaller than this

      true_xi_jets+=1;

      double genjet_eta = jet->Eta;
      double genjet_phi = jet->Phi;

      if(branchXiJet1->GetEntriesFast() == 0) continue;
      Jet *closest_jet;
      double closest_jet_R = 100000;
      //cout << "Number of Xi Jets in this event = " << branchXiJet1->GetEntriesFast() << endl;
      for(int k = 0; k < branchXiJet1->GetEntriesFast(); ++k){

        Jet *xi = (Jet*) branchXiJet1->At(k);
        double xi_eta = xi->Eta;
        double xi_phi = xi->Phi;

        double R = sqrt(pow(genjet_eta - xi_eta,2) + pow(genjet_phi - xi_phi,2));
        if(R < closest_jet_R) {
          closest_jet_R = R;
          closest_jet = (Jet*) xi;
        }

      }
      if(closest_jet_R > 0.05) continue;

      if(closest_jet->Tau[0] != 0){
        //std::cout << "Tau[2] = " << closest_jet->Tau[1] << std::endl;
        //std::cout << "Tau[1] = " << closest_jet->Tau[0] << std::endl;
  jettiness->Fill(closest_jet->Tau[1]/closest_jet->Tau[0]);
      }

      matched_reco_xi_jets+=1;

      true_reco_R->Fill(closest_jet_R);
      /*
      std::cout << "    We found a match for this True Xi Jet!" << std::endl;
      std::cout << "    DeltaR between GenJet and Reco Xi Jet = " << closest_jet_R << std::endl;
      std::cout << "    GenJet Pt = " << jet->PT << ", XiJet Pt = " << closest_jet->PT << std::endl;
      */

    }//loop over jets

  }//loop over events

  //double efficiency = reco_xi_jets/true_xi_jets;

  out->reco_xi_jets = reco_xi_jets;
  out->true_xi_jets = true_xi_jets;
  out->matched_reco_xi_jets = matched_reco_xi_jets;

  return out;
}//end of code

//------------------------------------------------------------------------------

void XiJets(const string filename)
{
  gSystem->Load("libDelphes");

  // Show resulting histograms
  int colors[8] = {1,2,3,4,6,8,40,47};

  TFile *outfile = new TFile(("/home/nastein/XiJets/Histograms/hists_" + filename).c_str(), "NEW");
 // TFile *outfile = new TFile(("test.root"), "NEW");


  std::cout << "Analyzing m = " << filename.substr(8,4) << " GeV." << std::endl;
  Output *out = AnalyseEvents(("/scratch/physdoe_project_root/physdoe_project/nastein/Benchmark_points/root_files/" + filename).c_str());

  out->const_deltaR->Write();
  out->true_recoR->Write();
  out->jettiness->Write();

  TVectorD v(4);
  v[0] = atof((filename.substr(8,4).c_str()));
  v[1] = out->true_xi_jets;
  v[2] = out->reco_xi_jets;
  v[3] = out->matched_reco_xi_jets;

  v.Write("Info");

}
