
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif


//------------------------------------------------------------------------------

struct Output{
  TH1D *const_deltaR;
  TH1D *true_recoR;
  double efficiency;
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

  TH1D* photon_const_R = new TH1D("photon_const_R", "photon_const_R", 40, 0.00, 0.05);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Photon *photon;
  GenParticle *Phot_Const1;
  GenParticle *Phot_Const2;
  double true_xi_jets = 0;
  double reco_xi_jets = 0;
  Long64_t entry;
  double truth_deltaR;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    //cout << "-- New Event" << endl;

    // Loop over all photons in event
    for(int i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      truth_deltaR = 0.0;

      photon = (Photon*) branchPhoton->At(i);
      if(photon->Particles.GetEntriesFast() != 2) continue;
      cout << "    -- New Photon" << endl;

      double deltaR = 0.0;
      double eta1, eta2;
      double phi1, phi2;

      for (int p = 0; p < photon->Particles.GetEntriesFast(); ++p) {
      Phot_Const1 = (GenParticle*) photon->Particles.At(p);
          //particle = (GenParticle*) object;
          cout << "        GenPart PID: " << Phot_Const1->PID << ", pt: " << Phot_Const1->PT << ", eta: " << Phot_Const1->Eta << ", phi: " << Phot_Const1->Phi <<", Status: " << Phot_Const1->Status << endl;
          if(p == 0) {
            eta1 = Phot_Const1->Eta; phi1 = Phot_Const1->Phi;
          }
          else {
            eta2 = Phot_Const1->Eta; phi2 = Phot_Const1->Phi;
          }
      }

      deltaR = sqrt(pow(eta1 - eta2,2) + pow(phi1 - phi2,2));
      std::cout << "        Delta R = " << deltaR << std::endl;


    }//loop over jets

  }//loop over events

  out->const_deltaR = photon_const_R;

  return out;
}//end of code

//------------------------------------------------------------------------------

void Photons(const char *input)
{
  gSystem->Load("libDelphes");

  Output *out = AnalyseEvents(input);
  TFile *outfile = new TFile("hists.root", "NEW");
  out->const_deltaR->Write();


}
