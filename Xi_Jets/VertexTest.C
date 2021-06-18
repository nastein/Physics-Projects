
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
};

//------------------------------------------------------------------------------

void VertexTest(const string filename)
{

  Output *out = new Output;

  TChain *chain = new TChain("Delphes");
  chain->Add(filename.c_str());

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchXiJet  = treeReader->UseBranch("XiJet2");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  TH1D *X = new TH1D("X", "X", 10, -100, 100);
  TH1D *Y = new TH1D("Y", "Y", 10, -100, 100);
  TH1D *Z = new TH1D("Z", "Z", 10, -100, 100);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Jet *jet;
  Photon *photon;
  Long64_t entry;
  GenParticle *part;
  TLorentzVector mom_sum;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {

    //std::cout << "New Event!" << std::endl;

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    for(int k = 0; k < branchParticle->GetEntriesFast(); ++k) {
      part = (GenParticle*)branchParticle->At(k);
      if(part->PID != 22) continue;
      //if (part->X < 1 || part->Y < 1) continue;
      X->Fill(part->X);
      Y->Fill(part->Y);
      Z->Fill(part->Z);
    }

  }//loop over events

  TCanvas *c1 = new TCanvas();
  gPad->SetLogy();
  X->Draw();

  TCanvas *c2 = new TCanvas();
  gPad->SetLogy();
  Y->Draw();

  TCanvas *c3 = new TCanvas();
  gPad->SetLogy();
  Z->Draw();

}//end of code

//----------------------------------------------------------------------------

