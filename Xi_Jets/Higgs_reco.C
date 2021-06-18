
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

std::vector<std::vector<std::pair<TLorentzVector, std::string>>> my_pow_set;
//------------------------------------------------------------------------------

struct Output{
  int num_xi_matches;
  int num_photon_matches;
  int num_xi_with_phot_matches;
  int num_only_xi_jets_matches;
  int num_other_comb_matches;
  int num_matches;
};

void findPowerSet(std::vector<std::pair<TLorentzVector, std::string>> cands, std::vector<std::pair<TLorentzVector, std::string>>& power_cands, int n)
{
  // if we have considered all elements
  if (n == 0)
  {
    my_pow_set.push_back(power_cands);
    return;
  }
  // consider nth element
  power_cands.push_back(cands[n - 1]);
  findPowerSet(cands, power_cands, n - 1);

  // or don't consider nth element
  power_cands.pop_back();
  findPowerSet(cands, power_cands, n - 1);
}

bool pass_higgs_mass(double minv, double mhiggs, double sigma) {
  return ((minv < mhiggs + sigma) && (minv > mhiggs - sigma));
}

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
  TClonesArray *branchXiJet  = treeReader->UseBranch("XiJet2");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  TH1D *Mass = new TH1D("Mass", "Mass", 15, 110, 140);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Jet *jet;
  Photon *photon;
  TObject *object;
  Long64_t entry;
  std::vector<std::pair<TLorentzVector,std::string>> momentums;
  std::vector<std::pair<TLorentzVector,std::string>> power_momentums;
  int num_matches = 0;
  int num_xi_matches = 0;
  int num_photon_matches = 0;
  int num_xi_with_phot_matches = 0;
  int num_only_xi_jets_matches = 0;
  int num_other_comb_matches = 0;

  int number_of_tracks;

  TLorentzVector mom_sum;

  bool a_weird_match;

  TLorentzVector best_mom_sum;
  std::vector<std::pair<TLorentzVector, std::string>> best_objects;
  double mass_diff;

  double mhiggs = 125.10;
  double sigma = 3.0;

  int num_xi_in_event;
  int num_phots_in_event;

  bool higgs_match;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {

    std::cout << "New Event!" << std::endl;

    higgs_match = false;

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);



    for(int k = 0; k < branchXiJet->GetEntriesFast(); ++k){
      number_of_tracks = 0;
      jet = (Jet*)branchXiJet->At(k);
      //Look at hadronic energy fraction
      if (log(jet->Ehad/(jet->Ehad + jet->Eem)) > -0.8) continue;

      //Look at number of tracks
      for(int j = 0; j < jet->Constituents.GetEntriesFast(); j++) {
        object = jet->Constituents.At(j);
        if(object == 0) continue;
        if(object->IsA() == Track::Class()) {
          number_of_tracks++;
        }
      }

      if(number_of_tracks > 0) continue;

      momentums.push_back(std::make_pair(jet->P4(),std::string("xi-jet")));
    }


    for(int k = 0; k < branchPhoton->GetEntriesFast(); ++k){
      photon = (Photon*)branchPhoton->At(k);
      momentums.push_back(std::make_pair(photon->P4(),"photon"));
    }

    int size = momentums.size();
    findPowerSet(momentums, power_momentums, size);

    best_mom_sum.SetPxPyPzE(0.0,0.0,0.0,0.0);
    mass_diff = 99999999999;

    //loop over the power set (each objs is a set of jets)
    for(auto objs = begin (my_pow_set); objs != end(my_pow_set); ++objs) {
      //so objs is a set of jets
      mom_sum.SetPxPyPzE(0.0,0.0,0.0,0.0);
      num_xi_in_event = 0;
      num_phots_in_event = 0;

      if(objs->size() > 4 || objs->size() == 0) continue;

      //Sum the four momenta of our objects
      for(int i = 0; i < objs->size(); i++) {
        mom_sum = mom_sum + (*objs)[i].first;
      }

      //If we fall into the higgs mass region, check if this is the best match, if so update best_object and best_mom_sum
      if(pass_higgs_mass(mom_sum.M(), mhiggs, sigma)) {
        higgs_match = true;
        num_matches++;
        a_weird_match = true;
        for(int i = 0; i < objs->size(); i++) {
          if ((*objs)[i].second == "xi-jet") num_xi_in_event++;
          else num_phots_in_event++;
        }
        if(
          (num_xi_in_event == 1 && num_phots_in_event == 1) ||
          (num_xi_in_event == 2 && num_phots_in_event == 1) ||
          (num_xi_in_event == 1 && num_phots_in_event == 2)
        ) {num_xi_with_phot_matches++; a_weird_match = false;}
        if(
          (num_xi_in_event == 2 && num_phots_in_event == 0)
        ) {num_only_xi_jets_matches++; a_weird_match = false;}
        if(
          (num_xi_in_event == 1 && num_phots_in_event == 1) ||
          (num_xi_in_event == 2 && num_phots_in_event == 1) ||
          (num_xi_in_event == 1 && num_phots_in_event == 2) ||
          (num_xi_in_event == 2 && num_phots_in_event == 0)
        ) {num_xi_matches++; a_weird_match = false;}
        if(
          (num_phots_in_event == 2 && num_xi_in_event == 0) ||
          (num_phots_in_event == 3 && num_xi_in_event == 0) ||
          (num_phots_in_event == 4 && num_xi_in_event == 0)
        ) {num_photon_matches++; a_weird_match = false;}
        if(a_weird_match == true)
          {num_other_comb_matches++;}
      }
    }


    momentums.clear();
    power_momentums.clear();
    my_pow_set.clear();

  }//loop over events

  out->num_matches = num_matches;
  out->num_xi_matches = num_xi_matches;
  out->num_xi_with_phot_matches = num_xi_with_phot_matches;
  out->num_photon_matches = num_photon_matches;
  out->num_only_xi_jets_matches = num_only_xi_jets_matches;
  out->num_other_comb_matches = num_other_comb_matches;

  return out;
}//end of code

//------------------------------------------------------------------------------

void Higgs_reco(const string filename)
{
  gSystem->Load("libDelphes");

  // Show resulting histograms
  int colors[8] = {1,2,3,4,6,8,40,47};

  TFile *outfile = new TFile(("/home/nastein/XiJets/HiggsReco/higgs_" + filename).c_str(), "NEW");

  Output *out = AnalyseEvents(("/scratch/physdoe_project_root/physdoe_project/nastein/Benchmark_points/root_files2/" + filename).c_str());

  TVectorD v(6);
  v[0] = out->num_matches;
  v[1] = out->num_xi_matches;
  v[2] = out->num_xi_with_phot_matches;
  v[3] = out->num_only_xi_jets_matches;
  v[4] = out->num_photon_matches;
  v[5] = out->num_other_comb_matches;

  v.Write("Info");

}


