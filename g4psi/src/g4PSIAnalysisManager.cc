#include "g4PSIAnalysisManager.hh"
#include "g4PSIDetectorParts.hh"
#include "g4PSIDetectorBase.hh"
//#include "g4PSIRunAction.hh"
#include "g4PSIRun.hh"

#include "G4UnitsTable.hh"
#include "G4SDManager.hh"

#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>

using namespace CLHEP;

g4PSIAnalysisManager* g4PSIAnalysisManager::fManager = 0;

g4PSIAnalysisManager* g4PSIAnalysisManager::getInstance() {
    if(!fManager) {
        fManager = new g4PSIAnalysisManager();
    }
    return fManager;
}


void g4PSIAnalysisManager::dispose() {
    delete fManager;
    fManager = NULL;
}

g4PSIAnalysisManager::g4PSIAnalysisManager() {
    
    // Resolution / smearing for scintillator and sci-fi detector
    scwall_tube_res = 0.0*ns;
    beamceren_tube_res = 0.0*ns;
    sf_plane_res = 0.0*ns;
    
    physics_bias_do_it_ = false;
    physics_bias_scale_ = 1.0;
    
    tree_ = NULL;
}

g4PSIAnalysisManager::~g4PSIAnalysisManager()
{
}

// =========================================================================

void g4PSIAnalysisManager::BeginOfRun(G4String root_file_name)
{

        event_id_ = 0;
        event_seed1_ = -1;
        event_seed2_ = -1;
        
        G4cout << "g4PSIAnalysisManager: Run has been started. Root file: " << root_file_name << G4endl;
        root_file_ = new TFile(root_file_name, "RECREATE", "MC Simulation for PSI Muon-Scattering Experiment");
        
        // prepare the root tree
        tree_ = new TTree("T", "MC Events");
        tree_->SetAutoSave();
        tree_->Branch("EventID", &event_id_, "EventID/I");
        tree_->Branch("EventSeed1", &event_seed1_, "EventSeed1/L");
        tree_->Branch("EventSeed2", &event_seed2_, "EventSeed2/L");
        
        g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
        DetectorParts->RootInitTree(tree_);
        
        tree_->Branch("MuonDecay", &MuonDecay_);
        tree_->Branch("MuonDecay_W", &MuonDecayW_);
        tree_->Branch("MuonDecay_PosX", &MuonDecayPosX_);
        tree_->Branch("MuonDecay_PosY", &MuonDecayPosY_);
        tree_->Branch("MuonDecay_PosZ", &MuonDecayPosZ_);
        tree_->Branch("MuonDecay_Mom", &MuonDecayMom_);
        tree_->Branch("MuonDecay_DirX", &MuonDecayDirX_);
        tree_->Branch("MuonDecay_DirY", &MuonDecayDirY_);
        tree_->Branch("MuonDecay_DirZ", &MuonDecayDirZ_);
        
    //    tree_->Branch("Target_event", &TargetEvent_);
    //    tree_->Branch("Target_w", &TargetW_);
    //    tree_->Branch("Target_index", &TargetIndex_);
    //    tree_->Branch("Target_vertex", &TargetVertex_);
    //    tree_->Branch("Target_p_in", &TargetMomentumIn_);
    //    tree_->Branch("Target_p_out", &TargetMomentumOut_);
        
        tree_->Branch("rf_time", &RF_time_);
        tree_->Branch("beam_particles", &beam_particles_);
        tree_->Branch("beam_id", &beam_id_);
        tree_->Branch("beam_p", &beam_p_);
        tree_->Branch("beam_theta", &beam_theta_);
        tree_->Branch("beam_phi", &beam_phi_);
        tree_->Branch("beam_x0", &beam_x0_);
        tree_->Branch("beam_y0", &beam_y0_);
        tree_->Branch("beam_z0", &beam_z0_);


}

// =========================================================================

void g4PSIAnalysisManager::EndOfRun(const G4Run* aRun) {

    // store in the output file information about the run and physics bias
    TVectorD v1(3);
    v1[0] = aRun->GetNumberOfEvent();
    v1[1] = physics_bias_do_it_ ? 1.0 : 0.0;
    v1[2] = physics_bias_scale_;
    v1.Write("RunControl");
    
    g4PSIDetectorParts::getInstance()->RootWrite();
    
    root_file_->Write();
    root_file_->Close();
    delete root_file_;

}

// =========================================================================

void g4PSIAnalysisManager::BeginOfEvent(const G4Event* evt)
{
    // Detector variables will be reset by the respective sensitive detectors
    // Here we need to reset only analysis-manager specific & event-based variables
    //
    // Beam variables will be reset by the primary action class, this BeginOfEvent
    // method will be called only after the primary particles are set.
    


    (void) evt;
    
    MuonDecay_ = false;
    MuonDecayW_ = 0.0;
    MuonDecayPosX_ = 0.0;
    MuonDecayPosY_ = 0.0;
    MuonDecayPosZ_ = 0.0;
    MuonDecayMom_ = 0.0;
    MuonDecayDirX_ = 0.0;
    MuonDecayDirY_ = 0.0;
    MuonDecayDirZ_ = 0.0;

    is_signal_ = false;

}

// =========================================================================

void g4PSIAnalysisManager::EndOfEvent(const G4Event* evt)
{
    G4int evtNb = evt->GetEventID();
    if (evtNb % 10000 == 0)  G4cout << "*=== Event number = " << evtNb << " ===*" << G4endl;
    
    if (g4PSIDetectorParts::getInstance()->HasTrigger()) {
        tree_->Fill();
    }
    
    /*
    if (is_signal_ == true) {
        tree_->Fill();
    }
    */
    
    
}

// =========================================================================

void g4PSIAnalysisManager::SetParticleGun(G4ParticleGun *pgun)
{
    G4ParticleDefinition* particle = pgun->GetParticleDefinition();
    G4double m0 = particle->GetPDGMass();
    G4double e = pgun->GetParticleEnergy() + m0;

    // as the kinetic energy of the particle is set in the gun, the momentum
    // is put to zero in G4ParticleGun.cc, so we need to read the energy
    // and calculate the momentum.
    
    beam_particles_++;
    beam_id_.push_back(particle->GetPDGEncoding());
    beam_theta_.push_back((pgun->GetParticleMomentumDirection()).theta());
    beam_phi_.push_back((pgun->GetParticleMomentumDirection()).phi());
    beam_p_.push_back(sqrt(e*e - m0*m0));
    beam_x0_.push_back((pgun->GetParticlePosition().x()));
    beam_y0_.push_back((pgun->GetParticlePosition().y()));
    beam_z0_.push_back((pgun->GetParticlePosition().z()));
}

void g4PSIAnalysisManager::ClearParticleGun()
{
    beam_particles_ = 0;
    beam_id_.clear();
    beam_theta_.clear();
    beam_phi_.clear();
    beam_p_.clear();
    beam_x0_.clear();
    beam_y0_.clear();
    beam_z0_.clear();
}