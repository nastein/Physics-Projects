#ifndef g4PSIAnalysisManager_h
#define g4PSIAnalysisManager_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIAnalysisManager
//!
//! Description: Singleton class to hold analysis parameters and build histograms.
//!              User cannot access to the constructor.
//!              The pointer of the only existing object can be got via
//!              g4PSIAnalysisManager::GetInstance() static method.
//!              The first invokation of this static method makes
//!              the singleton object.
//!
//! Reference:   ../geant4.9.2/examples/extended/radioactivedecay/exrdm
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "G4Event.hh"
#include "g4PSIRun.hh"
#include "G4ParticleGun.hh"

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>

//#include "g4PSIEvent.hh"

class g4PSIAnalysisManager {
    
public:
    static g4PSIAnalysisManager* getInstance();
    static void dispose();
    
public:
    void BeginOfRun(G4String root_file);
    void EndOfRun(const G4Run* aRun);
    
    void BeginOfEvent(const G4Event* evt);
    void EndOfEvent(const G4Event* evt);
    
    TTree* GetTree() {return tree_;};
    
    Int_t GetEventID() {return event_id_;};
    Int_t GetEventSeed1() {return event_seed1_;};
    Int_t GetEventSeed2() {return event_seed2_;};
    void SetEventID(Int_t n) {event_id_ = n;};
    void SetEventID(Int_t n, Long_t s1, Long_t s2) {event_id_ = n; event_seed1_ = s1; event_seed2_ = s2;};
    void SetParticleGun(G4ParticleGun* pgun);
    void ClearParticleGun();
    void SetMuonDecay(Double_t w, Double_t x, Double_t y, Double_t z, Double_t p, Double_t a, Double_t b, Double_t c) {
        MuonDecay_ = true;
        MuonDecayW_ = w;
        MuonDecayPosX_ = x; MuonDecayPosY_ = y; MuonDecayPosZ_ = z;
        MuonDecayMom_ = p;
        MuonDecayDirX_ = a; MuonDecayDirY_ = b; MuonDecayDirZ_ = c;
    };
    bool is_signal(Bool_t t) {
        is_signal_ = t;
    };
//    void SetTargetProcess(Double_t w, TVector3 vertex, TVector3 mom_in, TVector3 mom_out, int index) {
//        TargetEvent_++;
//        TargetW_ = w;
//        TargetVertex_ = vertex;
//        TargetMomentumIn_ = mom_in;
//        TargetMomentumOut_ = mom_out;
//        TargetIndex_ = index;
//    };
//    Double_t GetTargetVertexZ() {return TargetVertex_.z();};
    
    //  void SetSCPath(Double_t p) {sc_path_ = p;};
    //  void SetInSCActiveRegion(Int_t x) {inActiveRegion_ = x;};
    
    void SetResSCWall(Double_t value) {scwall_tube_res = value;};
    void SetResBeamCeren(Double_t value) {beamceren_tube_res = value;};
    void SetResSFPlane(Double_t value) {sf_plane_res = value;};
    void SetBiasingPhysics(Bool_t bias, Double_t scale) {physics_bias_do_it_ = bias; physics_bias_scale_ = scale;};
    Double_t GetBiasingPhysicsScale() {return physics_bias_scale_;};
    
private:
    
    TFile *root_file_;
    TTree *tree_;
    
    Double_t scwall_tube_res;
    Double_t beamceren_tube_res;
    Double_t sf_plane_res;
    
    Int_t event_id_;
    Long_t event_seed1_;
    Long_t event_seed2_;
    
    Double_t RF_time_;

    Int_t beam_particles_;
    std::vector<Int_t> beam_id_;
    std::vector<Double_t> beam_theta_;
    std::vector<Double_t> beam_phi_;
    std::vector<Double_t> beam_p_;
    std::vector<Double_t> beam_x0_;
    std::vector<Double_t> beam_y0_;
    std::vector<Double_t> beam_z0_;
    
//    typedef std::vector<double> VectorDouble;
//    typedef std::vector<int> VectorInt;
    
    Bool_t MuonDecay_;
    Double_t MuonDecayW_;
    Double_t MuonDecayPosX_;
    Double_t MuonDecayPosY_;
    Double_t MuonDecayPosZ_;
    Double_t MuonDecayMom_;
    Double_t MuonDecayDirX_;
    Double_t MuonDecayDirY_;
    Double_t MuonDecayDirZ_;

    Bool_t is_signal_;
    
//    Int_t TargetEvent_;
//    Double_t TargetW_;
//    TVector3 TargetVertex_;
//    TVector3 TargetMomentumIn_;
//    TVector3 TargetMomentumOut_;
//    Int_t TargetIndex_;
    
    Bool_t physics_bias_do_it_;
    Double_t physics_bias_scale_;
    
private:
    static g4PSIAnalysisManager* fManager;
    g4PSIAnalysisManager();    // private constructor
    ~g4PSIAnalysisManager();
    g4PSIAnalysisManager(const g4PSIAnalysisManager&);   // Prevent copy-construction
    g4PSIAnalysisManager& operator=(const g4PSIAnalysisManager&);  // Prevent assignment}
};

#endif
