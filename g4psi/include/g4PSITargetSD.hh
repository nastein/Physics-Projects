#ifndef g4PSITargetSD_h
#define g4PSITargetSD_h 1

#include "g4PSITargetHit.hh"
#include "TTree.h"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class g4PSITargetSD : public G4VSensitiveDetector
{
    
public:
    g4PSITargetSD(G4String name,
                   G4String colName,
                   G4double min_p,
                   G4double min_theta);
    
    ~g4PSITargetSD();
    
    /// Initialize() method is invoked at the beginning of each event.
    void Initialize(G4HCofThisEvent*HCE);
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    void EndOfEvent(G4HCofThisEvent*HCE);
    void clear();
    void DrawAll();
    void PrintAll();
    void InitTree(TTree *T, std::string suffix = "");
    void DeleteEventData();
    G4bool SDHasHit() {return hit_;};
    
private:
    g4PSITargetHitsCollection *trackerCollection_;
    G4int HCID_;
    G4bool hit_;
    
//    G4double p_threshold_;
//    G4bool isTestPlane_;
    G4double min_p_;
    G4double min_theta_;
    
    typedef std::vector<G4int> VectorInt;
    typedef std::vector<G4double> VectorDouble;
    typedef std::vector<G4bool> VectorBool;
    
    G4int* TargetHit_;
    VectorInt* TargetCopyID_;
    VectorInt* TargetParticleID_;
    VectorInt* TargetTrackID_;
    VectorInt* TargetProcessID_;
    VectorDouble* Targetw_;
    VectorDouble* TargetVertexX_;
    VectorDouble* TargetVertexY_;
    VectorDouble* TargetVertexZ_;
    VectorDouble* TargetTheta_;
    VectorDouble* MultipleScatterTheta_;
    VectorDouble* TargetPIn_;
    VectorDouble* TargetPOut_;
    VectorBool* TargetInScatter_;

};

#endif

