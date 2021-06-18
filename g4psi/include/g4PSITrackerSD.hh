#ifndef g4PSITrackerSD_h
#define g4PSITrackerSD_h 1

#include "g4PSITrackerHit.hh"
#include "TTree.h"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class g4PSITrackerSD : public G4VSensitiveDetector
{
    
public:
    g4PSITrackerSD(G4String name,
                   G4String colName,
                   G4bool isTestPlane = false,
                   G4double pth = 0 * CLHEP::MeV);  // if test_plane
    
    ~g4PSITrackerSD();
    
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
    void SetDepth(G4int d) {depth_ = d;};
    
private:
    g4PSITrackerHitsCollection *trackerCollection_;
    G4int HCID_;
    G4bool hit_;
    G4int depth_;
    
    G4double p_threshold_;
    G4bool isTestPlane_;
    
    typedef std::vector<G4int> VectorInt;
    typedef std::vector<G4double> VectorDouble;
    
    G4int* TrackHit_;
    VectorInt* TrackCopyID_;
    VectorInt* TrackParticleID_;
    VectorInt* TrackTrackID_;
    VectorDouble* Trackw_;
    VectorDouble* Trackp_;
    VectorDouble* TrackEdep_;
    VectorDouble* TrackPosHitX_;
    VectorDouble* TrackPosHitY_;
    VectorDouble* TrackPosHitZ_;
    VectorDouble* TrackDirHitX_;
    VectorDouble* TrackDirHitY_;
    VectorDouble* TrackDirHitZ_;
};

#endif

