#ifndef g4PSIScintillatorSD_h
#define g4PSIScintillatorSD_h 1

#include "g4PSIScintillatorHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "TTree.h"

/// \todo what is the meaning of the weight for this sensitive detector, when potentially the SD of the same copyID is hit by tracks of different weights?T

class g4PSIScintillatorSD : public G4VSensitiveDetector
{
    
public:
    g4PSIScintillatorSD(G4String name, G4int nCells, G4String colName);
    g4PSIScintillatorSD(G4String name, G4int nCells, G4String colname, bool collect_PID);
    ~g4PSIScintillatorSD();
    
    /// Initialize() method is invoked at the beginning of each event.
    void Initialize(G4HCofThisEvent*HCE);
    
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    void EndOfEvent(G4HCofThisEvent*HCE);
    void clear();
    void DrawAll();
    void PrintAll();
    G4int GetProcessCode(G4String process_name);
    void InitTree(TTree *T);
    void DeleteEventData();
    G4bool SDHasHit() {return hit_;};
    
private:
    g4PSIScintillatorHitsCollection *CalCollection;
    
    int* CellID;
    int numberOfCells;
    int HCID;
    bool hit_;
    bool collect_PID_ = false;
    
    typedef std::vector<G4int> VectorInt;
    typedef std::vector<G4double> VectorDouble;
    
    G4int* SCHit_;           // number of copies of this SD with hit.
    VectorInt* SCCopyID_;    // copyID of each hit copy.
    VectorDouble* SCW_;
    VectorInt* SCParticle_;
    VectorDouble* SCEdep_;
    VectorDouble* SCTime_;
    VectorDouble* SCPosX_;
    VectorDouble* SCPosY_;
    VectorDouble* SCPosZ_;
    VectorDouble* SCPosHitX_;
    VectorDouble* SCPosHitY_;
    VectorDouble* SCPosHitZ_;
    VectorDouble* SCDirHitX_;
    VectorDouble* SCDirHitY_;
    VectorDouble* SCDirHitZ_;
    VectorDouble* SCPosOriginX_;
    VectorDouble* SCPosOriginY_;
    VectorDouble* SCPosOriginZ_;
    VectorDouble* SCMuonEdep_;
    VectorDouble* SCElectronEdep_;
    VectorDouble* SCPhotonEdep_;
    //    VectorDouble* SCProcess_;
};

#endif

