#ifndef g4PSIVeto_h
#define g4PSIVeto_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIVeto
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

class g4PSIVeto : public g4PSIDetectorBase {
    
public:
    g4PSIVeto(G4String label,
                G4double width, G4double z, G4int nSeg, G4double r1, G4double r2, G4int nPlane);
    g4PSIVeto(G4String label,
                G4double width, G4double z, G4int n, G4double r1, G4double r2);
    
    ~g4PSIVeto();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return sc_label_;};
    void Write();
    void Info();
    
    void InitTree(TTree *T);
    void DeleteEventData();
    void ReadEventData(G4HCofThisEvent *HCE);
    
private:
    G4String sc_label_;
    G4double sc_width_;
    G4double sc_r1_;
    G4double sc_r2_;
    G4double sc_z_;
    G4int sc_no_segments_;
    G4int sc_no_planes_;
    
    G4LogicalVolume *sc_log_;
    g4PSIScintillatorSD* SD_;
};

#endif
