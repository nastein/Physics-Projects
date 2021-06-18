#ifndef g4PSINaI_h
#define g4PSINaI_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSINaI
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

class g4PSINaI : public g4PSIDetectorBase {
    
public:
    g4PSINaI(G4String label, G4double z, G4double t_pb);
    ~g4PSINaI();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return nai_label_;};
    void Write();
    void Info();
    
    void InitTree(TTree *T);
    void DeleteEventData();
    void ReadEventData(G4HCofThisEvent *HCE);
    
private:
    G4String nai_label_;
    G4double nai_z_;
    G4int nai_n_;
    G4double nai_t_pb_;
    G4double bar_x_;
    G4double bar_y_;
    G4double bar_z_;
    
    G4LogicalVolume *nai_log_;
    g4PSIScintillatorSD* SD_;
};

#endif
