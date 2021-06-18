#ifndef g4PSIBeamMonitor_h
#define g4PSIBeamMonitor_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIBeamMonitor
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

class g4PSIBeamMonitor : public g4PSIDetectorBase {

public:
    g4PSIBeamMonitor(G4String label, G4double z);
    ~g4PSIBeamMonitor();

public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return label_;};
    void Write();
    void Info();

    void InitTree(TTree *T);
    void DeleteEventData();
    void ReadEventData(G4HCofThisEvent *HCE);

private:
    const G4double sipm_height_ = 4*CLHEP::mm;
    const G4double sipm_width_ = 3*CLHEP::mm;
    
    G4String label_;
    G4double posz_;
    G4double gap_plate_;
    G4double gap_bar_;
    G4double sizx_;
    G4double sizy_;
    G4double sizz_;
    G4double sizw_;
    G4int N_;
    G4int nSC_vertical_bars_;
    G4int nSC_horizontal_bars_;
    G4double len_v_;
    G4double len_h_;
    
    G4double shiftX_;
    G4double shiftZ_;
    
    G4Material* m_;
    g4PSIScintillatorSD* SD_;
    G4LogicalVolume* bm_log1_;
    G4LogicalVolume* bm_log2_;
    G4LogicalVolume* bm_log3_;
    G4LogicalVolume* bm_sipm_log_;

};

#endif
