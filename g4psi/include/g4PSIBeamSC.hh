#ifndef g4PSIBeamSC_h
#define g4PSIBeamSC_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIBeamSC
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

class g4PSIBeamSC : public g4PSIDetectorBase {
    
public:
    g4PSIBeamSC(G4String, G4double, G4double, G4double, G4double, G4double, G4double);
    g4PSIBeamSC(G4String name,
                G4double angle,
                G4double fiber_sz,  // diameter, including fiber cladding; < 0 for square fiber
                G4double fiber_cladding,  // thickness in radial direction
                G4double fiber_spacing,
                G4double wx,
                G4double wy,
                G4double offsetx,
                G4double offsety,
                G4double posZ);
    
    ~g4PSIBeamSC();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return label_;};
    void Write();
    void Info();
    
    G4bool HasHit() {return SD_->SDHasHit();};
    
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
    
    G4String label_;
    G4double angle_;
    G4double fiber_sz_;
    G4double fiber_cladding_;
    G4double fiber_spacing_;
    G4double wx_;
    G4double wy_;
    G4double offsetx_;
    G4double offsety_;
    G4double posz_;
    G4int NFibers_;
    G4bool circular_fiber_;
    
    // G4String sc_label_;
    // G4double sc_angle_;
    // G4double sc_fiber_size_;
    // G4double sc_z_pos_;
    // G4double sc_width_;
    // G4double sc_height_;
    // G4double sc_detector_x_;
    // G4double sc_detector_y_;
    // G4double sc_offset_;
    // G4int sc_N_;
    
    G4LogicalVolume *sc_log_;
    g4PSIScintillatorSD* SD_;
    
};

#endif
