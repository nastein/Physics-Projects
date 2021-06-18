#ifndef g4PSIGEM_h
#define g4PSIGEM_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIGEM
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITrackerSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"

class g4PSIGEM : public g4PSIDetectorBase {
    
public:
    g4PSIGEM(G4String, G4double, G4double, G4double);
    g4PSIGEM(G4String, G4double, G4double, G4double, G4double, G4double, G4double);
    ~g4PSIGEM();
    
public:
    void Placement();
    G4LogicalVolume* GetLog() {return gem_detector_log_;};
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return gem_label_;};
    //  void SetColour(G4Colour c) {gem_colour_ = c;};
    void Write();
    void Info();
    double GetZPos();
    
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
    G4String gem_label_;
    G4double gem_dx_;
    G4double gem_dy_;
    G4double gem_dx_frame_;
    G4double gem_dy_frame_;
    G4double gem_r_;
    G4double gem_running_z_;
    G4double gem_angle_;
    
    G4LogicalVolume *gem_detector_log_;
    G4LogicalVolume *gem_assembly_log_;
    g4PSITrackerSD* SD_;
    
    G4LogicalVolume *gem_layer(G4Material *material,
                               G4double thickness_weight,
                               G4double dz,
                               G4String name, 
                               G4LogicalVolume *mother_volume,
                               G4Colour col);
};

#endif
