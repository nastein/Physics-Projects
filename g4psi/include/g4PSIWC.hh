#ifndef g4PSIWC_h
#define g4PSIWC_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIWC
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

class g4PSIWC : public g4PSIDetectorBase {

public:
    g4PSIWC(G4String label, G4double angle, G4double r, G4double x, G4double y);
    ~g4PSIWC();

public:
    void Placement();
    G4LogicalVolume* GetLog() {return wc_detector_log_;};
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return wc_label_;};
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    bool hitsWall(double vx, double vy, double vz, 
		  double x0, double y0, double z0);

private:
    G4String wc_label_;
    G4double wc_r_;
    G4double wc_z_;
    G4double wc_angle_;
    G4double wc_active_x_;
    G4double wc_active_y_;
    G4double wc_running_z_;

    G4LogicalVolume *wc_detector_log_;
    G4LogicalVolume *wc_assembly_log_;
    g4PSITrackerSD* SD_;

    G4LogicalVolume* wc_layer(G4Material *material, 
			      G4double dz, G4String name,
			      G4Colour col);

};

#endif
