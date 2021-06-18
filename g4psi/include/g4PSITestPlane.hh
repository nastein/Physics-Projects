#ifndef g4PSITestPlane_h
#define g4PSITestPlane_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSITestPlane
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

class g4PSITestPlane : public g4PSIDetectorBase {

public:
    g4PSITestPlane(G4String label, G4double angle, G4double r, G4double x, G4double y, G4LogicalVolume *m = NULL);
    ~g4PSITestPlane();

public:
    void Placement();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return testplane_label_;};
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();

private:
    G4String testplane_label_;
    G4double testplane_r_;
    G4double testplane_angle_;
    G4double testplane_active_x_;
    G4double testplane_active_y_;

    g4PSITrackerSD* SD_;

};

#endif
