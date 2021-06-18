#ifndef g4PSISPS_h
#define g4PSISPS_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSISPS
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

class g4PSISPS : public g4PSIDetectorBase {
    
    
public:
    enum SPSWallType {front, rear};
    
    g4PSISPS(G4String label,
             SPSWallType type,
             G4double width,
             G4double gap_small,
             G4double gap_wide,
             G4double thickness,
             G4double height,
             G4double distance,
             G4double theta,
             G4double shift,
             G4int Nunits);
    
    ~g4PSISPS();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4double GetThetaCenter() {return sc_theta_center_;};
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    G4bool HasHit() {return SD_->SDHasHit();};
    
private:
    SPSWallType wall_type_;
    G4double sc_width_;
    G4double sc_thickness_;
    G4double sc_gap_small_;
    G4double sc_gap_wide_;
    G4double sc_distance_to_front_;
    G4double sps_r_rear_of_backing_;
    G4double upstream_shift_;
    G4double sps_w1_;
    G4double sps_w2_;
    G4double sps_wall_dx_;
    G4double sps_wall_dy_;
    G4double sps_wall_dz_;
    G4double backing_length_;
    G4double backing_al_length_;
    G4double backing_al_thickness_;
    G4double backing_al_width_;
    G4double backing_rc_length_;
    G4double backing_rc_thickness_;
    G4double backing_rc_width_;

    G4RotationMatrix *sps_rot_;
    G4double sign_;
    
    G4double sc_theta_center_;
    G4double sc_bar_length_;
    G4int sc_Nunits_;
    
    G4LogicalVolume *sc_log_;
    g4PSIScintillatorSD* SD_;
    
    enum orientation {X, Y, Z};
    void place_extrusion(G4double L, orientation o, G4double x, G4double y, G4double z, G4double dy = 40 * CLHEP::mm);
};

#endif
