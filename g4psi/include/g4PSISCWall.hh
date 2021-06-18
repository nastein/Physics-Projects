#ifndef g4PSISCWall_h
#define g4PSISCWall_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSISCWall
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

class g4PSISCWall : public g4PSIDetectorBase {
    
public:
    g4PSISCWall(G4String label,
                G4double width, G4double thickness, G4double length,
                G4double posX, G4double posY, G4RotationMatrix* rot = NULL,
                G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"));
    
    g4PSISCWall(G4String label,
                G4double width, G4double gap, G4double thickness,
                G4double distance, G4double target_z, G4double target_r,
                G4double safety_u1, G4double safety_d1, G4double safety_y1,
                G4double thetaXmin, G4double thetaXmax,
                G4double thetaYmin, G4double thetaYmax);
    
    g4PSISCWall(G4String label,
                G4double width, G4double gap, G4double thickness, G4double height,
                G4double distance, G4double theta, G4int Nup, G4int Ndown,
                G4double al_thickness = 0);
    
    ~g4PSISCWall();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4double GetThetaCenter() {return sc_theta_center_;};
    void Write();
    
    bool hitsWall(double vx, double vy, double vz);
    void Info();
    
    void InitTree(TTree *T);
    void DeleteEventData();
    G4bool HasHit() {return SD_->SDHasHit();};
    
private:
    G4double sc_width_;
    G4double sc_thickness_;
    G4double sc_gap_;
    G4double sc_distance_to_front_;
    G4double sc_dx_up_;
    G4double sc_dx_down_;
    G4double sc_thetaXmin_;
    G4double sc_thetaXmax_;
    G4double sc_thetaYmin_;
    G4double sc_thetaYmax_;
    G4double sc_al_thickness_;
    G4RotationMatrix *sc_rot_;
    G4Material* sc_material_;
    
    G4double sc_theta_center_;
    G4double sc_bar_length_;
    G4int sc_N_up_;
    G4int sc_N_down_;
    G4double sc_posX_;
    G4double sc_posZ_;
    G4int sc_mode_;
    
    G4LogicalVolume *sc_log_;
    g4PSIScintillatorSD* SD_;
    
};

#endif
