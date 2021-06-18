#ifndef g4PSIChamberCylinder_h
#define g4PSIChamberCylinder_h 1

#include "globals.hh"
#include "g4PSIChamberBase.hh"
#include "g4PSITargetSD.hh"


class g4PSIChamberCylinder : public g4PSIChamberBase
{
public:
    
    g4PSIChamberCylinder(G4String label);
    g4PSIChamberCylinder(G4String label, G4Material *mat, G4double chamber_r, G4double chamber_t, G4double front_z = 0, G4double win_height = 0);
    
    void Placement();
    double GetZPos() {return 0;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:

    void init();
    
    g4PSITargetSD* chamber_sd_;

    G4Material* chamber_material_;

    G4double chamber_outer_radius_;
    G4double chamber_wall_thickness_;
    G4double chamber_lid_height_;
    G4double chamber_height_;
    
    G4double exit_window_max_angle_;
    G4double exit_window_height_;
    G4double chamber_exit_window_thickness_;
    
    G4double entrance_window_zpos_;
    G4double entrance_tube_inner_radius_;
    G4double entrance_tube_outer_radius_;
    G4double entrance_window_thickness_;
    
    G4LogicalVolume* full_solid_chamber_log_;
    G4LogicalVolume* full_vacuum_chamber_log_;
    G4LogicalVolume* chamber_lid_log_;
    G4LogicalVolume* chamber_exit_window_log_;
    G4LogicalVolume* chamber_entrance_window_log_;
};

#endif
