#ifndef g4PSITargetCylinder_h
#define g4PSITargetCylinder_h 1

#include "globals.hh"
#include "g4PSITargetBase.hh"
#include "g4PSITargetSD.hh"


class g4PSITargetCylinder : public g4PSITargetBase
{
public:
    
    g4PSITargetCylinder(G4String label);
    g4PSITargetCylinder(G4String label, g4PSIChamberBase* chamber);
    
    void Placement();
    double GetZPos() {return 0;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
    
    void init();
    
    G4LogicalVolume* target_log_;
    G4LogicalVolume* target_wall_log_;
    G4LogicalVolume* target_vacuum_log_;
    G4LogicalVolume* target_si_log_;
    
    g4PSITargetSD* target_sd_;

    G4double target_dx_;
    G4double target_dy_;
    G4double target_dz_;

    G4double target_entrance_cap_thickness_;
    G4double target_exit_cap_thickness_;
    G4double target_wall_thickness_;
    
    G4double target_mylar_thickness_;
    G4double target_flask_mylar_gap_;
    
};

#endif
