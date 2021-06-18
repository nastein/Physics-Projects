#ifndef g4PSITargetUM_h
#define g4PSITargetUM_h 1

#include "globals.hh"
#include "g4PSITargetBase.hh"
#include "g4PSITargetSD.hh"


class g4PSITargetUM : public g4PSITargetBase
{
public:
    
    g4PSITargetUM(G4String label, g4PSIChamberBase* chamber = NULL);
    
    void Placement();
    double GetZPos() {return 0;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
        
    G4LogicalVolume* target_log_;
    G4LogicalVolume* target_base_log_;
    G4LogicalVolume* target_cell_film_log_;
    G4LogicalVolume* target_si_log_;
    G4LogicalVolume* radial_stub_cu_log_;
    G4LogicalVolume* radial_stub_ss_log_;
    G4LogicalVolume* support_tube_cu_log_;
    G4LogicalVolume* support_tube_ss_log_;
    G4LogicalVolume* cu_tube_log_;
    
    g4PSITargetSD* target_sd_;
    
    G4double tgt_R_;
    G4double tgt_H_;
    
    G4double base_r1_;
    G4double base_r2_;
    G4double base_r3_;
    
    G4double base_h1_;
    G4double base_h3_;
    G4double tgt_wall_t_;  // Dany Horovitz e-mail, 03/13/16
    
    G4double tgt_cell_gap_;
    G4double tgt_cell_height_;
    
    G4double cu_tube_od_;
    G4double cu_tube_t_;
    G4double cu_tube_l_;
    G4double support_tube_cu_od_;
    G4double support_tube_cu_t_;
    G4double support_tube_cu_l_;
    G4double support_tube_ss_od_;
    G4double support_tube_ss_t_;
    G4double support_tube_ss_l_;
    G4double radial_stub_od_;
    G4double radial_stub_t_;
    G4double radial_stub_l_;
    G4double insulation_t_;
    
    G4Material* base_material_;
    G4Material* tgt_wall_material_;

};

#endif
