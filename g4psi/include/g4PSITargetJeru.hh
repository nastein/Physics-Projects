#ifndef g4PSITargetJeru_h
#define g4PSITargetJeru_h 1

#include "globals.hh"
#include "g4PSITargetBase.hh"
#include "g4PSITargetSD.hh"


class g4PSITargetJeru : public g4PSITargetBase
{
public:
    
    g4PSITargetJeru(G4String label, G4double tgt_R, G4double tgt_H, G4double base_r1, G4double base_r2, G4double base_r3, G4double base_h1, G4double base_h2, G4double base_h3, G4Material* base_material, G4double tgt_wall_t, G4Material* tgt_wall_material, g4PSIChamberBase* chamber = NULL);
    
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
    G4LogicalVolume* target_nipple_log_;
    G4LogicalVolume* target_pipe_h_log_;
    G4LogicalVolume* target_pipe_v_log_;
    G4LogicalVolume* target_probe_p_log_;
    G4LogicalVolume* target_probe_t_log_;
    
    g4PSITargetSD* target_sd_;
    
    G4double tgt_R_;
    G4double tgt_H_;
    
    G4double base_r1_;
    G4double base_r2_;
    G4double base_r3_;
    
    G4double base_h1_;
    G4double base_h2_;
    G4double base_h3_;
    G4double tgt_wall_t_;  // Dany Horovitz e-mail, 03/13/16
    
    G4Material* base_material_;
    G4Material* tgt_wall_material_;

};

#endif
