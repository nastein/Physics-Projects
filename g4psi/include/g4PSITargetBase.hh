#ifndef g4PSITargetBase_h
#define g4PSITargetBase_h 1

#include "globals.hh"
#include "G4UserLimits.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITargetSD.hh"
#include "g4PSIChamberBase.hh"


class g4PSITargetBase : public g4PSIDetectorBase {
    
public:
    
    g4PSITargetBase(G4String label);
    g4PSITargetBase(G4String label, g4PSIChamberBase* chamber);
    
    void SetTargetChamber(g4PSIChamberBase* chamber);

    // target position relative to the logical volume in
    // which the target will be placed
    void SetTargetPosition(G4ThreeVector pos);
    void SetTargetRotation(G4RotationMatrix* rot);
    
    G4ThreeVector GetTargetPosition();
    G4RotationMatrix* GetTargetRotation();
    G4LogicalVolume* GetTargetMotherLV();
    G4Material* GetMaterialForTargetMotherLV();
    
protected:
    
    G4bool is_in_chamber();
    g4PSIChamberBase* target_chamber_;
    G4UserLimits *target_step_limit_;
    G4double min_p_;
    G4double min_theta_;
    
private:
    
    void init();
    G4ThreeVector target_position_;
    G4RotationMatrix *target_rotation_;
};

#endif
