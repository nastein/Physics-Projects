#ifndef g4PSIChamberBase_h
#define g4PSIChamberBase_h 1

#include "globals.hh"
#include "G4UserLimits.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITargetSD.hh"


class g4PSIChamberBase : public g4PSIDetectorBase
{
public:
    
    g4PSIChamberBase(G4String label);
    G4LogicalVolume* GetMotherLVforTarget();
    G4ThreeVector GetChamberBeamPosition();
    G4RotationMatrix* GetChamberBeamRotation();
    G4Material* GetMaterialMotherLVforTarget();
    
protected:

    void SetMotherLVforTarget(G4LogicalVolume* lv);
    void SetChamberBeamPosition(G4ThreeVector pos);
    void SetChamberBeamRotation(G4RotationMatrix* rot);
    void SetChamberSDLimits(G4double p, G4double theta);

    G4UserLimits* GetChamberStepLimit();
    G4double GetChamberSDpMin();
    G4double GetChamberSDthetaMin();
    
private:
    
    G4LogicalVolume* mother_for_target_lv_;
    G4ThreeVector chamber_beam_position_;
    G4RotationMatrix *chamber_beam_rotation_;
    G4UserLimits* chamber_step_limit_;
    G4double min_p_;
    G4double min_theta_;

};

#endif
