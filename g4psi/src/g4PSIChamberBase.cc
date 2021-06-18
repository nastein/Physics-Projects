#include "g4PSIChamberBase.hh"

g4PSIChamberBase::g4PSIChamberBase(G4String label) : g4PSIDetectorBase(label)
{
    
    min_p_ = 10 * CLHEP::MeV;
    min_theta_ = 10 * CLHEP::deg;
    
    mother_for_target_lv_ = NULL;
    chamber_beam_position_ = G4ThreeVector(0,0,0);  // position of the target in the mother volume
    chamber_beam_rotation_ = NULL;
    
    chamber_step_limit_ = new G4UserLimits(0.5*CLHEP::mm);
}


G4LogicalVolume* g4PSIChamberBase::GetMotherLVforTarget() {
    return mother_for_target_lv_;
}


G4Material* g4PSIChamberBase::GetMaterialMotherLVforTarget() {
    if (mother_for_target_lv_) {
        return mother_for_target_lv_->GetMaterial();
    } else {
        return NULL;
    }
}


G4ThreeVector g4PSIChamberBase::GetChamberBeamPosition() {
    return chamber_beam_position_;
}


G4RotationMatrix* g4PSIChamberBase::GetChamberBeamRotation() {
    return chamber_beam_rotation_;
}


G4UserLimits* g4PSIChamberBase::GetChamberStepLimit() {
    return chamber_step_limit_;
}


void g4PSIChamberBase::SetMotherLVforTarget(G4LogicalVolume *lv) {
    mother_for_target_lv_ = lv;
}


void g4PSIChamberBase::SetChamberBeamPosition(G4ThreeVector pos) {
    chamber_beam_position_ = pos;
}


void g4PSIChamberBase::SetChamberBeamRotation(G4RotationMatrix *rot) {
    chamber_beam_rotation_ = rot;
}


void g4PSIChamberBase::SetChamberSDLimits(G4double p, G4double theta) {
    min_p_ = p;
    min_theta_ = theta;
}


G4double g4PSIChamberBase::GetChamberSDpMin() {return min_p_;}


G4double g4PSIChamberBase::GetChamberSDthetaMin() {return min_theta_;}
