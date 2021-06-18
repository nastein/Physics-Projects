#include "g4PSITargetBase.hh"


g4PSITargetBase::g4PSITargetBase(G4String label) : g4PSIDetectorBase(label)
{
    init();
}


g4PSITargetBase::g4PSITargetBase(G4String label, g4PSIChamberBase* chamber) : g4PSIDetectorBase(label)
{
    init();
    target_chamber_ = chamber;
}


void g4PSITargetBase::init() {
    min_p_ = 10 * CLHEP::MeV;
    min_theta_ = 10 * CLHEP::deg;
    
    target_chamber_ = NULL;
    target_position_ = G4ThreeVector(0,0,0);
    target_rotation_ = NULL;
    
    target_step_limit_ = new G4UserLimits(0.5*CLHEP::mm);
}


void g4PSITargetBase::SetTargetPosition(G4ThreeVector pos) {
    target_position_ = pos;
}


void g4PSITargetBase::SetTargetRotation(G4RotationMatrix *rot) {
    target_rotation_ = rot;
}


void g4PSITargetBase::SetTargetChamber(g4PSIChamberBase *chamber) {
    target_chamber_ = chamber;
}


G4LogicalVolume* g4PSITargetBase::GetTargetMotherLV() {
    if (target_chamber_) {
        return target_chamber_->GetMotherLVforTarget();
    } else {
        return mother_volume_;
    }
}


G4Material* g4PSITargetBase::GetMaterialForTargetMotherLV() {
    if (target_chamber_) {
        return target_chamber_->GetMaterialMotherLVforTarget();
    } else {
        return NULL;
    }
}


G4bool g4PSITargetBase::is_in_chamber() {
    return target_chamber_ != NULL;
}


G4ThreeVector g4PSITargetBase::GetTargetPosition() {
    if (target_chamber_) {
        return target_position_ + target_chamber_->GetChamberBeamPosition();
    } else {
        return target_position_;
    }
}


G4RotationMatrix* g4PSITargetBase::GetTargetRotation() {
    
    // todo: needs to be tested
    if (target_chamber_ && target_chamber_->GetChamberBeamRotation()) {
        if (target_rotation_) {
            G4RotationMatrix* rot = new G4RotationMatrix(*target_rotation_);
            *rot *= *target_chamber_->GetChamberBeamRotation();
            return rot;
        } else {
            return target_chamber_->GetChamberBeamRotation();
        }
    } else {
        return target_rotation_;
    }
}