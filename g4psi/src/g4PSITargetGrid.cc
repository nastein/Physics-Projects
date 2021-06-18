#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"

#include "g4PSIDetectorParts.hh"
#include "g4PSITargetGrid.hh"
#include "TVectorD.h"


using namespace CLHEP;


g4PSITargetGrid::g4PSITargetGrid(G4String label, g4PSIChamberBase* chamber) : g4PSITargetBase(label, chamber) {
    init();
}


void g4PSITargetGrid::init() {

    target_sd_ = NULL;
}


void g4PSITargetGrid::Info() {
    InfoTitle("Target (grid)");
}


void g4PSITargetGrid::Placement() {
    
    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    G4Material* plastic = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");
    G4Material* empty_mat = GetMaterialForTargetMotherLV();
    if (!empty_mat) {
        empty_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    }
    
    // --- geometry ---
    
    const G4double inches = 2.54 * cm;

    G4int ncell = 3;
    G4double grid_dx = 2.0 * cm;
    G4double grid_dy = 2.0 * cm;
    G4double grid_dz = (1./16.) * inches;
    G4double grid_w = 0.5 * mm;
    G4double grid_x = ncell * grid_dx + grid_w;
    G4double grid_y = ncell * grid_dy + grid_w;

    // target volume

    G4Box* target_assembly = new G4Box("grid_target", grid_x/2., grid_y/2., grid_dz/2.);
    target_log_ = new G4LogicalVolume(target_assembly, empty_mat, "grid_target_log");
    G4Box* target_v = new G4Box("grid_target_v", grid_w/2, grid_y/2., grid_dz/2.);
    target_v_log_ = new G4LogicalVolume(target_v, plastic, "grid_target_v_log");
    G4Box* target_h = new G4Box("grid_target_h", (grid_dx-grid_w)/2, grid_w/2., grid_dz/2.);
    target_h_log_ = new G4LogicalVolume(target_h, plastic, "grid_target_v_log");


    // placements
    // ----------
    
    G4LogicalVolume* mother_lv = GetTargetMotherLV();
    G4ThreeVector chamber_shift = GetTargetPosition();
    G4RotationMatrix* chamber_rot = GetTargetRotation();

    new G4PVPlacement(chamber_rot, G4ThreeVector(0, 0, 0) + chamber_shift, target_log_, "grid_target", mother_lv, false, 0);
    
    for (int ix = 0; ix <= ncell; ix++) {
        new G4PVPlacement(NULL, G4ThreeVector(-grid_x/2.+grid_w/2.+ix*grid_dx,0,0), target_v_log_, "grid_target_v", target_log_, false, 0);
        if (ix > 0) {
            for (int iy = 0; iy <= ncell; iy++) {
                double x = -grid_x/2.+grid_w/2.+(ix-0.5)*grid_dx;
                double y = -grid_y/2.+grid_w/2.+iy*grid_dy;
                new G4PVPlacement(NULL, G4ThreeVector(x, y, 0), target_h_log_, "grid_target_h", target_log_, false, 0);
            }
        }
    }
    
    //
    G4Colour target_color;
    G4Colour::GetColour("target", target_color);
    target_v_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_h_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_log_->SetVisAttributes(G4VisAttributes::Invisible);
    
    target_v_log_->SetUserLimits(target_step_limit_);
    target_h_log_->SetUserLimits(target_step_limit_);
    
    detector_parts->AttachBiasingOperator(target_v_log_);
    detector_parts->AttachBiasingOperator(target_h_log_);
}


void g4PSITargetGrid::SetSD(G4SDManager *SDman) {
    
    // if the target is placed into a target chamber, the sensitive detector
    // for the target is the one of the chamber if it is found
    
    G4String TargetSDname = "g4PSI/" + (is_in_chamber() ? target_chamber_->GetDetectorName() : label_);
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(TargetSDname);
    if (sd == NULL) {
        sd = target_sd_ = new g4PSITargetSD(TargetSDname, label_ + "_Collection", min_p_, min_theta_);
        SDman->AddNewDetector( sd );
    };
    
    if (target_v_log_) target_v_log_->SetSensitiveDetector( sd );
    if (target_h_log_) target_h_log_->SetSensitiveDetector( sd );
}


void g4PSITargetGrid::Write() {
    TVectorD v1(6);
    v1[0] = g4PSIDetectorParts::getInstance()->GetTargetState();
    v1.Write(label_);
}


void g4PSITargetGrid::InitTree(TTree *T) {
    if (target_sd_) target_sd_->InitTree(T);
}


void g4PSITargetGrid::DeleteEventData() {
    if (target_sd_) target_sd_->DeleteEventData();
}

