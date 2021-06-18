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
#include "g4PSITargetCylinder.hh"
#include "TVectorD.h"


using namespace CLHEP;


g4PSITargetCylinder::g4PSITargetCylinder(G4String label) : g4PSITargetBase(label) {
    init();
}


g4PSITargetCylinder::g4PSITargetCylinder(G4String label, g4PSIChamberBase* chamber) : g4PSITargetBase(label, chamber) {
    init();
}


void g4PSITargetCylinder::init() {

    target_sd_ = NULL;

    target_log_ = NULL;
    target_wall_log_ = NULL;
    target_vacuum_log_ = NULL;
    target_si_log_ = NULL;
    
    target_dx_ = 3.*cm;
    target_dy_ = 3.*cm;
    target_dz_ = 4.*cm;
    
    target_entrance_cap_thickness_ = 0.0625 * mm;
    target_exit_cap_thickness_ = 0.0625 * mm;
    target_wall_thickness_ = 0.0625 * mm;
    
    target_mylar_thickness_ = 0.050 * mm;
    target_flask_mylar_gap_ = 1 * cm;
}


void g4PSITargetCylinder::Info() {
    InfoTitle("Target (cylinder)");
}


void g4PSITargetCylinder::Placement() {
    
    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    // --- material
    
    G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
    
    G4Material* empty_mat = GetMaterialForTargetMotherLV();
    if (!empty_mat) {
        empty_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    }
    
    G4Material* target_material = NULL;
    if (detector_parts->GetTargetState() == g4PSIDetectorParts::kEmpty) {
        std::cout << "empty\n";
        target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
    } else if (detector_parts->GetTargetState() == g4PSIDetectorParts::kFull) {
        std::cout << "full\n";
        target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
    } else {
        G4Exception("g4PSITargetCylinder",
                    "Error: unknown target state.",
                    FatalException, "");
    };
    std::cout << "Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;
    
    // --- geometry
    
    G4EllipticalTube* target_si_geom = new G4EllipticalTube("target_si.geom",target_dx_ + target_wall_thickness_ + target_flask_mylar_gap_ + target_mylar_thickness_,target_dy_ + target_wall_thickness_ + target_flask_mylar_gap_ + target_mylar_thickness_,(target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/2 + target_flask_mylar_gap_ + target_mylar_thickness_);
    G4EllipticalTube* target_vacuum_geom = new G4EllipticalTube("target_vacuum.geom",target_dx_ + target_wall_thickness_ + target_flask_mylar_gap_,target_dy_ + target_wall_thickness_ + target_flask_mylar_gap_,(target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/2 +target_flask_mylar_gap_);
    G4EllipticalTube* target_wall_geom = new G4EllipticalTube("target_wall_geom",target_dx_ + target_wall_thickness_,target_dy_ + target_wall_thickness_,(target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/ 2.0);
    G4EllipticalTube* target_cell_geom = new G4EllipticalTube("target_geom",target_dx_,target_dy_,target_dz_ / 2.0);
    
    target_log_ = new G4LogicalVolume(target_cell_geom, target_material, "target_log", 0,0,0);
    target_wall_log_ = new G4LogicalVolume(target_wall_geom, Kapton, "target_wall_log", 0,0,0);
    target_vacuum_log_ = new G4LogicalVolume(target_vacuum_geom, empty_mat, "target_vacuum_log", 0,0,0);
    target_si_log_ = new G4LogicalVolume(target_si_geom, Mylar, "target_si_log", 0,0,0);
    
    G4Colour target_color;
    G4Colour::GetColour("target", target_color);
    G4Colour target_si_color;
    G4Colour::GetColour("target_si", target_si_color);
    G4Colour test_color;
    G4Colour::GetColour("white", test_color);
    target_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_wall_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_vacuum_log_->SetVisAttributes(G4VisAttributes::Invisible);
    target_si_log_->SetVisAttributes(new G4VisAttributes(target_si_color));
    
    new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                      target_log_, "target", target_wall_log_, false, 10);
    new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                      target_wall_log_, "target_wall", target_vacuum_log_, false, 11);
    new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                      target_vacuum_log_, "target_vacuum", target_si_log_, false, 12);
    new G4PVPlacement(GetTargetRotation(), GetTargetPosition(),
                      target_si_log_, "target_si", GetTargetMotherLV(), false, 13);
    
    detector_parts->AttachBiasingOperator(target_log_);
    detector_parts->AttachBiasingOperator(target_wall_log_);
    detector_parts->AttachBiasingOperator(target_vacuum_log_);
    detector_parts->AttachBiasingOperator(target_si_log_);
    
    target_log_->SetUserLimits(target_step_limit_);
    
    detector_parts->AttachBiasingOperator(target_log_);
    detector_parts->AttachBiasingOperator(target_wall_log_);
    detector_parts->AttachBiasingOperator(target_si_log_);
}


void g4PSITargetCylinder::SetSD(G4SDManager *SDman) {
    
    // if the target is placed into a target chamber, the sensitive detector
    // for the target is the one of the chamber if it is found
    
    G4String TargetSDname = "g4PSI/" + (is_in_chamber() ? target_chamber_->GetDetectorName() : label_);
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(TargetSDname);
    if (sd == NULL) {
        sd = target_sd_ = new g4PSITargetSD(TargetSDname, label_ + "_Collection", min_p_, min_theta_);
        SDman->AddNewDetector( sd );
    };
    
    if (target_log_) target_log_->SetSensitiveDetector( sd );
    if (target_wall_log_) target_wall_log_->SetSensitiveDetector( sd );
    if (target_vacuum_log_) target_vacuum_log_->SetSensitiveDetector( sd );
    if (target_si_log_) target_si_log_->SetSensitiveDetector( sd );
}


void g4PSITargetCylinder::Write() {
    TVectorD v1(6);
    v1[0] = g4PSIDetectorParts::getInstance()->GetTargetState();
    v1.Write(label_);
}


void g4PSITargetCylinder::InitTree(TTree *T) {
    if (target_sd_) target_sd_->InitTree(T);
}


void g4PSITargetCylinder::DeleteEventData() {
    if (target_sd_) target_sd_->DeleteEventData();
}

