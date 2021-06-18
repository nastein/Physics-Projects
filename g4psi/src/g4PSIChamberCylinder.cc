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
#include "g4PSIChamberCylinder.hh"
#include "TVectorD.h"


g4PSIChamberCylinder::g4PSIChamberCylinder(G4String label, G4Material *mat, G4double chamber_r, G4double chamber_t, G4double front_z, G4double win_height) : g4PSIChamberBase(label)
{
    init();
    
    chamber_material_ = mat;
    chamber_outer_radius_ = chamber_r;
    chamber_wall_thickness_ = chamber_t;
    entrance_window_zpos_ = front_z;
    exit_window_height_ = win_height;
}


g4PSIChamberCylinder::g4PSIChamberCylinder(G4String label) : g4PSIChamberBase(label)
{
    init();
}


void g4PSIChamberCylinder::init() {
    using namespace CLHEP;
    
    chamber_material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    
    chamber_sd_ = NULL;
    full_solid_chamber_log_ = NULL;
    full_vacuum_chamber_log_ = NULL;
    chamber_lid_log_ = NULL;
    chamber_exit_window_log_ = NULL;
    chamber_entrance_window_log_ = NULL;
    
    chamber_wall_thickness_ = 1 * cm;
    chamber_outer_radius_ = 12 * cm;
    chamber_lid_height_ = 1 * cm;
    chamber_height_ = 60 * cm;
    
    exit_window_max_angle_ = 110.*deg;
    exit_window_height_ = 15*cm;
    chamber_exit_window_thickness_ = 200*um;
    
    entrance_window_zpos_ = -15*cm;
    entrance_tube_inner_radius_ = 3*cm;
    entrance_tube_outer_radius_ = 4*cm;
    entrance_window_thickness_ = 200*um;
}


void g4PSIChamberCylinder::Info() {
    InfoTitle("Chamber (cylinder)");
}


void g4PSIChamberCylinder::Placement() {
    
    using namespace CLHEP;

    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    G4double chamber_lid_radius = chamber_outer_radius_ + 1 * cm;
    G4Tubs* chamber_lid = new G4Tubs("scat_lid", 0.0, chamber_lid_radius, chamber_lid_height_/2.0, 0, twopi);
    
    // Solid volume, including solid chamber walls and entrance flange
    G4RotationMatrix* RotUnion = new G4RotationMatrix;
    RotUnion->rotateY(90*deg);
    
    G4VSolid* full_solid_chamber = NULL;
    G4VSolid* full_vacuum_chamber = NULL;
    
    full_solid_chamber = new G4Tubs("scat_wall", 0.0, chamber_outer_radius_, chamber_height_/2.0, 0, twopi);
    full_vacuum_chamber = new G4Tubs("scat_vacuum", 0.0, chamber_outer_radius_ - chamber_wall_thickness_, chamber_height_/2.0, 0, twopi);
    
    
    if (entrance_window_zpos_ < -chamber_outer_radius_) {
        
        // mount entrance window and flange
        
        G4double entrance_tube_height = -entrance_window_zpos_ - sqrt(pow(chamber_outer_radius_-chamber_wall_thickness_,2) - pow(entrance_tube_inner_radius_, 2));
        G4Tubs* entrance_wall = new G4Tubs("scat_et_wall", 0.0, entrance_tube_outer_radius_, entrance_tube_height/2.0, 0, twopi);
        full_solid_chamber = new G4UnionSolid("chamber_wall", full_solid_chamber, entrance_wall, RotUnion, G4ThreeVector(entrance_window_zpos_ + entrance_tube_height/2,0,0));
        
        G4Tubs* entrance_inner_volume = new G4Tubs("scat_et_vacuum", 0.0, entrance_tube_inner_radius_, entrance_tube_height/2.0, 0, twopi);
        full_vacuum_chamber = new G4UnionSolid("chamber_vacuum", full_vacuum_chamber, entrance_inner_volume, RotUnion, G4ThreeVector(entrance_window_zpos_ + entrance_tube_height/2,0,0));
        
        G4Tubs* chamber_entrance_window = new G4Tubs("scat_entrance_window", 0, (entrance_tube_inner_radius_+entrance_tube_outer_radius_)/2, entrance_window_thickness_/2., 0, twopi);
        
        chamber_entrance_window_log_ = new G4LogicalVolume(chamber_entrance_window, Kapton, "chamber_entrance_window.log");
    }

    if (exit_window_height_ > 0) {
        
        // mount exit window
        
        G4Tubs* chamber_exitgap = new G4Tubs("scat_exit_vacuum", chamber_outer_radius_ - chamber_wall_thickness_, chamber_outer_radius_, exit_window_height_/2., -exit_window_max_angle_, 2.*exit_window_max_angle_);
        full_vacuum_chamber = new G4UnionSolid("chamber_plus_gap_vacuum", full_vacuum_chamber, chamber_exitgap);
        G4Tubs* chamber_exit_window = new G4Tubs("scat_exit_window", chamber_outer_radius_, chamber_outer_radius_ + chamber_exit_window_thickness_, exit_window_height_/2., -exit_window_max_angle_, 2.*exit_window_max_angle_);
        chamber_exit_window_log_ = new G4LogicalVolume(chamber_exit_window, Kapton, "chamber_exit_window.log");
    }
    
    
    // Chamber entrance and exit windows

    
    full_solid_chamber_log_ = new G4LogicalVolume(full_solid_chamber, chamber_material_, "chamber_wall.log");
    full_vacuum_chamber_log_ = new G4LogicalVolume(full_vacuum_chamber, Vacuum, "chamber_vacuum.log");
    chamber_lid_log_ = new G4LogicalVolume(chamber_lid, chamber_material_, "chamber_lid.log");


    
    G4RotationMatrix* rot = new G4RotationMatrix;
    rot->rotateX(90*deg);
    rot->rotateZ(90*deg);
    
    new G4PVPlacement(rot, G4ThreeVector(0,0,0), full_solid_chamber_log_, "full_al_chamber.vol", mother_volume_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0,0,0), full_vacuum_chamber_log_, "full_vacuum_chamber.vol", full_solid_chamber_log_, false, 1);
    if (chamber_exit_window_log_) {
        new G4PVPlacement(rot, G4ThreeVector(0,0,0), chamber_exit_window_log_, "scat_exit_window.vol", mother_volume_, false, 2);
    }
    if (chamber_entrance_window_log_) {
        new G4PVPlacement(0, G4ThreeVector(0,0, entrance_window_zpos_ - entrance_window_thickness_/2), chamber_entrance_window_log_, "scat_entrance_window.vol", mother_volume_, false, 3);
    }
    new G4PVPlacement(rot, G4ThreeVector(0, chamber_height_/2 + chamber_lid_height_/2,0), chamber_lid_log_, "scat_top_lid.vol", mother_volume_, false, 4);
    new G4PVPlacement(rot, G4ThreeVector(0, -chamber_height_/2 - chamber_lid_height_/2,0), chamber_lid_log_, "scat_bottom_lid.vol", mother_volume_, false, 5);
    
    G4Colour colour_chamber;
    G4Colour colour_chamber_window;
    G4Colour colour_chamber_air;
    G4Colour::GetColour("scattering_chamber", colour_chamber);
    G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
    G4Colour::GetColour("scattering_chamber_airgap", colour_chamber_air);
    G4VisAttributes* scatVisAtt = new G4VisAttributes(colour_chamber);
    G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
    G4VisAttributes* scatGVisAtt = new G4VisAttributes(colour_chamber_air);
    
    scatVisAtt->SetVisibility(true);
    scatVisAtt->SetForceLineSegmentsPerCircle(36);
    scatWVisAtt->SetVisibility(true);
    scatWVisAtt->SetForceLineSegmentsPerCircle(50);
    scatGVisAtt->SetVisibility(true);
    
    full_solid_chamber_log_->SetVisAttributes(scatVisAtt);
    chamber_lid_log_->SetVisAttributes(scatVisAtt);
    detector_parts->AttachBiasingOperator(full_solid_chamber_log_);
    detector_parts->AttachBiasingOperator(chamber_lid_log_);

    full_vacuum_chamber_log_->SetVisAttributes(scatGVisAtt);
    full_vacuum_chamber_log_->SetUserLimits(GetChamberStepLimit());
    
    if (chamber_entrance_window_log_) {
        chamber_entrance_window_log_->SetVisAttributes(scatWVisAtt);
        detector_parts->AttachBiasingOperator(chamber_entrance_window_log_);
    }
    if (chamber_exit_window_log_) {
        chamber_exit_window_log_->SetVisAttributes(scatWVisAtt);
        detector_parts->AttachBiasingOperator(chamber_exit_window_log_);
    }
    

    // determine where and how the target is to be placed:
    // --------------------------------------------------
    
    SetMotherLVforTarget(full_vacuum_chamber_log_);
    SetChamberBeamPosition(G4ThreeVector(0,0,0));
    SetChamberBeamRotation(new G4RotationMatrix(G4ThreeVector(0,0,1), G4ThreeVector(1,0,0), G4ThreeVector(0,1,0)));
}


void g4PSIChamberCylinder::SetSD(G4SDManager *SDman) {
    
    G4String ChamberSDname = "g4PSI/" + label_;
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(ChamberSDname);
    if (sd == NULL) {
        sd = chamber_sd_ = new g4PSITargetSD(ChamberSDname, label_ + "_Collection", GetChamberSDpMin(), GetChamberSDthetaMin());
        SDman->AddNewDetector( sd );
    };
    
    if (full_solid_chamber_log_) full_solid_chamber_log_->SetSensitiveDetector( sd );
    if (full_vacuum_chamber_log_) full_vacuum_chamber_log_->SetSensitiveDetector( sd );
    if (chamber_lid_log_) chamber_lid_log_->SetSensitiveDetector( sd );
    if (chamber_exit_window_log_) chamber_exit_window_log_->SetSensitiveDetector( sd );
    if (chamber_entrance_window_log_) chamber_entrance_window_log_->SetSensitiveDetector( sd );
}


void g4PSIChamberCylinder::Write() {
    TVectorD v1(1);
    v1[0] = 0;
    v1.Write(label_);
}


void g4PSIChamberCylinder::InitTree(TTree *T) {
    if (chamber_sd_) chamber_sd_->InitTree(T);
}


void g4PSIChamberCylinder::DeleteEventData() {
    if (chamber_sd_) chamber_sd_->DeleteEventData();
}

