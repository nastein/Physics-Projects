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
#include "g4PSITargetJeru.hh"
#include "TVectorD.h"


using namespace CLHEP;


g4PSITargetJeru::g4PSITargetJeru(G4String label, G4double tgt_R, G4double tgt_H, G4double base_r1, G4double base_r2, G4double base_r3, G4double base_h1, G4double base_h2, G4double base_h3, G4Material* base_material, G4double tgt_wall_t, G4Material* tgt_wall_material, g4PSIChamberBase* chamber) : g4PSITargetBase(label, chamber) {
    
    target_sd_ = NULL;
    
    tgt_R_ = tgt_R;
    tgt_H_ = tgt_H;
    base_r1_ = base_r1;
    base_r2_ = base_r2;
    base_r3_ = base_r3;
    base_h1_ = base_h1;
    base_h2_ = base_h2;
    base_h3_ = base_h3;
    base_material_ = base_material;
    
    tgt_wall_t_ = tgt_wall_t;
    tgt_wall_material_ = tgt_wall_material;
}


void g4PSITargetJeru::Info() {
    InfoTitle("Target (cylinder)");
    InfoParDouble("Target-cylinder height", tgt_H_ / mm, " mm");
    InfoParDouble("Target-cylinder radius", tgt_R_/ mm, " mm");
    InfoParDouble("Target wall thickness", tgt_wall_t_ / um, " um");
    InfoParString("Target wall material", tgt_wall_material_->GetName());
    InfoParString("Target cap material", base_material_->GetName());
}


void g4PSITargetJeru::Placement() {
    
    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    
    G4Material* empty_mat = GetMaterialForTargetMotherLV();
    if (!empty_mat) {
        empty_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    }
    
    G4Material* target_material = NULL;
    if (detector_parts->build[g4PSIDetectorParts::kTarget_TypeJeru]) {
        G4cout << "Full Target: Filling target with liquid hydrogen\n" << G4endl;
        target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
    } else {
        G4cout << "Empty Target: Filling target with hydryogen gas\n" << G4endl;
        target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
    };
    std::cout << "Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;
    
    // target volume
    
    G4double target_full_height = tgt_H_/2. + base_h2_ + base_h3_;
    G4double target_inner[] = {0, 0, 0, 0, 0, 0};
    G4double target_outer[] = {base_r2_, base_r2_, tgt_R_, tgt_R_, base_r2_, base_r2_};
    G4double target_z[] = {-target_full_height, -tgt_H_/2., -tgt_H_/2., tgt_H_/2., tgt_H_/2., target_full_height};
    G4Polycone* target_volume = new G4Polycone("target_cell_volume", 0, 2*pi, 6, target_z, target_inner, target_outer);
    
    target_log_ = new G4LogicalVolume(target_volume, target_material, "target_cell_volume.log", 0, 0, 0);
    
    // base
    
    G4Tubs* lip = new G4Tubs("target_cell_base_1", base_r2_, base_r3_, base_h3_/2, 0., twopi);
    G4Tubs* cap = new G4Tubs("target_cell_base_3",        0, base_r1_, base_h1_/2, 0., twopi);
    
    G4VSolid* base = lip;
    if (base_h2_ > 0) {
        G4Tubs* ring = new G4Tubs("target_cell_base_2", base_r2_, base_r1_, base_h2_/2, 0., twopi);
        base = new G4UnionSolid("target_cell_base12", base, ring, NULL, G4ThreeVector(0., 0., (base_h2_+base_h3_)/2));
    }
    base = new G4UnionSolid("target_cell_base_top", base, cap, NULL, G4ThreeVector(0., 0., (base_h1_+base_h3_)*mm/2 + base_h2_));

    
    // for now do not include holes in the base, assume solid base
    
    //        G4Tubs* base_hole = new G4Tubs("target_base_hole", 0, 6.8*mm/2, 8.6*mm/2, 0., twopi);
    //        for (G4int i = 0; i < 2; i++) {
    //            G4double r = base_inner_radius + 8.6*mm/2;
    //            G4RotationMatrix* Rot1 = new G4RotationMatrix;
    //            Rot1->rotateY(90*deg);
    //            Rot1->rotateX(90*deg + angle1);
    //            G4RotationMatrix* Rot2 = new G4RotationMatrix;
    //            Rot2->rotateY(90*deg);
    //            Rot2->rotateX(90*deg + angle2);
    //            target_base_top = new G4SubtractionSolid("target_cell_base.log", target_base_top, base_hole, Rot1, G4ThreeVector(r*sin(angle1),r*cos(angle1), 8.5*mm));
    //            target_base_bottom = new G4SubtractionSolid("target_cell_base.log", target_base_bottom, base_hole, Rot2, G4ThreeVector(r*sin(angle2),r*cos(angle2),8.5*mm));
    //            angle1 *= -1;
    //            angle2 *= -1;
    //        }
    
    target_base_log_ = new G4LogicalVolume(base, base_material_, "target_cell_base.log", 0, 0, 0);
    
    // film
    
    G4Tubs* target_cell_film = new G4Tubs("target_cell_film", tgt_R_, tgt_R_ + tgt_wall_t_, tgt_H_/2, 0, 2*pi);
    target_cell_film_log_ = new G4LogicalVolume(target_cell_film, tgt_wall_material_, "target_cell_film.log", 0, 0, 0);
    
    // nipple
    
    G4Tubs* nipple = new G4Tubs("target_nipple", 0., 25*mm/2, 20*mm/2, 0, 2*pi);
    target_nipple_log_ = new G4LogicalVolume(nipple, Kapton, "target_nipple.log", 0, 0, 0);
    
    // pipes
    
    G4double pipe_length_h = 70.5*mm - 61*mm/2 - 25*mm/2;
    G4Tubs* pipe_h = new G4Tubs("target_pipe_h", 0., 10*mm/2, pipe_length_h/2, 0, 2*pi);
    target_pipe_h_log_ = new G4LogicalVolume(pipe_h, Kapton, "target_pipe.log", 0, 0, 0);
    G4double pipe_length_v1 = 104*mm;
    G4Tubs* pipe_v1 = new G4Tubs("target_pipe_v1", 0., 10*mm/2, pipe_length_v1/2, 0, 2*pi);
    target_pipe_v_log_ = new G4LogicalVolume(pipe_v1, Kapton, "target_pipe.log", 0, 0, 0);
    
    // pressure probe
    
    const G4double probe_p_length1 = 8.5 * mm;
    const G4double probe_p_length2 = 20 * mm;
    const G4double probe_p_r1 = 10 * mm / 2;
    const G4double probe_p_r2 = 15 * mm / 2;
    
    const G4double probe_p_inner[] = {0, 0, 0, 0};
    const G4double probe_p_outer[] = {probe_p_r1, probe_p_r1, probe_p_r2, probe_p_r2};
    const G4double probe_p_z[] = {0, probe_p_length1, probe_p_length1, probe_p_length2};
    G4Polycone* probe_p_volume = new G4Polycone("probe_p_cell_volume", 0, 2*pi, 4, probe_p_z, probe_p_inner, probe_p_outer);
    
    target_probe_p_log_ = new G4LogicalVolume(probe_p_volume, Kapton, "probe_p_cell_volume.log", 0, 0, 0);
    
    // temperature probe
    
    const G4double probe_t_length1 = 8.5 * mm;
    const G4double probe_t_length2 = 30 * mm;
    const G4double probe_t_r1 = 10 * mm / 2;
    const G4double probe_t_r2 = 15 * mm / 2;
    
    const G4double probe_t_inner[] = {0, 0, 0, 0};
    const G4double probe_t_outer[] = {probe_t_r1, probe_t_r1, probe_t_r2, probe_t_r2};
    const G4double probe_t_z[] = {0, probe_t_length1, probe_t_length1, probe_t_length2};
    G4Polycone* probe_t_volume = new G4Polycone("probe_t_cell_volume", 0, 2*pi, 4, probe_t_z, probe_t_inner, probe_t_outer);
    
    target_probe_t_log_ = new G4LogicalVolume(probe_t_volume, Kapton, "probe_t_cell_volume.log", 0, 0, 0);
    
    // logical volume for entire target assembly
    
    G4double assembly_r = 90 * mm;
    G4double assembly_h = 180 * mm;
    G4double assembly_shift = -115 * mm;
    G4Tubs* assembly = new G4Tubs("target_assembly_tub", 0., assembly_r, assembly_h, 0, 2*pi);
    G4LogicalVolume* assembly_log = new G4LogicalVolume(assembly, empty_mat, "target_assembly_log");
    
    // placements
    // ----------
    
    G4LogicalVolume* mother_lv = GetTargetMotherLV();
    G4ThreeVector chamber_shift = GetTargetPosition();
    G4RotationMatrix* chamber_rot = NULL;
    if (GetTargetRotation()) {
        chamber_rot = new G4RotationMatrix(*GetTargetRotation());
    } else {
        chamber_rot = new G4RotationMatrix(); // unit matrix
    }
    G4RotationMatrix* rot_assembly = new G4RotationMatrix(G4ThreeVector(0,-1,0),G4ThreeVector(0,0,1),G4ThreeVector(-1,0,0));
    
    *rot_assembly *= *chamber_rot;
    
    G4ThreeVector tmp = (*rot_assembly) * G4ThreeVector(0, -assembly_shift, 0);
    
    new G4PVPlacement(rot_assembly, (*chamber_rot).invert() * G4ThreeVector(0, -assembly_shift, 0) + chamber_shift, assembly_log, "target_assembly", mother_lv, false, 0);
    
    const G4int NCELL = 2;
    G4double angle1 = 55*deg + 90*deg;
    G4double angle2 = 35*deg - 90*deg;
    G4RotationMatrix* Rot = new G4RotationMatrix;
    Rot->rotateX(180*deg);
    Rot->rotateZ(180*deg);
    G4RotationMatrix* Rot1A = new G4RotationMatrix;
    Rot1A->rotateY(90*deg);
    Rot1A->rotateX(90*deg + angle1);
    G4RotationMatrix* Rot1B = new G4RotationMatrix;
    Rot1B->rotateY(90*deg);
    Rot1B->rotateX(-90*deg - angle1);
    G4RotationMatrix* Rot2A = new G4RotationMatrix;
    Rot2A->rotateY(90*deg);
    Rot2A->rotateX(90*deg - angle2);
    G4RotationMatrix* Rot2B = new G4RotationMatrix;
    Rot2B->rotateY(90*deg);
    Rot2B->rotateX(-90*deg + angle2);
    for (G4int icell = 0; icell < NCELL; icell++) {
        G4double dz = icell * 124 * mm + assembly_shift;
        G4int id = (icell + 1) * 10; // target cell copy IDs are 10 and 20
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_log_,
                          Form("target_cell%d_volume", icell), assembly_log, false, id);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_cell_film_log_,
                          Form("target_cell%d_film", icell), assembly_log, false, id+1);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz+tgt_H_/2+base_h3_/2), target_base_log_,
                          Form("target_cell%d_base", icell), assembly_log, false, id+2);
        new G4PVPlacement(Rot, G4ThreeVector(0,0, dz-tgt_H_/2-base_h3_/2), target_base_log_,
                          Form("target_cell%d_base", icell), assembly_log, false, id+2);
        
        // place nipples
        
        G4double r = 70.5 * mm;
        
        new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm), target_nipple_log_, Form("target_nipple%d", 2*icell), assembly_log, false, id+3);
        new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm), target_nipple_log_, Form("target_nipple%d", 2*icell+1), assembly_log, false, id+3);
        
        // place pipes
        
        new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), assembly_log, false, id+4);
        new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), assembly_log, false, id+4);
        
        r = 61*mm/2 + pipe_length_h / 2;
        new G4PVPlacement(Rot1A, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell), assembly_log, false, id+5);
        new G4PVPlacement(Rot2A, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell+1), assembly_log, false, id+5);
        
        // place probes
        
        r = 61*mm/2;
        new G4PVPlacement(Rot1B, G4ThreeVector(r*sin(angle1), -r*cos(angle1), dz+40*mm+11*mm), target_probe_t_log_, Form("target_probe_t%d", 2*icell), assembly_log, false, id+6);
        new G4PVPlacement(Rot2B, G4ThreeVector(-r*sin(angle2), -r*cos(angle2), dz-40*mm-11*mm), target_probe_p_log_, Form("target_probe_p%d", 2*icell), assembly_log, false, id+6);
        
    };
    
    //
    G4Colour target_color;
    G4Colour::GetColour("target", target_color);
    G4Colour target_si_color;
    G4Colour::GetColour("target_si", target_si_color);
    G4Colour target_base_color;
    G4Colour::GetColour("target_base", target_base_color);
    G4Colour target_film_color;
    G4Colour::GetColour("target_film", target_film_color);
    G4Colour target_pipe_color;
    G4Colour::GetColour("target_pipe", target_pipe_color);
    assembly_log->SetVisAttributes(G4VisAttributes::Invisible);
    target_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_base_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
    target_cell_film_log_->SetVisAttributes(new G4VisAttributes(target_film_color));
    target_nipple_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    target_pipe_h_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    target_pipe_v_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    target_probe_p_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    target_probe_t_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    
    target_log_->SetUserLimits(target_step_limit_);
    target_base_log_->SetUserLimits(target_step_limit_);
    
    detector_parts->AttachBiasingOperator(target_log_);
    detector_parts->AttachBiasingOperator(target_base_log_);
    detector_parts->AttachBiasingOperator(target_cell_film_log_);
    detector_parts->AttachBiasingOperator(target_nipple_log_);
    detector_parts->AttachBiasingOperator(target_pipe_h_log_);
    detector_parts->AttachBiasingOperator(target_pipe_v_log_);
    detector_parts->AttachBiasingOperator(target_probe_p_log_);
    detector_parts->AttachBiasingOperator(target_probe_t_log_);

}


void g4PSITargetJeru::SetSD(G4SDManager *SDman) {
    
    // if the target is placed into a target chamber, the sensitive detector
    // for the target is the one of the chamber if it is found
    
    G4String TargetSDname = "g4PSI/" + (is_in_chamber() ? target_chamber_->GetDetectorName() : label_);
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(TargetSDname);
    if (sd == NULL) {
        sd = target_sd_ = new g4PSITargetSD(TargetSDname, label_ + "_Collection", min_p_, min_theta_);
        SDman->AddNewDetector( sd );
    };
    
    if (target_log_) target_log_->SetSensitiveDetector( sd );
    if (target_log_) target_log_->SetSensitiveDetector( sd );
    if (target_cell_film_log_) target_cell_film_log_->SetSensitiveDetector( sd );
    if (target_base_log_) target_base_log_->SetSensitiveDetector( sd );
    if (target_nipple_log_) target_nipple_log_->SetSensitiveDetector( sd );
    if (target_pipe_v_log_) target_pipe_v_log_->SetSensitiveDetector( sd );
    if (target_pipe_h_log_) target_pipe_h_log_->SetSensitiveDetector( sd );
    if (target_probe_t_log_) target_probe_t_log_->SetSensitiveDetector( sd );
    if (target_probe_p_log_) target_probe_p_log_->SetSensitiveDetector( sd );
}


void g4PSITargetJeru::Write() {
    TVectorD v1(6);
    v1[0] = g4PSIDetectorParts::getInstance()->GetTargetState();
    v1.Write(label_);
}


void g4PSITargetJeru::InitTree(TTree *T) {
    if (target_sd_) target_sd_->InitTree(T);
}


void g4PSITargetJeru::DeleteEventData() {
    if (target_sd_) target_sd_->DeleteEventData();
}

