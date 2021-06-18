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
#include "g4PSITargetUM.hh"
#include "TVectorD.h"


using namespace CLHEP;


g4PSITargetUM::g4PSITargetUM(G4String label, g4PSIChamberBase* chamber) : g4PSITargetBase(label, chamber) {
    
    target_sd_ = NULL;
    
    tgt_R_ = 30*mm;
    tgt_H_ = 110*mm;
    base_r1_ = 6.38*cm/2;
    base_r2_ = 60*mm/2 - 0.15*cm;
    base_r3_ = 60*mm/2;
    base_h1_ = 5*mm;
    base_h3_ = 10*mm;
    
    tgt_wall_t_ = 120*um;
    tgt_cell_gap_ = 1*cm;
    tgt_cell_height_ = tgt_H_ + 2*(base_h1_ + base_h3_);
    
    radial_stub_od_ = 0.313*cm;
    radial_stub_t_ = 0.81*mm;

    support_tube_cu_od_ = 0.635*cm;
    support_tube_cu_t_ = 0.51*mm;
    support_tube_ss_od_ = 0.635*cm;
    support_tube_ss_t_ = 0.51*mm;
    
    cu_tube_od_ = 1*cm;
    cu_tube_t_ = 1*mm;
    cu_tube_l_ = 21.17*cm;
    
    support_tube_cu_l_ = cu_tube_l_ + tgt_cell_height_ + tgt_cell_gap_ + radial_stub_od_;
    support_tube_ss_l_ = 14*cm;
    radial_stub_l_ = 5.842*cm - support_tube_cu_od_/2 - base_r1_;

    insulation_t_ = .03125*mm;  //~5 layers of superinsulation

    
    base_material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
    tgt_wall_material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
}


void g4PSITargetUM::Info() {
    InfoTitle("Target (cylinder)");
    InfoParDouble("Target-cylinder height", tgt_H_ / mm, " mm");
    InfoParDouble("Target-cylinder radius", tgt_R_/ mm, " mm");
    InfoParDouble("Target wall thickness", tgt_wall_t_ / um, " um");
    InfoParString("Target wall material", tgt_wall_material_->GetName());
    InfoParString("Target cap material", base_material_->GetName());
    InfoParDouble("Target insulation thickness", insulation_t_/ mm, " mm");
}


void g4PSITargetUM::Placement() {
    
    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    G4Material* Copper = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
    G4Material* StainlessSteel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
    G4Material* AlMylar = new G4Material("AluminizedMylar", 1.762 * g/cm3, 2);
    AlMylar->AddMaterial(Mylar, 0.7217);
    AlMylar->AddMaterial(Al, 0.2783);
    
    G4Material* empty_mat = GetMaterialForTargetMotherLV();
    if (!empty_mat) {
        empty_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    }
    
    G4Material* target_mat = NULL;
    if (detector_parts->GetTargetState() == g4PSIDetectorParts::kFull) {
        G4cout << "Full Target: Filling target with liquid hydrogen\n" << G4endl;
        target_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
    } else if (detector_parts->GetTargetState() == g4PSIDetectorParts::kEmpty) {
        G4cout << "Empty Target: Filling target with hydryogen gas\n" << G4endl;
        target_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
    };
    std::cout << "Target density: " << target_mat->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;
    
    // target volume
    
    G4double target_full_height = tgt_H_/2. + base_h3_;
    G4double target_inner[] = {0, 0, 0, 0, 0, 0};
    G4double target_outer[] = {base_r2_, base_r2_, tgt_R_, tgt_R_, base_r2_, base_r2_};
    G4double target_z[] = {-target_full_height, -tgt_H_/2., -tgt_H_/2., tgt_H_/2., tgt_H_/2., target_full_height};
    G4Polycone* target_volume = new G4Polycone("target_cell_volume", 0, 2*pi, 6, target_z, target_inner, target_outer);
    
    target_log_ = new G4LogicalVolume(target_volume, target_mat, "target_cell_volume.log", 0, 0, 0);
    
    // base
    
    G4Tubs* lip = new G4Tubs("target_cell_base_1", base_r2_, base_r3_, base_h3_/2, 0., twopi);
    G4Tubs* cap = new G4Tubs("target_cell_base_3",        0, base_r1_, base_h1_/2, 0., twopi);
    
    G4VSolid* base = lip;
    base = new G4UnionSolid("target_cell_base_top", base, cap, NULL, G4ThreeVector(0., 0., (base_h1_+base_h3_)*mm/2));
    
    target_base_log_ = new G4LogicalVolume(base, base_material_, "target_cell_base.log", 0, 0, 0);
    
    // film
    
    G4Tubs* target_cell_film = new G4Tubs("target_cell_film", tgt_R_, tgt_R_ + tgt_wall_t_, tgt_H_/2., 0, 2*pi);
    target_cell_film_log_ = new G4LogicalVolume(target_cell_film, tgt_wall_material_, "target_cell_film.log", 0, 0, 0);
    
    // radial stubs
    
    G4Tubs* radial_stub = new G4Tubs("radial_stub", radial_stub_od_/2-radial_stub_t_, radial_stub_od_/2, radial_stub_l_/2, 0, 2*pi);
    radial_stub_cu_log_ = new G4LogicalVolume(radial_stub, Copper, "radial_stub_cu_log");
    radial_stub_ss_log_ = new G4LogicalVolume(radial_stub, StainlessSteel, "radial_stub_ss_log");
    
    // support tubes
    
    G4Tubs* support_tube_ss = new G4Tubs("support_tube_ss", support_tube_ss_od_/2-support_tube_ss_t_, support_tube_ss_od_/2, support_tube_ss_l_/2, 0, 2*pi);
    support_tube_ss_log_ = new G4LogicalVolume(support_tube_ss, StainlessSteel, "support_tube_ss_log");
    G4Tubs* support_tube_cu = new G4Tubs("support_tube_cu", support_tube_cu_od_/2-support_tube_cu_t_, support_tube_cu_od_/2, support_tube_cu_l_/2, 0, 2*pi);
    support_tube_cu_log_ = new G4LogicalVolume(support_tube_cu, Copper, "support_tube_cu_log");

    G4Tubs* cu_tube = new G4Tubs("support_tube_cu", cu_tube_od_/2-cu_tube_t_, cu_tube_od_/2, cu_tube_l_/2, 0, 2*pi);
    cu_tube_log_ = new G4LogicalVolume(cu_tube, Copper, "support_tube_cu_log");

    // super insulation
    
    G4Trd *Trd1 = new G4Trd("Trd", 2*cm, 5.5*cm, 8*cm, 8*cm, 4.5*cm);
    G4Trd *Trd2 = new G4Trd("Trd", 2*cm - insulation_t_, 5.5*cm - insulation_t_, 9*cm, 9*cm, 4.5*cm - insulation_t_);
    G4SubtractionSolid* Trd = new G4SubtractionSolid("Trd", Trd1, Trd2, NULL, G4ThreeVector(0., 0., 0.));
    target_si_log_ = new G4LogicalVolume(Trd, AlMylar, "target_si", 0,0,0);
    

    // logical volume for entire target assembly
    
    G4double assembly_r = 90 * mm;
    G4double assembly_h = 270 * mm;
    G4double assembly_shift = -30*mm;
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
    
    new G4PVPlacement(rot_assembly, (*chamber_rot).invert() * G4ThreeVector(0, -assembly_shift, -1*mm) + chamber_shift, assembly_log, "target_assembly", mother_lv, false, 0);
    
    const G4int NCELL = 2;
    G4double angle1 = 45*deg + 90*deg;
    G4double angle2 = 45*deg - 90*deg;
    G4RotationMatrix* Rot = new G4RotationMatrix;
    Rot->rotateX(180*deg);
    Rot->rotateZ(180*deg);
    G4RotationMatrix* Rot1 = new G4RotationMatrix;
    Rot1->rotateY(90*deg);
    Rot1->rotateX(90*deg + 45*deg + 90*deg);
    G4RotationMatrix* Rot2 = new G4RotationMatrix;
    Rot2->rotateY(90*deg);
    Rot2->rotateX(90*deg + 45*deg);

    for (G4int icell = 0; icell < NCELL; icell++) {
        G4double dz = -icell * (tgt_cell_height_ + tgt_cell_gap_) + assembly_shift;
        G4int id = (icell + 1) * 10; // target cell copy IDs are 10 and 20
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_log_,
                          Form("target_cell%d_volume", icell), assembly_log, false, id);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_cell_film_log_,
                          Form("target_cell%d_film", icell), assembly_log, false, id+1);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,dz+tgt_H_/2+base_h3_/2), target_base_log_,
                          Form("target_cell%d_base", icell), assembly_log, false, id+2);
        new G4PVPlacement(Rot, G4ThreeVector(0,0, dz-tgt_H_/2-base_h3_/2), target_base_log_,
                          Form("target_cell%d_base", icell), assembly_log, false, id+2);
        
        G4double r = base_r1_ + radial_stub_l_ / 2;
        new G4PVPlacement(Rot1, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+tgt_cell_height_/2-radial_stub_od_/2), radial_stub_cu_log_, Form("radial_stub_cu%d", 2*icell), assembly_log, false, id+5);
        new G4PVPlacement(Rot2, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz+tgt_cell_height_/2-radial_stub_od_/2), radial_stub_cu_log_, Form("radial_stub_cu%d", 2*icell+1), assembly_log, false, id+5);
    };
    
    
    G4double r = base_r1_ + radial_stub_l_ + support_tube_cu_od_/2;
    G4double z = assembly_shift + tgt_cell_height_/2 + cu_tube_l_/2;
    new G4PVPlacement(NULL, G4ThreeVector(0,0,z), cu_tube_log_, "cu_tube", assembly_log, false, 30);
    
    z = assembly_shift-tgt_cell_height_/2-tgt_cell_gap_-radial_stub_od_+ support_tube_cu_l_/2;
    new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), z), support_tube_cu_log_, "support_tube_cu1", assembly_log, false, 31);
    new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), z), support_tube_cu_log_, "support_tube_cu2", assembly_log, false, 31);
    
    z = assembly_shift - tgt_cell_height_/2 - tgt_cell_gap_ - radial_stub_od_ - support_tube_ss_l_/2;
    new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), z), support_tube_ss_log_, "support_tube_ss1", assembly_log, false, 32);
    new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), z), support_tube_ss_log_, "support_tube_ss2", assembly_log, false, 32);
    

    G4RotationMatrix* RotAl = new G4RotationMatrix;
    RotAl->rotateX(-90*deg);
    RotAl->rotateY(-90*deg);
    RotAl->rotateZ(0*deg);
    new G4PVPlacement(RotAl, G4ThreeVector(0,0, assembly_shift), target_si_log_, "target_insulation", assembly_log, false, 33);

    //
    G4Colour target_color;
    G4Colour::GetColour("target", target_color);
    G4Colour target_si_color;
    G4Colour::GetColour("target_si", target_si_color);
    G4Colour copper_color;
    G4Colour::GetColour("copper", copper_color);
    G4Colour kapton_color;
    G4Colour::GetColour("kapton", kapton_color);
    G4Colour ss_color;
    G4Colour::GetColour("aluminum", ss_color);
    
    assembly_log->SetVisAttributes(G4VisAttributes::Invisible);
    target_log_->SetVisAttributes(new G4VisAttributes(target_color));
    target_cell_film_log_->SetVisAttributes(new G4VisAttributes(kapton_color));
    target_base_log_->SetVisAttributes(new G4VisAttributes(copper_color));
    target_si_log_->SetVisAttributes(new G4VisAttributes(target_si_color));
    radial_stub_cu_log_->SetVisAttributes(new G4VisAttributes(copper_color));
    support_tube_cu_log_->SetVisAttributes(new G4VisAttributes(copper_color));
    support_tube_ss_log_->SetVisAttributes(new G4VisAttributes(ss_color));
    cu_tube_log_->SetVisAttributes(new G4VisAttributes(copper_color));
    
    
//    target_nipple_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
//    target_pipe_h_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
//    target_pipe_v_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
//    target_probe_p_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
//    target_probe_t_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
    
    assembly_log->SetUserLimits(target_step_limit_);
    target_log_->SetUserLimits(target_step_limit_);
    target_base_log_->SetUserLimits(target_step_limit_);
    
    detector_parts->AttachBiasingOperator(target_log_);
    detector_parts->AttachBiasingOperator(target_base_log_);
    detector_parts->AttachBiasingOperator(target_cell_film_log_);
    detector_parts->AttachBiasingOperator(target_si_log_);
    detector_parts->AttachBiasingOperator(radial_stub_cu_log_);
    detector_parts->AttachBiasingOperator(support_tube_cu_log_);
    detector_parts->AttachBiasingOperator(cu_tube_log_);

//    detector_parts->AttachBiasingOperator(target_nipple_log_);
//    detector_parts->AttachBiasingOperator(target_pipe_h_log_);
//    detector_parts->AttachBiasingOperator(target_pipe_v_log_);
//    detector_parts->AttachBiasingOperator(target_probe_p_log_);
//    detector_parts->AttachBiasingOperator(target_probe_t_log_);

}


void g4PSITargetUM::SetSD(G4SDManager *SDman) {
    
    // if the target is placed into a target chamber, the sensitive detector
    // for the target is the one of the chamber if it is found
    
    G4String TargetSDname = "g4PSI/" + (is_in_chamber() ? target_chamber_->GetDetectorName() : label_);
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(TargetSDname);
    if (sd == NULL) {
        sd = target_sd_ = new g4PSITargetSD(TargetSDname, label_ + "_Collection", min_p_, min_theta_);
        SDman->AddNewDetector( sd );
    };
    
    if (target_log_) target_log_->SetSensitiveDetector( sd );
    if (target_cell_film_log_) target_cell_film_log_->SetSensitiveDetector( sd );
    if (target_base_log_) target_base_log_->SetSensitiveDetector( sd );
    if (target_si_log_) target_si_log_->SetSensitiveDetector( sd );
    if (radial_stub_cu_log_) radial_stub_cu_log_->SetSensitiveDetector( sd );
    if (support_tube_cu_log_) support_tube_cu_log_->SetSensitiveDetector( sd );
    if (cu_tube_log_) cu_tube_log_->SetSensitiveDetector( sd );
}


void g4PSITargetUM::Write() {
    TVectorD v1(6);
    v1[0] = g4PSIDetectorParts::getInstance()->GetTargetState();
    v1.Write(label_);
}


void g4PSITargetUM::InitTree(TTree *T) {
    if (target_sd_) target_sd_->InitTree(T);
}


void g4PSITargetUM::DeleteEventData() {
    if (target_sd_) target_sd_->DeleteEventData();
}

