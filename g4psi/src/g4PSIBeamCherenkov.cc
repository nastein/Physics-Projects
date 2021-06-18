#include "g4PSIBeamCherenkov.hh"
#include "g4PSIScintillatorSD.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Trap.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSIBeamCherenkov::g4PSIBeamCherenkov(G4String label, G4double z) : g4PSIDetectorBase(label) {
    bcc_z_ = z;
    bcc_N1_ = 30;
    bcc_len_ = 80*mm;
    bcc_paddle_width1_ = 3*mm;
    bcc_paddle_thickness_ = 3*mm;
    bcc_angle_yaw_ = 45*deg;
    bcc_gap_ = 0;
    bcc_angle_roll_ = 0;
    bcc_mode_ = g4PSIBeamCherenkov::kBCC_vMode;
    bcc_mat_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    bcc_log1_ = NULL;
    bcc_log2_ = NULL;
    bcc_sipm1_log_ = NULL;
    bcc_sipm2_log_ = NULL;
}

g4PSIBeamCherenkov::g4PSIBeamCherenkov(G4String label, G4Material *mat, G4double z,
                                       G4double angle, G4double len,
                                       G4double sx, G4double sy) : g4PSIDetectorBase(label) {
    bcc_z_ = z;
    bcc_N1_ = 1;
    bcc_len_ = len;
    bcc_paddle_width1_ = sx;
    bcc_paddle_thickness_ = sy;
    bcc_angle_yaw_ = angle;
    bcc_mat_ = mat;
    bcc_gap_ = 0;
    bcc_angle_roll_ = 0;
    bcc_mode_ = g4PSIBeamCherenkov::kBCC_SinglePaddle;
    bcc_log1_ = NULL;
    bcc_log2_ = NULL;
    bcc_sipm1_log_ = NULL;
    bcc_sipm2_log_ = NULL;
}

g4PSIBeamCherenkov::g4PSIBeamCherenkov(G4String label, G4Material *mat, G4double z,
                                       G4double angle, G4double len,
                                       G4int N1, G4double width1, G4int N2, G4double width2, G4double thickness,
                                       G4double rot, G4double gap) : g4PSIDetectorBase(label) {
    bcc_z_ = z;
    bcc_N1_ = N1;
    bcc_N2_ = N2;
    bcc_len_ = len;
    bcc_paddle_width1_ = width1;
    bcc_paddle_width2_ = width2;
    bcc_paddle_thickness_ = thickness;
    bcc_angle_yaw_ = angle;
    bcc_mat_ = mat;
    bcc_gap_ = gap;
    bcc_angle_roll_ = rot;
    if (N2 == 0) {
        bcc_mode_ = g4PSIBeamCherenkov::kBCC_MultiplePaddles;
    } else {
        bcc_mode_ = g4PSIBeamCherenkov::kBCC_SISH;
    }
    bcc_log1_ = NULL;
    bcc_log2_ = NULL;
    bcc_sipm1_log_ = NULL;
    bcc_sipm2_log_ = NULL;
}


g4PSIBeamCherenkov::~g4PSIBeamCherenkov() {
}


double g4PSIBeamCherenkov::GetZPos() {return bcc_z_;}

void g4PSIBeamCherenkov::Info() {
    
    InfoTitle("Beam Cherenkov Counter");
    InfoParString("Material", bcc_mat_->GetName());
    InfoParDouble("Angle", bcc_angle_yaw_/deg, "deg");
    InfoParDouble("Position z", bcc_z_/cm, "cm");
    if (bcc_mode_ == kBCC_vMode) {
        InfoParString("Design", "vMode");
        InfoParInt("Number of fibers", bcc_N1_);
        InfoParDouble("Size x", bcc_paddle_width1_/cm, "cm");
        InfoParDouble("Size y", bcc_paddle_thickness_/cm, "cm");
        InfoParDouble("Length", bcc_len_/cm, "cm");
    } else if (bcc_mode_ == kBCC_SISH) {
        InfoParString("Design", "Multiple Planes");
        InfoParInt("Number of inner paddles", bcc_N1_);
        InfoParInt("Number of outer paddles", bcc_N2_);
        InfoParDouble("Inner paddle width", bcc_paddle_width1_/mm, "mm");
        InfoParDouble("Outer paddle width", bcc_paddle_width2_/mm, "mm");
        InfoParDouble("Paddle length", bcc_len_/mm, "mm");
        InfoParDouble("Paddle thickness", bcc_paddle_thickness_/mm, "mm");
        InfoParDouble("Gap between paddles", bcc_gap_/mm, "mm");
        InfoParDouble("Plane orientation (roll)", bcc_angle_roll_/deg, "deg");
    } else if (bcc_mode_ == kBCC_MultiplePaddles) {
        InfoParString("Design", "Multiple Planes");
        InfoParInt("Number of segments", bcc_N1_);
        InfoParDouble("Segment width", bcc_paddle_width1_/cm, "cm");
        InfoParDouble("Segment length", bcc_len_/cm, "cm");
        InfoParDouble("Segment thickness", bcc_paddle_thickness_/cm, "cm");
        InfoParDouble("Gap between segments", bcc_gap_/mm, "mm");
        InfoParDouble("Plane orientation (roll)", bcc_angle_roll_/deg, "deg");
    } else if (bcc_mode_ == kBCC_SinglePaddle) {
        InfoParString("Design", "Single Plane");
        InfoParDouble("Size x", bcc_paddle_width1_/cm, "cm");
        InfoParDouble("Size y", bcc_paddle_thickness_/cm, "cm");
        InfoParDouble("Length", bcc_len_/cm, "cm");
    } else {
        InfoParString("Design", "unknown");
    }
}


void g4PSIBeamCherenkov::Placement() {
    
    // Material
    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* silicon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4Material* aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    if (bcc_mode_ == g4PSIBeamCherenkov::kBCC_vMode) {
        G4double delta = bcc_paddle_width1_ / tan(bcc_angle_yaw_);
        
        G4Trap* bcc_group = new G4Trap("bcc_group", bcc_N1_*bcc_paddle_width1_, bcc_paddle_width1_, bcc_len_, bcc_len_ - delta);
        G4LogicalVolume *bcc_up_log = new G4LogicalVolume(bcc_group, air, (label_+"_up_log").c_str(), 0, 0, 0);
        G4LogicalVolume *bcc_down_log = new G4LogicalVolume(bcc_group, air, (label_+"_down_log").c_str(), 0, 0, 0);
        
        // Place lower half of the BCC
        G4RotationMatrix* zRot1 = new G4RotationMatrix;
        zRot1->rotateY(90*deg);
        zRot1->rotateZ(180*deg + bcc_angle_yaw_);
        zRot1->rotateX(180*deg);
        new G4PVPlacement(zRot1, G4ThreeVector(0, -(bcc_len_  - delta/2)*sin(bcc_angle_yaw_) / 2, bcc_z_),
                          bcc_down_log, (label_+"_lower_phys").c_str(), mother_volume_, false, 0);
        
        // Place upper half of the BCC
        G4RotationMatrix* zRot2 = new G4RotationMatrix;
        zRot2->rotateY(90*deg);
        zRot2->rotateZ(180*deg - bcc_angle_yaw_);
        new G4PVPlacement(zRot2, G4ThreeVector(0, +(bcc_len_  - delta/2)*sin(bcc_angle_yaw_) / 2, bcc_z_),
                          bcc_up_log, (label_+"_upper_phys").c_str(), mother_volume_, false, 0);
        
        bcc_up_log->SetVisAttributes (G4VisAttributes::Invisible);
        bcc_down_log->SetVisAttributes (G4VisAttributes::Invisible);
        
        // Place single BCC crystals
        G4Trap* bcc_crystal = new G4Trap("bcc_crystal", bcc_paddle_width1_, bcc_paddle_width1_, bcc_len_, bcc_len_ - delta);
        bcc_log1_ = new G4LogicalVolume(bcc_crystal, bcc_mat_, (label_+"_crystal_log").c_str(), 0, 0, 0);
        for (G4int i = 1; i <= 2*bcc_N1_; i++) {
            char pName[60];
            sprintf(pName, "%s_%02d", label_.c_str(), i);
            new G4PVPlacement(0, G4ThreeVector(0,0,((i-1)/2-bcc_N1_/2)*bcc_paddle_width1_), bcc_log1_, pName, (i-1)%2 == 0 ? bcc_up_log : bcc_down_log, false, i-1);
        }
    } else if (bcc_mode_ == g4PSIBeamCherenkov::kBCC_MultiplePaddles) {
        
        // \todo fill gap with aluminized mylar and air; presently solid aluminum
        
        G4double bcc_group_dx = (bcc_N1_*bcc_paddle_width1_ + (bcc_N1_-1) * bcc_gap_)/2 + bcc_gap_;
        G4double bcc_group_dy = fmax(bcc_paddle_thickness_/2, sipm_width_/2) + bcc_gap_;
        G4double bcc_group_dz = bcc_len_/2 + sipm_height_ + bcc_gap_;
        G4Box* bcc_group = new G4Box("bcc_group", bcc_group_dx, bcc_group_dy, bcc_group_dz);
        G4LogicalVolume *bcc_log = new G4LogicalVolume(bcc_group, aluminum, (label_+"_log").c_str(), 0, 0, 0);
        
        G4RotationMatrix* zRot1 = new G4RotationMatrix;
        zRot1->rotateY(90.*deg);
        zRot1->rotateX(90.*deg - bcc_angle_roll_);  // bcc_angle_roll_ is the +rotation of the detector about the z axis.
        zRot1->rotateZ(bcc_angle_yaw_);
        new G4PVPlacement(zRot1, G4ThreeVector(0, 0, bcc_z_),
                          bcc_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
        
        //        G4VisAttributes* logVisAtt = new G4VisAttributes();
        //        logVisAtt->SetForceWireframe(true);
        //        bcc_log->SetVisAttributes (logVisAtt);
        bcc_log->SetVisAttributes (G4VisAttributes::Invisible);
        
        //        // Place single BCC crystals
        G4Box* bcc_crystal = new G4Box("bcc_crystal", bcc_paddle_width1_/2, bcc_paddle_thickness_/2, bcc_len_/2);
        bcc_log1_ = new G4LogicalVolume(bcc_crystal, bcc_mat_, (label_+"_crystal_log").c_str(), 0, 0, 0);
        G4Box* bcc_sipm = new G4Box("bcc_sipm", bcc_paddle_width1_/2, sipm_width_/2, sipm_height_/2);
        bcc_sipm1_log_ = new G4LogicalVolume(bcc_sipm, silicon, (label_+"_sipm_log").c_str(), 0, 0, 0);
        G4Colour bcc_col;
        G4Colour::GetColour("pmt", bcc_col);
        bcc_sipm1_log_->SetVisAttributes(new G4VisAttributes(bcc_col));
        for (G4int i = 1; i <= bcc_N1_; i++) {
            char pName[60];
            sprintf(pName, "%s_%02d", label_.c_str(), i);
            G4double x = (2*i-1-bcc_N1_)*(bcc_paddle_width1_+bcc_gap_)/2.;
            new G4PVPlacement(0, G4ThreeVector(x,0,0), bcc_log1_, pName, bcc_log, false, i-1);
            G4double y = bcc_len_/2+sipm_height_/2;
            sprintf(pName, "%s_sipm1_%02d", label_.c_str(), i);
            new G4PVPlacement(0, G4ThreeVector(x,0,+y), bcc_sipm1_log_, pName, bcc_log, false, bcc_N1_ + i-1);
            sprintf(pName, "%s_sipm2_%02d", label_.c_str(), i);
            new G4PVPlacement(0, G4ThreeVector(x,0,-y), bcc_sipm1_log_, pName, bcc_log, false, 2*bcc_N1_ + i-1);
        }
    } else if (bcc_mode_ == g4PSIBeamCherenkov::kBCC_SISH) {
        
        // \todo fill gap with aluminized mylar and air; presently air
        
        G4double bcc_group_dx = (bcc_N1_*bcc_paddle_width1_ + (bcc_N1_-1) * bcc_gap_)/2 + bcc_N2_*bcc_paddle_width2_ + (bcc_N2_-1) * bcc_gap_;
        if (bcc_N1_ > 0) bcc_group_dx += bcc_gap_;
        if (bcc_N2_ > 0) bcc_group_dx += bcc_gap_;
        G4double bcc_group_dy = fmax(bcc_paddle_thickness_/2, sipm_width_/2) + bcc_gap_;
        G4double bcc_group_dz = bcc_len_/2 + sipm_height_ + bcc_gap_;
        G4Box* bcc_group = new G4Box("bcc_group", bcc_group_dx, bcc_group_dy, bcc_group_dz);
        G4LogicalVolume *bcc_log = new G4LogicalVolume(bcc_group, air, (label_+"_log").c_str(), 0, 0, 0);
        
        G4RotationMatrix* zRot1 = new G4RotationMatrix;
        zRot1->rotateY(90.*deg);
        zRot1->rotateX(90.*deg - bcc_angle_roll_);  // bcc_angle_roll_ is the +rotation of the detector about the z axis.
        zRot1->rotateZ(bcc_angle_yaw_);
        new G4PVPlacement(zRot1, G4ThreeVector(0, 0, bcc_z_),
                          bcc_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
        
//        G4VisAttributes* logVisAtt = new G4VisAttributes();
//        logVisAtt->SetForceWireframe(true);
//        bcc_log->SetVisAttributes (logVisAtt);
        bcc_log->SetVisAttributes (G4VisAttributes::Invisible);

        
        G4Box* bcc_crystal1 = new G4Box("bcc_crystal1", bcc_paddle_width1_/2, bcc_paddle_thickness_/2, bcc_len_/2);
        bcc_log1_ = new G4LogicalVolume(bcc_crystal1, bcc_mat_, (label_+"_crystal1_log").c_str(), 0, 0, 0);
        G4Box* bcc_crystal2 = new G4Box("bcc_crystal2", bcc_paddle_width2_/2, bcc_paddle_thickness_/2, bcc_len_/2);
        bcc_log2_ = new G4LogicalVolume(bcc_crystal2, bcc_mat_, (label_+"_crystal2_log").c_str(), 0, 0, 0);
        G4Box* bcc_sipm1 = new G4Box("bcc_sipm1", bcc_paddle_width1_/2, sipm_width_/2, sipm_height_/2);
        bcc_sipm1_log_ = new G4LogicalVolume(bcc_sipm1, silicon, (label_+"_sipm1_log").c_str(), 0, 0, 0);
        G4Box* bcc_sipm2 = new G4Box("bcc_sipm2", bcc_paddle_width2_/2, sipm_width_/2, sipm_height_/2);
        bcc_sipm2_log_ = new G4LogicalVolume(bcc_sipm2, silicon, (label_+"_sipm2_log").c_str(), 0, 0, 0);
        G4Colour bcc_col;
        G4Colour::GetColour("pmt", bcc_col);
        bcc_sipm1_log_->SetVisAttributes(new G4VisAttributes(bcc_col));
        bcc_sipm2_log_->SetVisAttributes(new G4VisAttributes(bcc_col));
        
        // copy ID:
        //
        //   XX : inner paddle
        //   YY : outer paddle
        //
        //          0    1  N2-1  N2  N2+1 N2+N1-1 N1+N2 ... N1+2*N2-1
        //         XX   XX   XX   YY   YY   YY   XX   XX   XX
        //    -------------------------+-------------------------------> X
        
        
        // Place first outer paddles
        
        int copyID = 0;
        double w1 = bcc_N1_ * bcc_paddle_width1_ + (bcc_N1_ - 1) * bcc_gap_;
        double w2 = bcc_N2_ * bcc_paddle_width2_ + (bcc_N2_ - 1) * bcc_gap_;
        double x = -(w1 / 2. + 3*bcc_gap_/2. + w2);
        double step1 = (bcc_paddle_width1_ + bcc_gap_) / 2.0;
        double step2 = (bcc_paddle_width2_ + bcc_gap_) / 2.0;
        for (G4int i = 1; i <= bcc_N2_; i++) {
            place_paddle(bcc_log2_, bcc_sipm2_log_, bcc_log, x, copyID, step2);
        }
        for (G4int i = 1; i <= bcc_N1_; i++) {
            place_paddle(bcc_log1_, bcc_sipm1_log_, bcc_log, x, copyID, step1);
        }
        for (G4int i = 1; i <= bcc_N2_; i++) {
            place_paddle(bcc_log2_, bcc_sipm2_log_, bcc_log, x, copyID, step2);
        }
                
    } else if (bcc_mode_ == g4PSIBeamCherenkov::kBCC_SinglePaddle) {
        G4Box* bcc_group = new G4Box("bcc_group", bcc_paddle_width1_/2, bcc_paddle_thickness_/2, bcc_len_/2);
        bcc_log1_ = new G4LogicalVolume(bcc_group, bcc_mat_, (label_+"_log").c_str(), 0, 0, 0);
        
        G4RotationMatrix* zRot1 = new G4RotationMatrix;
        //	zRot1->rotateY(90*deg);
        zRot1->rotateY(bcc_angle_yaw_);
        zRot1->rotateZ(90*deg);
        new G4PVPlacement(zRot1, G4ThreeVector(0, 0, bcc_z_),
                          bcc_log1_, (label_+"_phys").c_str(), mother_volume_, false, 0);
    } else {
        
    }
    
    G4Colour bcc_col;
    G4Colour::GetColour("bcc", bcc_col);
    bcc_log1_->SetVisAttributes(new G4VisAttributes(bcc_col));
    bcc_log2_->SetVisAttributes(new G4VisAttributes(bcc_col));
}


void g4PSIBeamCherenkov::SetSD(G4SDManager *SDman) {
    G4int nCells = 0;
    if (bcc_log1_) nCells += bcc_N1_;
    if (bcc_log2_) nCells += 2*bcc_N2_;
    if (bcc_sipm1_log_) nCells += 2*bcc_N1_;
    if (bcc_sipm2_log_) nCells += 4*bcc_N2_;
    if (nCells > 0) {
        G4String BeamCherenkovSDname = "g4PSI/BeamCherenkov/" + label_;
        G4VSensitiveDetector* bcc_SD = SDman->FindSensitiveDetector(BeamCherenkovSDname);
        if (bcc_SD == NULL) {
            bcc_SD = SD_ = new g4PSIScintillatorSD(BeamCherenkovSDname, nCells, label_ + "_Collection" );
            SDman->AddNewDetector( bcc_SD );
        };
        if (bcc_log1_) bcc_log1_->SetSensitiveDetector( bcc_SD );
        if (bcc_log2_) bcc_log2_->SetSensitiveDetector( bcc_SD );
        if (bcc_sipm1_log_) bcc_sipm1_log_->SetSensitiveDetector( bcc_SD );
        if (bcc_sipm2_log_) bcc_sipm2_log_->SetSensitiveDetector( bcc_SD );
    }
}


void g4PSIBeamCherenkov::Write() {
    TVectorD v1(12);
    v1[0] = bcc_z_;
    v1[1] = bcc_N1_;
    v1[2] = bcc_len_;
    v1[3] = bcc_paddle_width1_;
    v1[4] = bcc_paddle_thickness_;
    v1[5] = bcc_angle_yaw_;
    v1[6] = bcc_angle_roll_;
    v1[7] = bcc_gap_;
    v1[8] = 0;
    v1[9] = bcc_mode_;
    v1[10] = bcc_N2_;
    v1[11] = bcc_paddle_width2_;
    v1.Write(label_);
}

void g4PSIBeamCherenkov::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSIBeamCherenkov::DeleteEventData() {
    SD_->DeleteEventData();
}


void g4PSIBeamCherenkov::place_paddle(G4LogicalVolume *paddle, G4LogicalVolume* sipm, G4LogicalVolume *mother, double &x, int &copyID, double step) {
    char pName[60];
    sprintf(pName, "%s_%02d", label_.c_str(), copyID+1);

    x += step;
    new G4PVPlacement(0, G4ThreeVector(x,0,0), paddle, pName, mother, false, copyID);
    //    std::cout << "Placing " << pName << " at x = " << x << " with copyID " << copyID << "\n";

    G4double y = bcc_len_/2+sipm_height_/2;
    sprintf(pName, "%s_sipm1_%02d", label_.c_str(), copyID+1);
    new G4PVPlacement(0, G4ThreeVector(x,0,+y), sipm, pName, mother, false, bcc_N1_+2*bcc_N2_ + copyID);
    sprintf(pName, "%s_sipm2_%02d", label_.c_str(), copyID+1);
    new G4PVPlacement(0, G4ThreeVector(x,0,-y), sipm, pName, mother, false, 2*bcc_N1_+2*bcc_N2_ + copyID);
    x += step;
    copyID++;
}
