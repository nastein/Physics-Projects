#include "g4PSICal.hh"
#include "g4PSIScintillatorSD.hh"
#include "g4PSIDetectorParts.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSICal::g4PSICal(G4String label, CalType ct, G4double z, G4double t_pb, G4double t_sc, G4int n)
: g4PSIDetectorBase(label) {
    Cal_label_ = label;
    Cal_z_ = z;

    cal_sc_log_ = NULL;
    cal_pb_log_ = NULL;
    
    cal_type_ = ct;
    
    switch (cal_type_) {
        case kPbSC:
            cal_nx_ = 1;
            cal_ny_ = 1;
            cal_nz_ = n;
            cal_t_pb_ = t_pb;
            cal_t_sc_ = t_sc;
            cal_dx_ = 50 * cm;
            cal_dy_ = 50 * cm;
            cal_mat_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");
            break;
            
        case kPbGlass:
            cal_nx_ = 20;
            cal_ny_ = 20;
            cal_nz_ = n;
            cal_t_pb_ = t_pb;
            cal_t_sc_ = t_sc;
            cal_dx_ = 2 * cm;
            cal_dy_ = 2 * cm;
            cal_mat_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_LEAD");
            break;
            
        default:
            break;
    }
    // 2 cm x 2 cm x 20 cm - LdGlass
    
}


g4PSICal::~g4PSICal() {
}


void g4PSICal::Info() {
    InfoTitle("Cal Detector");
    InfoParDouble("Detector (front face) z", Cal_z_/cm, "cm");
}

double g4PSICal::GetZPos() {return Cal_z_;}

void g4PSICal::Placement() {
    
    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* Lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

    // Detector group
    
    G4double box_x = cal_nx_ * cal_dx_;
    G4double box_y = cal_ny_ * cal_dy_;
    G4double box_z = cal_nz_ * (cal_t_pb_ + cal_t_sc_);

    G4Box *box = new G4Box((label_+"_box").c_str(), box_x/2., box_y/2., box_z/2.);
    G4LogicalVolume *box_log = new G4LogicalVolume(box, Air, (Cal_label_+"_box_log").c_str(), 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(0,0,Cal_z_ + box_z/2.), box_log, (Cal_label_+"_box_phys").c_str(), mother_volume_, false, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    VisAtt->SetForceWireframe(true);
    box_log->SetVisAttributes(VisAtt);
    box_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    // calorimeter pb + sc sheets
    
    G4Box *cal_sc = new G4Box((label_+"_cal_sc").c_str(), cal_dx_/2., cal_dy_/2.,cal_t_sc_/2.);
    cal_sc_log_ = new G4LogicalVolume(cal_sc, cal_mat_, (Cal_label_+"_cal_sc_log").c_str(), 0, 0, 0);
    if (cal_t_pb_ > 0) {
        G4Box *cal_pb = new G4Box((label_+"_cal_pb").c_str(), cal_dx_/2., cal_dy_/2., cal_t_pb_/2.);
        cal_pb_log_ = new G4LogicalVolume(cal_pb, Lead, (Cal_label_+"_cal_pb_log").c_str(), 0, 0, 0);
    } else {
        cal_pb_log_ = NULL;
    }
    int ncopy = 0;
    for (int ix = 0; ix < cal_nx_; ix++) {
        for (int iy = 0; iy < cal_ny_; iy++) {
            for (int iz = 0; iz < cal_nz_; iz++) {
                G4double x = -box_x/2. + (ix + 0.5) * cal_dx_;
                G4double y = -box_y/2. + (iy + 0.5) * cal_dy_;
                G4double z = -box_z/2. + iz * (cal_t_pb_ + cal_t_sc_) + cal_t_pb_ / 2.;
                if (cal_pb_log_)
                    new G4PVPlacement(0, G4ThreeVector(x, y, z), cal_pb_log_, (Cal_label_+"_cal_pb_phys").c_str(), box_log, false, ncopy);
                z = -box_z/2. + iz * (cal_t_pb_ + cal_t_sc_) + cal_t_pb_ + cal_t_sc_/2.;
                new G4PVPlacement(0, G4ThreeVector(x, y, z), cal_sc_log_, (Cal_label_+"_cal_sc_phys").c_str(), box_log, false, ncopy);
                ncopy++;
            }
        }
    }
    G4Colour col;
    G4Colour::GetColour("scwall", col);
    G4VisAttributes* cal_sc_vis = new G4VisAttributes(col);
    cal_sc_log_->SetVisAttributes(cal_sc_vis);
    if (cal_pb_log_) {
        G4Colour::GetColour("lead", col);
        G4VisAttributes* cal_pb_vis = new G4VisAttributes(col);
        cal_pb_log_->SetVisAttributes(cal_pb_vis);
    }
}


void g4PSICal::SetSD(G4SDManager *SDman) {
    G4String RingSCSDname = "g4PSI/Cal/" + Cal_label_;
    G4VSensitiveDetector* CalSD = SDman->FindSensitiveDetector(RingSCSDname);
    if (CalSD == NULL) {
        CalSD = SD_ = new g4PSIScintillatorSD(RingSCSDname, cal_nx_ * cal_ny_ * cal_nz_, Cal_label_ + "_Collection" );
        SDman->AddNewDetector( CalSD );
    };
    cal_sc_log_->SetSensitiveDetector( CalSD );
}


void g4PSICal::Write() {
    TVectorD v1(8);
    v1[0] = Cal_z_;
    v1[1] = cal_nx_;
    v1[2] = cal_ny_;
    v1[3] = cal_nz_;
    v1[4] = cal_t_pb_;
    v1[5] = cal_t_sc_;
    v1[6] = cal_dx_;
    v1[7] = cal_dy_;
    v1.Write(Cal_label_);
}

void g4PSICal::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSICal::DeleteEventData() {
    SD_->DeleteEventData();
}
