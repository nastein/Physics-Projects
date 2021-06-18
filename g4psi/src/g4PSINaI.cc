#include "g4PSINaI.hh"
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

g4PSINaI::g4PSINaI(G4String label, G4double z, G4double t)
: g4PSIDetectorBase(label) {
    nai_label_ = label;
    nai_z_ = z;
    nai_n_ = 7;
    nai_t_pb_ = t;
    
    const G4double inch = 2.54 * cm;
    bar_x_ = 16 * inch;
    bar_y_ =  2 * inch;
    bar_z_ =  4 * inch;
}


g4PSINaI::~g4PSINaI() {
}


void g4PSINaI::Info() {
    InfoTitle("NaI Detector");
    InfoParDouble("Detector (front face) z", nai_z_/cm, "cm");
}

double g4PSINaI::GetZPos() {return nai_z_;}

void g4PSINaI::Placement() {
    
    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* NaI = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    G4Material* Lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

    // Detector group
    


    G4double box_x = bar_x_;
    G4double box_y = nai_n_ * bar_y_;
    G4double box_z = bar_z_ + nai_t_pb_;

    G4Box *box = new G4Box((label_+"_box").c_str(), box_x/2., box_y/2., box_z/2.);
    G4LogicalVolume *box_log = new G4LogicalVolume(box, Air, (nai_label_+"_box_log").c_str(), 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(0,0,nai_z_ + box_z/2.), box_log, (nai_label_+"_box_phys").c_str(), mother_volume_, false, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    VisAtt->SetForceWireframe(true);
    box_log->SetVisAttributes(VisAtt);
    box_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    // NaI bars
    
    G4Box *nai = new G4Box((label_+"_nai").c_str(), bar_x_/2., bar_y_/2., bar_z_/2.);
    nai_log_ = new G4LogicalVolume(nai, NaI, (nai_label_+"_nai_log").c_str(), 0, 0, 0);
    for (int i = 0; i < nai_n_; i++) {
        G4double y = (i + 0.5 - nai_n_/2.0)*bar_y_;
        new G4PVPlacement(0, G4ThreeVector(0,y, nai_t_pb_/2.), nai_log_, (nai_label_+"_nai_phys").c_str(), box_log, false, i);
    }
    G4Colour col;
    G4Colour::GetColour("scwall", col);
    G4VisAttributes* nai_vis = new G4VisAttributes(col);
    nai_log_->SetVisAttributes(nai_vis);
    
    // Pb shield
    
    if (nai_t_pb_ > 0) {
        G4Box *pb = new G4Box((label_+"_pb").c_str(), bar_x_/2., box_y/2., nai_t_pb_/2.);
        G4LogicalVolume *pb_log = new G4LogicalVolume(pb, Lead, (nai_label_+"_pb_log").c_str(), 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0,0,-bar_z_/2.), pb_log, (nai_label_+"_pb_phys").c_str(), box_log, false, 0);
        G4Colour::GetColour("lead", col);
        G4VisAttributes* pb_vis = new G4VisAttributes(col);
        pb_log->SetVisAttributes(pb_vis);
    }
}


void g4PSINaI::SetSD(G4SDManager *SDman) {
    G4String RingSCSDname = "g4PSI/NaI/" + nai_label_;
    G4VSensitiveDetector* NaISD = SDman->FindSensitiveDetector(RingSCSDname);
    if (NaISD == NULL) {
        NaISD = SD_ = new g4PSIScintillatorSD(RingSCSDname, nai_n_, nai_label_ + "_Collection" );
        SDman->AddNewDetector( NaISD );
    };
    nai_log_->SetSensitiveDetector( NaISD );
}


void g4PSINaI::Write() {
    TVectorD v1(6);
    v1[0] = nai_z_;
    v1[1] = nai_n_;
    v1[2] = nai_t_pb_;
    v1[3] = bar_x_;
    v1[4] = bar_y_;
    v1[5] = bar_z_;
    v1.Write(nai_label_);
}

void g4PSINaI::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSINaI::DeleteEventData() {
    SD_->DeleteEventData();
}
