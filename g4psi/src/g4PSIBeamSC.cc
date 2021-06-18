#include "g4PSIBeamSC.hh"
#include "g4PSIScintillatorSD.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSIBeamSC::g4PSIBeamSC(G4String name,
                         G4double angle,
                         G4double fiber_sz,
                         G4double fiber_cladding,
                         G4double fiber_spacing,
                         G4double wx,
                         G4double wy,
                         G4double offsetx,
                         G4double offsety,
                         G4double posz)
: g4PSIDetectorBase(name) {
    
    label_ = name;
    angle_ = angle;
    fiber_sz_ = fabs(fiber_sz);
    fiber_cladding_ = fiber_cladding;
    fiber_spacing_ = fiber_spacing;
    wx_ = wx;
    wy_ = wy;
    offsetx_ = offsetx;
    offsety_ = offsety;
    posz_ = posz;
    
    circular_fiber_ = fiber_sz > 0.0;
    
    NFibers_ = ceil( (wx + fiber_spacing_) / (fiber_sz_ + fiber_spacing_) );
    wx_ = NFibers_ * (fiber_sz_ + fiber_spacing_) - fiber_spacing_;
}


g4PSIBeamSC::g4PSIBeamSC(G4String name, G4double angle,
                         G4double fiber_sz, G4double wx,
                         G4double wy, G4double posz, G4double offset)
: g4PSIDetectorBase(name) {
    
    label_ = name;
    angle_ = angle;
    fiber_sz_ = fabs(fiber_sz);
    fiber_cladding_ = 0.0;
    fiber_spacing_ = 0.0;
    wx_ = wx;
    wy_ = wy;
    offsetx_ = offset;
    offsety_ = 0.0;
    posz_ = posz;
    
    circular_fiber_ = fiber_sz > 0.0;
    
    NFibers_ = ceil( (wx + fiber_spacing_) / (fiber_sz_ + fiber_spacing_) );
    wx_ = NFibers_ * (fiber_sz_ + fiber_spacing_) - fiber_spacing_;
}

g4PSIBeamSC::~g4PSIBeamSC() {
    
}


double g4PSIBeamSC::GetZPos() {return posz_;}

void g4PSIBeamSC::Info() {
    InfoTitle("Beam Scintillating Fiber");
    InfoParInt("Number of fibers", NFibers_);
    InfoParString("Fiber Shape", (circular_fiber_ ? "circular" : "square"));
    if (circular_fiber_) {
        InfoParDouble("Fiber diameter (incl. cladding)", fiber_sz_/mm, "mm");
        InfoParDouble("Radial thickness of fiber cladding", fiber_cladding_/mm, "mm");
        InfoPar2Double("Physical size (diameter x length)", fiber_sz_/cm, wy_/cm, "cm");
    } else {
        InfoPar3Double("Physical size", fiber_sz_/cm, fiber_sz_/cm, wy_/cm, "cm");
    }
    InfoParDouble("Detector angle", angle_/deg, "deg");
    InfoParDouble("Detector width x", wx_/cm, "cm");
    InfoParDouble("Detector width y", wy_/cm, "cm");
    InfoParDouble("Detector offset x", offsetx_/mm, "mm");
    InfoParDouble("Detector offset y", offsety_/mm, "mm");
    InfoParDouble("Position z", posz_/cm, "cm");
    
    
//    G4cout << "Beam Scintillating Fiber:" << G4endl;
//    G4cout << "========================" << G4endl;
//    G4cout << G4endl;
//    G4cout << "                                     name = " << label_ << G4endl;
//    G4cout << "                         number of fibers = " << NFibers_ << G4endl;
//    G4cout << "                              fiber shape = " << (circular_fiber_ ? "circular" : "square") << G4endl;
//    if (circular_fiber_) {
//        G4cout << "          fiber diameter (incl. cladding) = " << fiber_sz_/mm << " mm" << G4endl;
//        G4cout << "           radial thickness of fiber cladding = " << fiber_cladding_/mm << " mm" << G4endl;
//        
//        G4cout << "                          physical length = "
//	       << fiber_sz_/cm << " cm (diam.) x " << wy_/cm << " cm" << G4endl;
//    } else {
//        G4cout << "                          physical length = "
//	       << fiber_sz_/cm << " cm x " << fiber_sz_/cm << " cm x " << wy_/cm << " cm" << G4endl;
//    }
//    G4cout << "                           detector angle = " << angle_/deg << " deg" << G4endl;
//    G4cout << "                         detector width x = " << wx_/cm << " cm" << G4endl;
//    G4cout << "                         detector width y = " << wy_/cm << " cm" << G4endl;
//    G4cout << "                        detector offset x = " << offsetx_/cm << " cm" << G4endl;
//    G4cout << "                        detector offset y = " << offsety_/cm << " cm" << G4endl;
//    G4cout << "                               Position z = " << posz_/cm << " cm" << G4endl;
//    G4cout << G4endl;
    
}



void g4PSIBeamSC::Placement() {
    
    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* BCF10Fiber = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");
    
    //------------------------- Detector
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateZ(angle_);
    G4Box* plane = new G4Box((label_+"_box").c_str(), wx_/2, wy_/2, fiber_sz_/2);
    G4LogicalVolume *plane_log = new G4LogicalVolume(plane, Air, (label_+"_log").c_str(), 0, 0, 0);
    new G4PVPlacement(zRot,
                      G4ThreeVector(offsetx_, offsety_, posz_),
                      plane_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
    plane_log->SetVisAttributes (G4VisAttributes::Invisible);
    // G4VisAttributes* logVisAtt = new G4VisAttributes();
    // logVisAtt->SetForceWireframe(true);
    // plane_log->SetVisAttributes(logVisAtt);
    
    //------------------------- SC Fibers
    G4VSolid* sc_fiber = NULL;
    G4VSolid* sc_cladding = NULL;
    if (circular_fiber_) {
        sc_fiber = new G4Tubs("sc_fiber", 0, fiber_sz_/2. - fiber_cladding_, wy_/2, 0, twopi);
        sc_cladding = new G4Tubs("sc_cladding",
                                 fiber_sz_/2. - fiber_cladding_, fiber_sz_/2., wy_/2, 0, twopi);
    } else {
        sc_fiber = new G4Box("sc_fiber", fiber_sz_/2. - fiber_cladding_,
                             fiber_sz_/2. - fiber_cladding_, wy_/2.);
    };
    
    sc_log_ = new G4LogicalVolume(sc_fiber, BCF10Fiber, (label_+"_fib_log").c_str(), 0, 0, 0);
    G4LogicalVolume *sc_cladding_log = NULL;
    if (sc_cladding) {
        sc_cladding_log = new G4LogicalVolume(sc_cladding, BCF10Fiber, (label_+"_cladding_log").c_str(), 0, 0, 0);
    };
    G4RotationMatrix* zRotFib = new G4RotationMatrix;
    zRotFib->rotateX(90*deg);
    for (G4int i = 0; i < NFibers_; i++) {
        char pName[60];
        sprintf(pName, "%s_%02d_phys", label_.c_str(), i + 1);
        new G4PVPlacement(zRotFib,
                          G4ThreeVector(fiber_sz_/2. + i*(fiber_sz_ + fiber_spacing_) - wx_ / 2., 0, 0),
                          sc_log_, pName, plane_log, false, i);
        if (sc_cladding_log) {
            sprintf(pName, "%s_%02d_clad_phys", label_.c_str(), i + 1);
            new G4PVPlacement(zRotFib,
                              G4ThreeVector(fiber_sz_/2. + i*(fiber_sz_ + fiber_spacing_) - wx_ / 2., 0, 0),
                              sc_cladding_log, pName, plane_log, false, i);
        }
    }
    
    G4Colour sc_col;
    G4Colour::GetColour("scifi", sc_col);
    sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
    if (sc_cladding_log) {
        G4Colour::GetColour("scifi_cladding", sc_col);
        sc_cladding_log->SetVisAttributes(new G4VisAttributes(sc_col));
    }
}


void g4PSIBeamSC::SetSD(G4SDManager *SDman) {
    if (NFibers_ >= 0) {
        G4String BeamSCSDname = "g4PSI/BeamSC/" + label_;
        G4VSensitiveDetector* BeamSCSD = SDman->FindSensitiveDetector(BeamSCSDname);
        if (BeamSCSD == NULL) {
            BeamSCSD = SD_ = new g4PSIScintillatorSD(BeamSCSDname, NFibers_, label_ + "_Collection" );
            SDman->AddNewDetector( BeamSCSD );
        };
        sc_log_->SetSensitiveDetector( BeamSCSD );
    }
}


void g4PSIBeamSC::Write() {
    TVectorD v1(10);
    v1[0] = NFibers_;
    v1[1] = fiber_sz_;
    v1[2] = fiber_spacing_;
    v1[3] = wx_;
    v1[4] = wy_;
    v1[5] = offsetx_;
    v1[6] = offsety_;
    v1[7] = posz_;
    v1[8] = angle_;
    v1[9] = circular_fiber_ ? 1 : 0;
    v1.Write(label_);
}

void g4PSIBeamSC::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSIBeamSC::DeleteEventData() {
    SD_->DeleteEventData();
}
