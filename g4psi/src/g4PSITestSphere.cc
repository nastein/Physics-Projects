//
//  g4PSITestSphere.cpp
//  g4PSI
//
//  Created by Steffen Strauch on 6/5/15.
//
//

#include "g4PSITestSphere.hh"
#include "g4PSITrackerSD.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSITestSphere::g4PSITestSphere(G4String label,
                                 G4double r,
                                 G4double p_threshold,
                                 G4LogicalVolume *mv) : g4PSIDetectorBase(label)  {
    if (mv) {
        this->SetMotherVolume(mv);
    };
    r_ = r;
    p_threshold_ = p_threshold;
}

g4PSITestSphere::~g4PSITestSphere() {
}


void g4PSITestSphere::Info() {
    G4cout << "Testsphere:" << G4endl;
    G4cout << "=========" << G4endl;
    G4cout << G4endl;
    G4cout << "                                     name = " << label_ << G4endl;
    G4cout << "                                   radius = " << r_/cm << " cm" << G4endl;
    G4cout << G4endl;
    
}

void g4PSITestSphere::Placement() {
    
    // Material
    G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    G4Colour testsphere_colour;
    G4Colour::GetColour("testplane", testsphere_colour);
    
    G4double rmin = r_;
    G4double rmax = rmin + 0.1 * CLHEP::mm;
    
    G4Sphere* testplane_detector = new G4Sphere("testsphere_detector", rmin, rmax, 0, 360*CLHEP::deg, 0*CLHEP::deg, 180*CLHEP::deg);
    
    testsphere_assembly_log_ = new G4LogicalVolume(testplane_detector, vacuum, (label_+"_log").c_str(), 0, 0, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0),
                      testsphere_assembly_log_, (label_+"_phys").c_str(), mother_volume_, false, 0,
                      false);
    G4VisAttributes* testplane_assembly_logVisAtt = new G4VisAttributes(testsphere_colour);
    testplane_assembly_logVisAtt->SetForceWireframe(true);
    testsphere_assembly_log_->SetVisAttributes(testplane_assembly_logVisAtt);
    
    G4Colour col;
    G4Colour::GetColour("stc_frame", col);
    G4VisAttributes* vis_att = new G4VisAttributes(col);
    testsphere_assembly_log_->SetVisAttributes (vis_att);
    
    //    testsphere_assembly_log_->SetVisAttributes (G4VisAttributes::Invisible);
    testsphere_detector_log_ = testsphere_assembly_log_;
}


void g4PSITestSphere::SetSD(G4SDManager *SDman) {
    G4String TestSphereSDname = "g4PSI/TestSphere/" + label_;
    G4VSensitiveDetector* TestSphereSD = SDman->FindSensitiveDetector(TestSphereSDname);
    if (TestSphereSD == NULL) {
        TestSphereSD = SD_ = new g4PSITrackerSD(TestSphereSDname, label_ + "_Collection", true, p_threshold_);
        SDman->AddNewDetector( TestSphereSD );
    };
    testsphere_assembly_log_->SetSensitiveDetector( TestSphereSD );
}


void g4PSITestSphere::Write() {
    TVectorD v1(1);
    v1[0] = r_;
    v1.Write(label_);
}

void g4PSITestSphere::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSITestSphere::DeleteEventData() {
    SD_->DeleteEventData();
}
