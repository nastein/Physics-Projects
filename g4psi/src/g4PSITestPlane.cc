#include "g4PSITestPlane.hh"
#include "g4PSITrackerSD.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"   
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSITestPlane::g4PSITestPlane(G4String label, G4double angle, 
			       G4double r, G4double x, G4double y,
			       G4LogicalVolume *mv) : g4PSIDetectorBase(label)  {
    if (mv) {
	this->SetMotherVolume(mv);
    };
    testplane_label_ = label;
    testplane_active_x_ = x;
    testplane_active_y_ = y;
    testplane_angle_ = angle;
    testplane_r_ = r;
}

g4PSITestPlane::~g4PSITestPlane() {
}


void g4PSITestPlane::Info() {
    InfoTitle("Testplane");
    InfoParDouble("Position at r", testplane_r_/cm, "cm");
    InfoPar2Double("Active area", testplane_active_x_/cm, testplane_active_y_/cm, "cm");
    InfoParDouble("Angle", testplane_angle_/deg, "deg");
}

void g4PSITestPlane::Placement() {

    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");  

    G4Colour testplane_colour;
    G4Colour::GetColour("testplane", testplane_colour);

    G4double testplane_z =  0.1*mm;
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateY(-testplane_angle_);
    G4Box* testplane_detector = new G4Box("testplane_detector", testplane_active_x_/2, testplane_active_y_/2, testplane_z/2);
  
    detector_log_ = new G4LogicalVolume(testplane_detector, Air, (testplane_label_+"_log").c_str(), 0, 0, 0);
    new G4PVPlacement(zRot, G4ThreeVector(0, 0, testplane_r_ * cos(testplane_angle_)),
		      detector_log_, (testplane_label_+"_phys").c_str(), mother_volume_, false, 0,
		      false);
    
    G4VisAttributes* vis_attributes = new G4VisAttributes(testplane_colour);
    vis_attributes->SetForceWireframe(true);
    detector_log_->SetVisAttributes(vis_attributes);

    detector_log_->SetVisAttributes (G4VisAttributes::Invisible);
}


void g4PSITestPlane::SetSD(G4SDManager *SDman) {
    G4String TestPlaneSDname = "g4PSI/TestPlane/" + testplane_label_;
    G4VSensitiveDetector* TestPlaneSD = SDman->FindSensitiveDetector(TestPlaneSDname);
    if (TestPlaneSD == NULL) {
	TestPlaneSD = SD_ = new g4PSITrackerSD(TestPlaneSDname, testplane_label_ + "_Collection", true);
	SDman->AddNewDetector( TestPlaneSD );
    };
    detector_log_->SetSensitiveDetector( TestPlaneSD );
}


void g4PSITestPlane::Write() {
    TVectorD v1(4);
    v1[0] = testplane_r_;
    v1[1] = testplane_angle_;
    v1[2] = testplane_active_x_;
    v1[3] = testplane_active_y_;
    v1.Write(testplane_label_);
}

void g4PSITestPlane::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSITestPlane::DeleteEventData() {
    SD_->DeleteEventData();
}
