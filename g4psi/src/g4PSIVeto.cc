#include "g4PSIVeto.hh"
#include "g4PSIScintillatorSD.hh"
#include "g4PSIDetectorParts.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSIVeto::g4PSIVeto(G4String label,
                         G4double width, G4double z, G4int nSeg, G4double r1, G4double r2, G4int nPlane)
: g4PSIDetectorBase(label) {
    sc_label_ = label;
    sc_width_ = width;
    sc_z_ = z;
    sc_r1_ = r1;
    sc_r2_ = r2;
    sc_no_segments_ = nSeg;
    sc_no_planes_ = nPlane;
}


g4PSIVeto::~g4PSIVeto() {
}



void g4PSIVeto::Info() {
    InfoTitle("Veto Detector");
    InfoParInt("Number of segments", sc_no_segments_);
    InfoParInt("Number of planes", sc_no_planes_);
    InfoParDouble("Detector (center) z", sc_z_/cm, "cm");
    InfoParDouble("Detector thickness dz", sc_width_/cm, "cm");
    InfoParDouble("r_min", sc_r1_/cm, "cm");
    InfoParDouble("r_max", sc_r2_/cm, "cm");
//    G4cout << "Scintillator Ring:" << G4endl;
//    G4cout << "=================" << G4endl;
//    G4cout << G4endl;
//    G4cout << "                                     name = " << sc_label_ << G4endl;
//    G4cout << "                                        n = " << sc_no_segments_ << " pieces" << G4endl;
//    G4cout << "                                        z = " << sc_z_/cm << " cm" << G4endl;
//    G4cout << "                                       dz = " << sc_width_/cm << " cm" << G4endl;
//    G4cout << "                                    r_min = " << sc_r1_/cm << " cm" << G4endl;
//    G4cout << "                                    r_max = " << sc_r2_/cm << " cm" << G4endl;
//    G4cout << G4endl;
}

double g4PSIVeto::GetZPos() {return sc_z_;}

void g4PSIVeto::Placement() {
    
    // Material
    G4Material* Scinti = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    if (sc_no_segments_ == 0) {
        // circular ring detector
        G4Tubs* sc_tub = new G4Tubs((sc_label_+"_geom").c_str(), sc_r1_, sc_r2_, sc_width_/2., 0, twopi);
        sc_log_ = new G4LogicalVolume(sc_tub, Scinti, (sc_label_+"_log").c_str() , 0,0,0);
        new G4PVPlacement(0,G4ThreeVector(0.,0., sc_z_ - sc_width_/2.), sc_log_, (sc_label_+"_phys").c_str(),
                          mother_volume_,false,0);
        
        G4Colour sc_col;
        G4Colour::GetColour("scwall", sc_col);
        sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
    } //else {
    // multi-trapezoidal ring detector
    //G4double phiStart = 0;
    //G4double phiTotal = twopi;
    //G4int numSide = sc_no_segments_;
    //const G4int numZPlanes = 2;
    //G4double rInner[numZPlanes] = {sc_r1_, sc_r1_};
    //G4double rOuter[numZPlanes] = {sc_r2_, sc_r2_};
    //G4double z[numZPlanes] = {0, sc_width_};
    //G4Polyhedra* sc_tub = new G4Polyhedra((sc_label_+"_geom").c_str(), phiStart, phiTotal,
    //numSide, numZPlanes, z, rInner, rOuter);
    //sc_log_ = new G4LogicalVolume(sc_tub, Scinti, (sc_label_+"_log").c_str() , 0,0,0);
    //new G4PVPlacement(0,G4ThreeVector(0.,0., sc_z_ - sc_width_/2.), sc_log_, (sc_label_+"_phys").c_str(),
    //mother_volume_,false,0);
    
    //G4Colour sc_col;
    //G4Colour::GetColour("scwall", sc_col);
    //sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
    //}
    else {
        ////Trapezoids (Wedges)
        double dx1 = (sc_r1_*sin(pi/sc_no_segments_));                        //half-length along x at bottom surface
        double dx2 = (sc_r2_*sin(pi/sc_no_segments_));                        //half-length along x at top surface
        double dy1 = sc_width_/2;                                   //half-length of wedge thickness along bottom surface
        double dy2 = sc_width_/2;                                   //half-length of wedge thickness along top surface
        double height = ((sc_r2_-sc_r1_)*(cos(pi/sc_no_segments_)))/2;        //half-length of height of wedge
        double rshift = 1*mm;
        
        
        G4Trd* sc_trd = new G4Trd((sc_label_+"_geom").c_str(), dx1, dx2, dy1, dy2, height);
        sc_log_ = new G4LogicalVolume(sc_trd, Scinti, (sc_label_+"_log").c_str(),0,0,0);
        G4Colour sc_col;
        G4Colour::GetColour("scwall", sc_col);
        sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
        
        int iCopy = 0;
        for (int i = 0; i < sc_no_segments_; i++) {
            for (int p = 0; p < sc_no_planes_; p++) {
                char cName[100];
                sprintf(cName, "%s_%d", label_.c_str(), iCopy+1);
                double n = i + 0.5 + (double) p / sc_no_planes_;
                double anglePos = (pi/sc_no_segments_ + n*(2*pi/sc_no_segments_));
                double angleRot = -((pi/2)-(pi/sc_no_segments_))+ n*(2*pi/sc_no_segments_);
                G4RotationMatrix* Rot = new G4RotationMatrix;
                Rot->rotateX(pi/2);
                Rot->rotateY(angleRot);
                double xpos = ((sc_r1_*cos(pi/sc_no_segments_) + height + rshift) * (cos(anglePos)));
                double ypos = ((sc_r1_*cos(pi/sc_no_segments_) + height + rshift) * (sin(anglePos)));
                new G4PVPlacement(Rot, G4ThreeVector(xpos, ypos, sc_z_ + p*sc_width_), sc_log_, (sc_label_+"_phys").c_str(), mother_volume_, 0, iCopy);
                iCopy++;
            }
        }
    }
}


void g4PSIVeto::SetSD(G4SDManager *SDman) {
    G4String RingSCSDname = "g4PSI/RingSC/" + sc_label_;
    G4VSensitiveDetector* RingSCSD = SDman->FindSensitiveDetector(RingSCSDname);
    if (RingSCSD == NULL) {
        RingSCSD = SD_ = new g4PSIScintillatorSD(RingSCSDname, sc_no_segments_ * sc_no_planes_, sc_label_ + "_Collection" );
        SDman->AddNewDetector( RingSCSD );
    };
    sc_log_->SetSensitiveDetector( RingSCSD );
}


void g4PSIVeto::Write() {
    TVectorD v1(4);
    v1[0] = sc_z_;
    v1[1] = sc_width_;
    v1[2] = sc_r1_;
    v1[3] = sc_r2_;
    v1.Write(sc_label_);
}

void g4PSIVeto::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSIVeto::DeleteEventData() {
    SD_->DeleteEventData();
}
