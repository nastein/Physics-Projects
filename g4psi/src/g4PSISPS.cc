#include "g4PSISPS.hh"
#include "g4PSIScintillatorSD.hh"
#include "g4PSIDetectorParts.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"
#include "G4UserLimits.hh"

using namespace CLHEP;


g4PSISPS::g4PSISPS(G4String label,
                   SPSWallType type,
                   G4double width,
                   G4double gap_small,
                   G4double gap_wide,
                   G4double thickness,
                   G4double height,
                   G4double distance,
                   G4double theta,
                   G4double shift,
                   G4int Nunits) : g4PSIDetectorBase(label) {
    
    wall_type_ = type;
    sc_width_ = width;
    sc_thickness_ = thickness;
    sc_gap_small_ = gap_small;
    sc_gap_wide_ = gap_wide;
    sc_distance_to_front_ = distance;
    sc_bar_length_ = height;
    upstream_shift_ = shift;
    
    sc_theta_center_ = theta;
    
    sc_Nunits_ = Nunits;
    
    const G4double inch = 2.54 * cm;

    if (wall_type_ == front) {
        backing_length_ = 166.8 * cm;
        backing_al_length_ = 23.4 * cm + 4 * inch; // check
        backing_al_thickness_ = 1./8. * inch;
        backing_al_width_ = 4. * inch;
        backing_rc_thickness_ = 1.5 * cm;
        backing_rc_width_ = 4 * inch;
        backing_rc_length_ = 118.0 * cm;
    } else if (wall_type_ == rear) {
        backing_length_ = 266.8 * cm;
        backing_al_length_ = backing_length_;
        backing_al_thickness_ = 1./8. * inch;
        backing_al_width_ = 4. * inch;
        backing_rc_thickness_ = 0 * cm;
        backing_rc_width_ = 0;
        backing_rc_length_ = 0 * cm;
    } else {
        G4Exception("g4PSISPS::Placement", "unknown wall type", FatalException, "");
    }
    
    sps_wall_dy_ = sc_thickness_ + backing_al_thickness_ + backing_rc_thickness_;
    sps_wall_dz_ = backing_length_;
    sps_wall_dx_ = sc_Nunits_ * (2 * sc_width_ + sc_gap_small_) + (sc_Nunits_ + 1) * sc_gap_wide_;
    
    sps_r_rear_of_backing_ = sc_distance_to_front_ + sps_wall_dy_;
    
    sps_w1_ = sps_wall_dx_ / 2 + upstream_shift_;
    sps_w2_ = sps_wall_dx_ / 2 - upstream_shift_;
    
    sps_rot_ = new G4RotationMatrix;
    sps_rot_->rotateX(90*deg);
    sps_rot_->rotateZ(-sc_theta_center_);
    sign_ = (sc_theta_center_ > 0) ? 1.0 : -1.0;
}


g4PSISPS::~g4PSISPS() {
}


double g4PSISPS::GetZPos() {return 0;}

void g4PSISPS::Info() {
    InfoTitle("Scintillator Wall");
    InfoParInt("Number of scintillator bars ", 2*sc_Nunits_);
    InfoParInt("Nunit", sc_Nunits_);
    InfoPar3Double("Physical length", sc_width_/cm, sc_thickness_/cm, sc_bar_length_/cm, "cm");
    InfoParDouble("Front-face distance to the target center", sc_distance_to_front_/cm, "cm");
    InfoParDouble("Gapwidth between scintillators inside units", sc_gap_small_/cm, "cm");
    InfoParDouble("Gapwidth between scintillator units", sc_gap_wide_/cm, "cm");
    InfoParDouble("Detector width", (sps_w1_+sps_w2_)/cm, "cm");
    InfoParDouble(" -- upstream part", sps_w2_/cm, "cm");
    InfoParDouble(" -- downstream part", sps_w1_/cm, "cm");
}


void g4PSISPS::Placement() {
    
    // Material
    /// \todo don't use Air, use the material of the world volume
    
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* Glass = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");
    G4Material* SCMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    // Place a full wall of scintillator bars
    // --------------------------------------



    //------------------------- Wall
    //    x0 and z0 are the coordinates of the center of the wall including
    //    scintillator and backing structures, not frames
    //    construct and place container
    
    //  N=3 units: *[].[]*[].[]*[].[]*
    
    G4double x0 = (sps_r_rear_of_backing_ - sps_wall_dy_ / 2.0) * sin(sc_theta_center_) + sign_ * upstream_shift_ * cos(sc_theta_center_);
    G4double z0 = (sps_r_rear_of_backing_ - sps_wall_dy_ / 2.0) * cos(sc_theta_center_) - sign_ * upstream_shift_ * sin(sc_theta_center_);
    
    G4Box* sc_wall = new G4Box((label_+"sc_wall").c_str(), sps_wall_dx_/2, sps_wall_dy_/2, sps_wall_dz_/2);
    G4LogicalVolume *sc_wall_log = new G4LogicalVolume(sc_wall, Air, (label_+"_m1_log").c_str(), 0, 0, 0);
    new G4PVPlacement(sps_rot_, G4ThreeVector(x0, 0, z0), sc_wall_log, (label_+"_phys").c_str(), mother_volume_, false, 0);

    G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(1,0,0));
    VisAtt->SetForceWireframe(true);
    sc_wall_log->SetVisAttributes(VisAtt);
    sc_wall_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    
    //------------------------- SC Bars
    
    G4Box* sc_bar = new G4Box("sc", sc_width_/2, sc_thickness_/2, sc_bar_length_/2);
    sc_log_ = new G4LogicalVolume(sc_bar, SCMat, (label_+"_bar_log").c_str(), 0, 0, 0);
    for (G4int i = 0; i < 2 * sc_Nunits_; i++) {
        char pName[60];
        sprintf(pName, "%s%02d", label_.c_str(), i + 1);
        //  N=3 units: *[].[]*[].[]*[].[]*
        //               0  1  2  3  4  5
        
        G4double x = -sps_wall_dx_/2. + sc_gap_wide_ + sc_width_ * (i + 0.5) + sc_gap_wide_ * int(i/2) + sc_gap_small_ * int((i+1)/2);
        new G4PVPlacement(0, G4ThreeVector(sign_*x,sps_wall_dy_/2-sc_thickness_/2,0), sc_log_, pName, sc_wall_log, false, i);
    }
    G4Colour sc_col;
    G4Colour::GetColour("scwall", sc_col);
    sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
    
    //------------------------- PMTs

    G4LogicalVolume *pmt_log = NULL;
    
    const double pmt_height = 134*mm;
    const double pmt_width = 52.6*mm;
    
    G4Colour pmt_col;
    G4Colour::GetColour("pmt", pmt_col);
    G4Tubs* pmt = new G4Tubs((label_ + "_pmt_geom").c_str(), 0, pmt_width/2, pmt_height/2, 0, 360*deg);
    pmt_log = new G4LogicalVolume(pmt, Glass, (label_ + "_pmt_log").c_str(), 0,0,0);
    pmt_log->SetVisAttributes(new G4VisAttributes(pmt_col));
    
    for (G4int i = 0; i < 2 * sc_Nunits_; i++) {
        char pName[60];
        
        G4double x = -sps_wall_dx_/2. + sc_gap_wide_ + sc_width_ * (i + 0.5) + sc_gap_wide_ * int(i/2) + sc_gap_small_ * int((i+1)/2);
        
        sprintf(pName, "%s_%02d_PMT_up", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(sign_*x,sps_wall_dy_/2-sc_thickness_/2, +(sc_bar_length_+pmt_height)/2.),
                          pmt_log, pName, sc_wall_log, false, 100+i);
        sprintf(pName, "%s_%02d_PMT_down", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(sign_*x,sps_wall_dy_/2-sc_thickness_/2, -(sc_bar_length_+pmt_height)/2.),
                          pmt_log, pName, sc_wall_log, false, 200+i);
    }
    
    //------------------------- Backing structure
    
    G4Box* al_back = new G4Box((label_+"al_back").c_str(), backing_al_width_/2, backing_al_thickness_/2, backing_al_length_/2);
    G4LogicalVolume *al_back_log = new G4LogicalVolume(al_back, Air, (label_+"_al_back_log").c_str(), 0, 0, 0);
    G4Box* rc_back = NULL;
    G4LogicalVolume* rc_back_log = NULL;
    if (wall_type_ == front) {
        rc_back = new G4Box((label_+"rc_back").c_str(), backing_rc_width_/2, backing_rc_thickness_/2, backing_rc_length_/2);
        rc_back_log = new G4LogicalVolume(rc_back, Air, (label_+"_rc_back_log").c_str(), 0, 0, 0);
    }
    for (G4int i = 0; i < sc_Nunits_; i++) {
        char pName[60];
        sprintf(pName, "%s_Backing%02d", label_.c_str(), i + 1);
        G4double x = -sps_wall_dx_/2. + sc_gap_wide_ + sc_width_ + sc_gap_small_ * 0.5 + (sc_gap_small_ + sc_gap_wide_ + 2*sc_width_) * i;
        G4double y = -sps_wall_dy_/2+backing_al_thickness_/2;
        if (wall_type_ == front) {
            G4double z = (backing_length_ - backing_al_length_) / 2.;
            G4double yrc = -sps_wall_dy_/2 + backing_al_thickness_ + backing_rc_thickness_/2;

            new G4PVPlacement(0, G4ThreeVector(sign_*x,y,+z), al_back_log, pName, sc_wall_log, false, i);
            new G4PVPlacement(0, G4ThreeVector(sign_*x,y,-z), al_back_log, pName, sc_wall_log, false, i);
            new G4PVPlacement(0, G4ThreeVector(sign_*x,yrc,0), rc_back_log, pName, sc_wall_log, false, i);
        } else {
            new G4PVPlacement(0, G4ThreeVector(sign_*x,y,0), al_back_log, pName, sc_wall_log, false, i);
        }
    }
    G4Colour backing_al_col;
    G4Colour::GetColour("aluminum", backing_al_col);
    al_back_log->SetVisAttributes(new G4VisAttributes(backing_al_col));
    G4Colour backing_rc_col;
    G4Colour::GetColour("carbon", backing_rc_col);
    if (rc_back_log) rc_back_log->SetVisAttributes(new G4VisAttributes(backing_rc_col));
    
    //------------------------- Frame
    
    G4double dd = 40 * mm;
    if (wall_type_ == front) {
        const G4double L0 = 209 * cm;
        const G4double Lx = 117 * cm;
        const G4double Ly = 257 * cm;
        
        const G4double y0 = (Ly - 42 * cm - 34 * cm - 8 * cm) / 2. + 42 * cm + 4 * cm;
        const G4double y1 = y0 - 42.0 * cm - 80 * mm;
        const G4double y2 = y0 - Ly + 34.0 * cm + 80 * mm;
        const G4double y3 = y0 - Ly;
        const G4double x1 = Lx / 2 + dd / 2;
        const G4double xs = 6 * cm;
        
        place_extrusion(L0, X, xs, y0, dd/2, 80 * mm);
        place_extrusion(Lx, X, xs, y1, dd/2, 80 * mm);
        place_extrusion(Lx, X, xs, y2, dd/2, 80 * mm);
        place_extrusion(Lx, X, xs, y3, dd/2, 80 * mm);
        place_extrusion(Ly, Y, xs+x1, 0, dd/2);
        place_extrusion(Ly, Y, xs-x1, 0, dd/2);
        
        const G4double dy = 5 * cm;
        place_extrusion(L0, X, xs, y0, dd/2 - dy, 80 * mm);
        place_extrusion(Lx, X, xs, y3, dd/2 - dy, 80 * mm);
        place_extrusion(Ly, Y, xs+x1, 0, dd/2 - dy);
        place_extrusion(Ly, Y, xs-x1, 0, dd/2 - dy);
    } else {
        const G4double L0 = 209 * cm;
        const G4double Lx = 172 * cm;
        const G4double Ly = 257 * cm + 8 * cm;
        
        const G4double y0 = 257 * cm / 2 + 8 * cm / 2;
        const G4double y1 = y0 - Ly / 3.;
        const G4double y2 = y0 - Ly / 3. * 2.;
        const G4double y3 = y0 - Ly;
        const G4double x1 = Lx / 2 + dd / 2;
        
        place_extrusion(L0, X, 0, y0, dd/2, 80 * mm);
        place_extrusion(Lx, X, 0, y1, dd/2, 80 * mm);
        place_extrusion(Lx, X, 0, y2, dd/2, 80 * mm);
        place_extrusion(Lx, X, 0, y3, dd/2, 80 * mm);
        place_extrusion(Ly, Y,  x1, -4 * cm, dd/2);
        place_extrusion(Ly, Y, -x1, -4 * cm, dd/2);
        
        const G4double dy = 11 * cm;
        place_extrusion(L0, X, 0, y0, dd/2 - dy, 80 * mm);
        place_extrusion(Lx, X, 0, y3, dd/2 - dy, 80 * mm);
        place_extrusion(Ly, Y,  x1, -4 * cm, dd/2 - dy);
        place_extrusion(Ly, Y, -x1, -4 * cm, dd/2 - dy);
    }
}


void g4PSISPS::SetSD(G4SDManager *SDman) {
    if (sc_Nunits_ > 0) {
        G4String SDname = "g4PSI/SPS/" + label_;
        G4VSensitiveDetector* SD = SDman->FindSensitiveDetector(SDname);
        if (SD == NULL) {
            SD = SD_ = new g4PSIScintillatorSD(SDname, 2 * sc_Nunits_, label_ + "_Collection" );
            SDman->AddNewDetector( SD );
        };
        sc_log_->SetSensitiveDetector( SD );
    }
}


void g4PSISPS::Write() {
    
    /// \todo cleanup; many parameters not in use; information about left/right missing
    
    TVectorD v1(6);
    v1[0] = 2*sc_Nunits_;
    v1[1] = 0;
    v1[2] = sc_distance_to_front_;
    v1[3] = sc_width_;
    v1[4] = sc_thickness_;
    v1[5] = sc_bar_length_;
    v1.Write(label_);
}


void g4PSISPS::InitTree(TTree *T) {
    SD_->InitTree(T);
}


void g4PSISPS::DeleteEventData() {
    SD_->DeleteEventData();
}


void g4PSISPS::place_extrusion(G4double L, orientation o, G4double x, G4double y, G4double z, G4double dy) {
    
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    G4Box* ex = NULL;
    G4double dd = 40 * mm;
    if (o == X) {
        ex = new G4Box((label_+"_ex").c_str(), L/2, dd/2, dy/2);
    } else if (o == Z) {
        ex = new G4Box((label_+"_ex").c_str(), dd/2, L/2, dd/2);
    } else if (o == Y) {
        ex = new G4Box((label_+"_ex").c_str(), dd/2, dd/2, L/2);
    }
    
    G4LogicalVolume *ex_log = new G4LogicalVolume(ex, Al, (label_+"_ex_log").c_str(), 0, 0, 0);

    G4double x0 = (sps_r_rear_of_backing_ + z) * sin(sc_theta_center_) + sign_ * (upstream_shift_ + x) * cos(sc_theta_center_);
    G4double z0 = (sps_r_rear_of_backing_ + z) * cos(sc_theta_center_) - sign_ * (upstream_shift_ + x) * sin(sc_theta_center_);
    
    new G4PVPlacement(sps_rot_, G4ThreeVector(x0, y, z0), ex_log, (label_+"_ex_phys").c_str(), mother_volume_, false, 0);
    
    G4Colour ex_col;
    G4Colour::GetColour("aluminum", ex_col);
    G4VisAttributes* VisAtt = new G4VisAttributes(ex_col);
    //VisAtt->SetForceWireframe(true);
    ex_log->SetVisAttributes(VisAtt);
    //    sc_wall_log->SetVisAttributes (G4VisAttributes::Invisible);
}
