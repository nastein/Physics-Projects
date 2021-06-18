#include "g4PSISCWall.hh"
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

g4PSISCWall::g4PSISCWall(G4String label,
                         G4double width, G4double thickness, G4double length,
                         G4double posX, G4double posZ, G4RotationMatrix *rot, G4Material* mat) : g4PSIDetectorBase(label) {
    
    sc_N_up_ = 0;
    sc_N_down_ = 0;
    sc_distance_to_front_ = 0;
    sc_thetaXmin_ = 0;
    sc_thetaXmax_ = 0;
    sc_thetaYmin_ = 0;
    sc_thetaYmax_ = 0;
    
    sc_width_ = width;
    sc_thickness_ = thickness;
    sc_bar_length_ = length;
    sc_posX_ = posX;
    sc_posZ_ = posZ;
    sc_mode_ = 0;
    sc_al_thickness_ = 0;
    sc_rot_ = rot;
    sc_material_ = mat;
}


g4PSISCWall::g4PSISCWall(G4String label,
                         G4double width, G4double gap, G4double thickness, G4double height,
                         G4double distance, G4double theta, G4int Nup, G4int Ndown,
                         G4double al_thickness) : g4PSIDetectorBase(label) {
    
    sc_width_ = width;
    sc_thickness_ = thickness;
    sc_gap_ = gap;
    sc_distance_to_front_ = distance;
    sc_bar_length_ = height;
    sc_thetaXmin_ = 0;
    sc_thetaXmax_ = 0;
    sc_thetaYmin_ = 0;
    sc_thetaYmax_ = 0;
    sc_al_thickness_ = al_thickness;
    sc_rot_ = NULL;
    
    sc_theta_center_ = theta;
    
    sc_N_up_ = Nup;
    sc_N_down_ = Ndown;
    
    sc_dx_up_ = sc_N_up_ * (width + gap) + width / 2. + gap;
    sc_dx_down_ = sc_N_down_ * (width + gap) + width / 2. + gap;
    sc_posX_ = 0.;
    sc_posZ_ = 0.;
    sc_mode_ = 1;
    sc_material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
}


g4PSISCWall::g4PSISCWall(G4String label,
                         G4double width, G4double gap, G4double thickness,
                         G4double distance, G4double target_z, G4double target_r,
                         G4double safety_u1, G4double safety_d1, G4double safety_y1,
                         G4double thetaXmin, G4double thetaXmax,
                         G4double thetaYmin, G4double thetaYmax): g4PSIDetectorBase(label) {
    
    // let's protect us from input mistakes about min and max angles.
    
    if (thetaXmin > thetaXmax) {
        G4double tmp = thetaXmin;
        thetaXmin = thetaXmax;
        thetaXmax = tmp;
    }
    if (thetaYmin > thetaYmax) {
        G4double tmp = thetaYmin;
        thetaYmin = thetaYmax;
        thetaYmax = tmp;
    }
    
    G4double theta1 = fabs(thetaXmax) > fabs(thetaXmin) ? fabs(thetaXmin) : fabs(thetaXmax);
    G4double theta2 = fabs(thetaXmax) > fabs(thetaXmin) ? fabs(thetaXmax) : fabs(thetaXmin);
    G4double thetaC = (theta1 + theta2) / 2.;
    
    G4double dx_tgt = target_z/2. * (sin(theta1)/sin(theta1+thetaC) + sin(theta2)/sin(90*deg+thetaC-theta2));
    G4double dy_tgt = target_r * 2;
    
    G4double dxu = (distance + thickness)*tan(fabs(theta2-theta1)/2.+safety_u1) + dx_tgt;  // downstream width
    G4double dxd = (distance + thickness)*tan(fabs(theta2-theta1)/2.+safety_d1) + dx_tgt;  // upstream width
    G4double dy = 2.*(distance + thickness)*tan(fabs(thetaYmax-thetaYmin)/2.+safety_y1) + dy_tgt;
    
    ///
    
    sc_width_ = width;
    sc_thickness_ = thickness;
    sc_gap_ = gap;
    sc_distance_to_front_ = distance;
    sc_bar_length_ = dy;
    sc_thetaXmin_ = thetaXmin;
    sc_thetaXmax_ = thetaXmax;
    sc_thetaYmin_ = thetaYmin;
    sc_thetaYmax_ = thetaYmax;
    sc_al_thickness_ = 0;
    
    sc_rot_ = NULL;
    
    sc_theta_center_ = (sc_thetaXmin_ + sc_thetaXmax_) / 2.0;
    
    sc_N_up_ = ceil((dxu - gap - width / 2.) / (width + gap));
    sc_N_down_ = ceil((dxd - gap - width / 2.) / (width + gap));
    
    sc_dx_up_ = sc_N_up_ * (width + gap) + width / 2. + gap;
    sc_dx_down_ = sc_N_down_ * (width + gap) + width / 2. + gap;
    sc_posX_ = 0.;
    sc_posZ_ = 0.;
    sc_mode_ = 1;
    
    sc_material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
}

g4PSISCWall::~g4PSISCWall() {
}


double g4PSISCWall::GetZPos() {return sc_posZ_;}

void g4PSISCWall::Info() {
    if (sc_mode_ == 1) {
        InfoTitle("Scintillator Wall");
        InfoParString("Material", sc_material_->GetName());
        InfoParInt("Number of scintillator bars Nup+Ndown+1", sc_N_up_ + sc_N_down_ + 1);
        InfoParInt("Nup", sc_N_up_);
        InfoParInt("Ndown", sc_N_down_);
        InfoPar3Double("Physical length", sc_width_/cm, sc_thickness_/cm, sc_bar_length_/cm, "cm");
        InfoParDouble("Front-face distance to the target center", sc_distance_to_front_/cm, "cm");
        InfoParDouble("Gapwidth between scintillator bars", sc_gap_/cm, "cm");
    } else {
        InfoTitle("Scintillator");
        InfoParString("Material", sc_material_->GetName());
        InfoPar3Double("Physical length", sc_width_/cm, sc_thickness_/cm, sc_bar_length_/cm, "cm");
    }
}


void g4PSISCWall::Placement() {
    
    // Material
    /// \todo don't use Air, use the material of the world volume
    
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    
    if (sc_mode_ == 1) {
        
        // Place a full wall of scintillator bars  (Mode == 1)
        // ---------------------------------------------------
        
        G4RotationMatrix* zRot = new G4RotationMatrix;
        zRot->rotateX(90*deg);
        zRot->rotateZ(-sc_theta_center_);
        
        //------------------------- Wall
        //    x0 and z0 are the coordinates of the center of the wall
        //    construct and place container
        
        G4double sign = (sc_theta_center_ > 0) ? 1.0 : -1.0;
        G4double xx = (sc_dx_up_ - sc_dx_down_) / 2.0;
        G4double x0 = (sc_distance_to_front_ + sc_thickness_ / 2.0) * sin(sc_theta_center_) + sign * xx * cos(sc_theta_center_);
        G4double z0 = (sc_distance_to_front_ + sc_thickness_ / 2.0) * cos(sc_theta_center_) - sign * xx * sin(sc_theta_center_);
        
        G4Box* sc_wall = new G4Box((label_+"sc_wall").c_str(), (sc_dx_up_ + sc_dx_down_)/2, sc_thickness_/2, sc_bar_length_/2);
        G4LogicalVolume *sc_wall_log = new G4LogicalVolume(sc_wall, Air, (label_+"_m1_log").c_str(), 0, 0, 0);
        new G4PVPlacement(zRot, G4ThreeVector(x0, 0, z0),
                          sc_wall_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
        if (sc_al_thickness_ > 0) {
            G4Box* sc_back = new G4Box((label_+"sc_back").c_str(), (sc_dx_up_ + sc_dx_down_)/2, sc_al_thickness_/2, sc_bar_length_/2 + 0.0*cm);
            G4LogicalVolume *sc_back_log = new G4LogicalVolume(sc_back, Air, (label_+"_back_log").c_str(), 0, 0, 0);
            G4double r = sc_distance_to_front_ + sc_thickness_ + sc_al_thickness_ / 2.0;
            G4double x1 = r * sin(sc_theta_center_) + sign * xx * cos(sc_theta_center_);
            G4double z1 = r * cos(sc_theta_center_) - sign * xx * sin(sc_theta_center_);
            new G4PVPlacement(zRot, G4ThreeVector(x1, 0, z1),
                              sc_back_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
            G4Colour sc_col;
            G4Colour::GetColour("aluminum", sc_col);
            sc_back_log->SetVisAttributes(new G4VisAttributes(sc_col));
        }
        
        // G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(1,0,0));
        // VisAtt->SetForceWireframe(true);
        // sc_wall_log->SetVisAttributes(VisAtt);
        sc_wall_log->SetVisAttributes (G4VisAttributes::Invisible);
        //------------------------- SC Bars
        
        G4double dx = sc_width_ + sc_gap_;
        G4Box* sc_bar = new G4Box("sc", sc_width_/2, sc_thickness_/2, sc_bar_length_/2);
        sc_log_ = new G4LogicalVolume(sc_bar, sc_material_, (label_+"_bar_log").c_str(), 0, 0, 0);
        
        // G4VPVParameterisation * SCParam = new g4PSISCWallParametrisation(sc_N_, dx);
        // new G4PVParameterised(label_.c_str(), sc_log_, sc_wall_log,
        // 			kXAxis, 2*sc_N_+1, SCParam);
        
        G4LogicalVolume *pmt_log = NULL;
        
        double pmt_height = 14*cm;
        
        if (g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kPMTs]) {
            double pmt_width = 5.08*cm;
            G4Colour pmt_col;
            G4Colour::GetColour("pmt", pmt_col);
            G4Tubs* pmt = new G4Tubs((label_ + "_pmt_geom").c_str(), 0, pmt_width/2, pmt_height/2, 0, 360*deg);
            pmt_log = new G4LogicalVolume(pmt, sc_material_, (label_ + "_pmt_log").c_str(), 0,0,0);
            pmt_log->SetVisAttributes(new G4VisAttributes(pmt_col));
        }
        
        for (G4int i = -sc_N_down_; i <= sc_N_up_; i++) {
            //-------------------------- SC bars
            
            char pName[60];
            sprintf(pName, "%s_%02d", label_.c_str(), sc_N_down_ + i + 1);
            new G4PVPlacement(0, G4ThreeVector(sign*(i*dx - xx),0,0), sc_log_, pName, sc_wall_log, false, sc_N_down_ + i);
            char pmtName[60];
            
            //-------------------------- PMTs
            
            if (pmt_log) {
                sprintf(pmtName, "%s_%02d_PMT1", label_.c_str(), sc_N_down_ + i + 1);
                new G4PVPlacement(0, G4ThreeVector(sign*(i*dx - xx), 0, (sc_bar_length_+pmt_height)/2.),
                                  pmt_log, pmtName, sc_wall_log, false, sc_N_down_ + i);
                sprintf(pmtName, "%s_%02d_PMT2", label_.c_str(), sc_N_down_ + i + 1);
                new G4PVPlacement(0, G4ThreeVector(sign*(i*dx - xx), 0, -(sc_bar_length_+pmt_height)/2.),
                                  pmt_log, pmtName, sc_wall_log, false, sc_N_down_ + i);
            }
        }
        
        
    } else {
        
        // Place a single bar  (Mode == 0)
        // -------------------------------
        
        G4Box* sc_bar = new G4Box("sc", sc_width_/2, sc_thickness_/2, sc_bar_length_/2);
        sc_log_ = new G4LogicalVolume(sc_bar, sc_material_, (label_+"_bar_log").c_str(), 0, 0, 0);
        new G4PVPlacement(sc_rot_, G4ThreeVector(sc_posX_, 0, sc_posZ_), sc_log_, (label_ + "_bar_phys").c_str(),
                          mother_volume_, false, 0);
        sc_log_->SetUserLimits(new G4UserLimits(0.5*mm));
    }
    G4Colour sc_col;
    G4Colour::GetColour("scwall", sc_col);
    sc_log_->SetVisAttributes(new G4VisAttributes(sc_col));
}


void g4PSISCWall::SetSD(G4SDManager *SDman) {
    if (sc_mode_ == 1) {
        if (sc_N_up_ + sc_N_down_ >= 0) {
            G4String SDname = "g4PSI/SCWall/" + label_;
            G4VSensitiveDetector* SD = SDman->FindSensitiveDetector(SDname);
            if (SD == NULL) {
                SD = SD_ = new g4PSIScintillatorSD(SDname, 1+sc_N_up_+sc_N_down_, label_ + "_Collection" );
                SDman->AddNewDetector( SD );
            };
            sc_log_->SetSensitiveDetector( SD );
        }
    } else {
        G4String SDname = "g4PSI/SCBar/" + label_;
        G4VSensitiveDetector* SD = SDman->FindSensitiveDetector(SDname);
        if (SD == NULL) {
            SD = SD_ = new g4PSIScintillatorSD(SDname, 1, label_ + "_Collection" );
            SDman->AddNewDetector( SD );
        };
        sc_log_->SetSensitiveDetector( SD );
    }
}


void g4PSISCWall::Write() {
    
    /// \todo cleanup; many parameters not in use; information about left/right missing
    
    TVectorD v1(15);
    v1[0] = sc_N_up_;
    v1[1] = sc_N_down_;
    v1[2] = sc_distance_to_front_;
    v1[3] = sc_width_;
    v1[4] = sc_thickness_;
    v1[5] = sc_bar_length_;
    v1[7] = sc_thetaXmin_;
    v1[8] = sc_thetaXmax_;
    v1[9] = sc_thetaYmin_;
    v1[10] = sc_thetaYmax_;
    v1[11] = sc_posX_;
    v1[12] = sc_posZ_;
    v1[13] = sc_mode_;
    v1[14] = sc_al_thickness_;
    
    v1.Write(label_);
}


bool g4PSISCWall::hitsWall(double vx, double vy, double vz) {
    
    G4double thetaX = 90.0 * deg;
    G4double thetaY = 0.0 * deg;
    
    if (vz != 0) thetaX = atan2(vx, vz);
    if (vx != 0) thetaY = atan2(vy, vx * sin(sc_theta_center_));
    
    return (thetaX > sc_thetaXmin_ && thetaX < sc_thetaXmax_ && thetaY > sc_thetaYmin_ && thetaY < sc_thetaYmax_);
}

// bool g4PSISCWall::hitsWall(double vx, double vy, double vz,
// 			   double x0, double y0, double z0,
// 			   double &x, double &y, double &path) {

//   double st = sin(sc_theta_center_);
//   double ct = cos(sc_theta_center_);

//   // lam1: direction scale in the horizontal direction (x,z)
//   // lam2: direction scale in the vertical direction (y)

//   double lam1 = 2.*(-vz*x0 + vx*z0 - sc_distance_to_front_*vx*ct + sc_distance_to_front_*vz*st) / (vz*ct + vx*st) / sc_active_length_;
//   double lam2 = 2.*(sc_distance_to_front_*vy + vz*y0*ct - vy*z0*ct - vy*x0*st + vx*y0*st) / (vz*ct + vx*st) / sc_active_length_;
//   double mu = -(-sc_distance_to_front_ + z0*ct + x0*st) / (vz*ct + vx*st);
//   x = lam1 * sc_active_length_ / 2.;
//   y = lam2 * sc_active_length_ / 2.;
//   path = mu / sqrt(vx*vx + vy*vy + vz*vz);
//   return fabs(lam1) < 1.0 && fabs(lam2) < 1.0 && mu > 0;
// }


// bool g4PSISCWall::hitsAnyWall(double vx, double vy, double vz,
// 			      double x0, double y0, double z0) {

//   double st = sin(theta0_);
//   double ct = cos(theta0_);

//   // lam1: direction scale in the horizontal direction (x,z)
//   // lam2: direction scale in the vertical direction (y)

//   // check first wall:
//   double lam1 = 2.*(-vz*x0 + vx*z0 - rMin1_*vx*ct + rMin1_*vz*st) / (vz*ct + vx*st) / active_length1_;
//   double lam2 = 2.*(rMin1_*vy + vz*y0*ct - vy*z0*ct - vy*x0*st + vx*y0*st) / (vz*ct + vx*st) / active_length1_;
//   double mu = -(-rMin1_ + z0*ct + x0*st) / (vz*ct + vx*st);
//   if (fabs(lam1) < 1.0 && fabs(lam2) < 1.0 && mu > 0) return true;

//   // check second wall:
//   lam1 = 2.*(-vz*x0 + vx*z0 - rMin2_*vx*ct + rMin2_*vz*st) / (vz*ct + vx*st) / active_length2_;
//   lam2 = 2.*(rMin2_*vy + vz*y0*ct - vy*z0*ct - vy*x0*st + vx*y0*st) / (vz*ct + vx*st) / active_length2_;
//   mu = -(-rMin2_ + z0*ct + x0*st) / (vz*ct + vx*st);
//   return fabs(lam1) < 1.0 && fabs(lam2) < 1.0 && mu > 0;
// }

void g4PSISCWall::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSISCWall::DeleteEventData() {
    SD_->DeleteEventData();
}

