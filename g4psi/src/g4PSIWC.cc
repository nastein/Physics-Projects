#include "g4PSIWC.hh"
#include "g4PSITrackerSD.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"   
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSIWC::g4PSIWC(G4String label, G4double angle, G4double r, G4double x, G4double y)
: g4PSIDetectorBase(label) {
  wc_label_ = label;
  wc_active_x_ = x;
  wc_active_y_ = y;
  wc_angle_ = angle;
  wc_r_ = r;
  wc_z_ =  9.06*cm;
}

g4PSIWC::~g4PSIWC() {
}


void g4PSIWC::Info() {
  G4cout << "Beam WC Detector:" << G4endl;
  G4cout << "==================" << G4endl;
  G4cout << G4endl;
  G4cout << "                                     name = " << wc_label_ << G4endl;
  G4cout << "                                    Angle = " << wc_angle_/deg << " deg" << G4endl;
  G4cout << "          Position of chamber center at r = " << wc_r_/cm << " cm" << G4endl;
  G4cout << "                                 Active x = " << wc_active_x_/cm << " cm" << G4endl;
  G4cout << "                                 Active y = " << wc_active_y_/cm << " cm" << G4endl;
  G4cout << G4endl;
  
}

G4LogicalVolume* g4PSIWC::wc_layer(G4Material *material, 
				   G4double dz, G4String name,
				   G4Colour col) {
  G4VSolid* solid = new G4Box((name+"_box").c_str(), wc_active_x_/2., wc_active_y_/2., dz/2.);
  G4LogicalVolume* lv = new G4LogicalVolume(solid, material, (name+"_box").c_str(), 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0, wc_running_z_ + 0.5*dz),
		    lv, (name+"_box").c_str(), wc_assembly_log_, false, 0);
  lv->SetVisAttributes(new G4VisAttributes(col));
  wc_running_z_ += dz;
  return lv;
}

void g4PSIWC::Placement() {

  // Material
  G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");  
  G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");  
//  G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");  
  G4Material* Argon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");  
//  G4Material* Copper = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");  
//  G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");  
  G4Material* CO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

  // WCGas at 2 atm & room temprature, 90% Argon and 10% CO2 gas

  G4Material* WCGas = G4Material::GetMaterial("WCGas", false);
  if(!WCGas) {
      G4double density = 2 * 1.678*mg/cm3;
      G4String name;
      G4int ncomponents;
      G4double fractionmass;
      WCGas = new G4Material(name="WCGas", density, ncomponents=2);
      WCGas->AddMaterial(Argon, fractionmass=0.9);
      WCGas->AddMaterial(CO2, fractionmass=0.1);
  } 

  G4Colour wc_colour;
  G4Colour gas_colour;
  G4Colour::GetColour("wc_detector", wc_colour);
  G4Colour::GetColour("wc_detector_gas", gas_colour);
  
  //------------------------- Detector
  G4RotationMatrix* zRot = new G4RotationMatrix;
  zRot->rotateY(-wc_angle_);
  G4Box* wc_detector = new G4Box("wc_detector", wc_active_x_/2, wc_active_y_/2, wc_z_/2);
  
  wc_assembly_log_ = new G4LogicalVolume(wc_detector, Air, (wc_label_+"_log").c_str(), 0, 0, 0);
  new G4PVPlacement(zRot, G4ThreeVector(wc_r_ * sin(wc_angle_), 0, wc_r_ * cos(wc_angle_)),
		    wc_assembly_log_, (wc_label_+"_phys").c_str(), mother_volume_, false, 0,
		    false);
  G4VisAttributes* wc_assembly_logVisAtt = new G4VisAttributes(wc_colour);
  wc_assembly_logVisAtt->SetForceWireframe(true);
  wc_assembly_log_->SetVisAttributes(wc_assembly_logVisAtt);
    
  //------------------------- The Detector Parts
  //
  //  the dimensions and materials: 
  //  * split the gas volume in two and attach to the second part the SD.
  //


  // Simplified straw-chamber setup
  //
  const int no_of_planes = 10;
  const double straw_spacing = 1.01*cm;
  const double straw_wall_thickness = 30.*um;
  
  const double t_mylar = pi * straw_wall_thickness * no_of_planes;    // <t> = (2 * pi * r * t) / (2 * r)
  const double t_arco2 = pi/2 * (straw_spacing/2) * no_of_planes;     // <t> = pi/2 * r

  wc_running_z_ = -wc_z_/2.;
  
  wc_layer(Mylar,   t_mylar/2.,  (wc_label_+"_mylarus1").c_str(), wc_colour);
  wc_layer(WCGas,   t_arco2/2.,  (wc_label_+"_gas1").c_str(), gas_colour);
  wc_detector_log_ = 
    wc_layer(WCGas, t_arco2/2.,  (wc_label_+"_gas2").c_str(), gas_colour);
  wc_layer(Mylar,   t_mylar/2.,  (wc_label_+"_mylards2").c_str(), wc_colour);
}


void g4PSIWC::SetSD(G4SDManager *SDman) {
  G4String WCSDname = "g4PSI/WC/" + wc_label_;
  G4VSensitiveDetector* WCSD = SDman->FindSensitiveDetector(WCSDname);
  if (WCSD == NULL) {
    WCSD = SD_ = new g4PSITrackerSD(WCSDname, wc_label_ + "_Collection" );
    SDman->AddNewDetector( WCSD );
  };
  wc_detector_log_->SetSensitiveDetector( WCSD ); 
}


bool g4PSIWC::hitsWall(double vx, double vy, double vz, 
		       double x0, double y0, double z0) {

  double st = sin(wc_angle_);
  double ct = cos(wc_angle_);

  // lam1: direction scale in the horizontal direction (x,z)
  // lam2: direction scale in the vertical direction (y)

//  double r = wc_r_ - wc_z_ / 2.0;
  double r = wc_r_;
  double h = vz * ct + vx * st;

  double mu = -(-r + z0*ct + x0*st) / h;
  double lam1 = (-vz*x0 + vx*z0 - r*vx*ct + r*vz*st) / h;
  double lam2 = (r*vy + (vz*y0 - vy*z0)*ct + (-vy*x0 + vx*y0)*st) / h;

  return fabs(lam1) < wc_active_x_/2.0 && fabs(lam2) < wc_active_y_/2.0 && mu > 0;
}


void g4PSIWC::Write() {
  TVectorD v1(4);
  v1[0] = wc_r_;
  v1[1] = wc_angle_;
  v1[2] = wc_active_x_;
  v1[3] = wc_active_y_;
  v1.Write(wc_label_);
}

void g4PSIWC::InitTree(TTree *T) {
  SD_->InitTree(T);
}

void g4PSIWC::DeleteEventData() {
  SD_->DeleteEventData();
}

