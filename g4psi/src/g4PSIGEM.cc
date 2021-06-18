#include "g4PSIGEM.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSIGEM::g4PSIGEM(G4String label, G4double dGEM, G4double dFRAME, G4double z)
    : g4PSIDetectorBase(label)
{
    gem_label_ = label;
    gem_dx_ = dGEM;
    gem_dy_ = dGEM;
    gem_dx_frame_ = dFRAME;
    gem_dy_frame_ = dFRAME;
    gem_r_ = z;
    gem_angle_ = 0.;
}


g4PSIGEM::g4PSIGEM(G4String label, G4double angle, G4double r, G4double x, G4double y, G4double xout, G4double yout) : g4PSIDetectorBase(label)
{
    gem_label_ = label;
    gem_dx_ = x;
    gem_dy_ = y;
    gem_dx_frame_ = xout;
    gem_dy_frame_ = yout;
    gem_r_ = r;
    gem_angle_ = angle;
}


g4PSIGEM::~g4PSIGEM()
{

}

double g4PSIGEM::GetZPos() {return gem_r_;}

G4LogicalVolume* g4PSIGEM::gem_layer(G4Material *material,
				     G4double thickness_weight,
                                     G4double dz, G4String name,
				     G4LogicalVolume *mv,
                                     G4Colour col)
{
    dz *= thickness_weight;
    G4VSolid* solid = new G4Box((name+"_box").c_str(), gem_dx_/2., gem_dy_/2., dz/2.);
    G4LogicalVolume* lv = new G4LogicalVolume(solid, material, (name+"_box").c_str(), 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, gem_running_z_ + 0.5*dz),
                      lv, (name+"_box").c_str(), mv, false, 0);
    lv->SetVisAttributes(new G4VisAttributes(col));
    gem_running_z_ += dz;
    return lv;
}


void g4PSIGEM::Info()
{
    InfoTitle("GEM Detector");
    InfoParDouble("Position z", gem_r_/cm, "cm");
    InfoPar2Double("Size", gem_dx_/cm, gem_dy_/cm, "cm");
    
//    G4cout << "Beam GEM Detector:" << G4endl;
//    G4cout << "==================" << G4endl;
//    G4cout << G4endl;
//    G4cout << "                                     name = " << gem_label_ << G4endl;
//    G4cout << "                            Position at z = " << gem_r_/cm << " cm" << G4endl;
//    G4cout << "                                     Size = " << gem_dx_/cm << " cm x "  << gem_dy_/cm << " cm" << G4endl;
//    G4cout << G4endl;

}


void g4PSIGEM::Placement()
{
    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Copper = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
    G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");

    G4String name;

    // GEMGas at 1 atm & room temprature, 70% Argon and 30% CO2 gas
    name = "GEMGas";
    G4Material *GEMGas = G4Material::GetMaterial(name, false);
    if (!GEMGas) {
	G4double density = 1.716*mg/cm3;
	G4int ncomponents = 2;
	G4double fractionmass;
	G4Material* Argon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
	G4Material* CO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
	GEMGas = new G4Material(name, density, ncomponents);
	GEMGas->AddMaterial(Argon, fractionmass=0.7);
	GEMGas->AddMaterial(CO2, fractionmass=0.3);
    }

    name = "NemaG10";
    G4Material *G10 = G4Material::GetMaterial(name, false);
    if (!G10) {
	G4double density = 1.700*g/cm3;
	G4int ncomponents = 4;
	G4int natoms;
	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4Element* elH = new G4Element("Hydrogen", "H", z=1., a=1.00794*g/mole);
	G4Element* elC = new G4Element("Carbon", "C", z=6., a=12.01*g/mole);
	G4Element* elO  = new G4Element("Oxygen"  , "O", z= 8., a = 16.00*g/mole);
	G4Element* elSi = new G4Element("Silicon", "Si", z= 14., a=28.09*g/mole);
	G10 = new G4Material("NemaG10", density, ncomponents);
	G10->AddElement(elSi, natoms=1);
	G10->AddElement(elO , natoms=2);
	G10->AddElement(elC , natoms=3);
	G10->AddElement(elH , natoms=3);
    }

    G4Colour gem_colour;
    G4Colour gem_frame_colour;
    G4Colour gas_colour;
    G4Colour::GetColour("gem_detector", gem_colour);
    G4Colour::GetColour("gem_detector_frame", gem_frame_colour);
    G4Colour::GetColour("gem_detector_gas", gas_colour);

    //------------------------- Detector
    G4double gem_z =  2.*cm;
    G4Box* gem_frame = new G4Box((gem_label_+"_box").c_str(), gem_dx_frame_/2, gem_dy_frame_/2, gem_z/2);
    gem_assembly_log_ = new G4LogicalVolume(gem_frame, G10, (gem_label_+"_log").c_str(), 0, 0, 0);
    G4Box* gem_detector = new G4Box((gem_label_+"_air").c_str(), gem_dx_/2, gem_dy_/2, gem_z/2);
    G4LogicalVolume *gem_inner_volume_log = new G4LogicalVolume(gem_detector, Air, (gem_label_+"_log").c_str(), 0, 0, 0);

    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateY(-gem_angle_);
    new G4PVPlacement(zRot, G4ThreeVector(gem_r_ * sin(gem_angle_), 0, gem_r_ * cos(gem_angle_)),
                      gem_assembly_log_, (gem_label_+"_frame_phys").c_str(), mother_volume_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                      gem_inner_volume_log, (gem_label_+"_detector_phys").c_str(), gem_assembly_log_, false, 0);

    G4VisAttributes* VisAtt = new G4VisAttributes(gem_frame_colour);
    //  VisAtt->SetForceWireframe(true);
    gem_assembly_log_->SetVisAttributes(VisAtt);

    //------------------------- The Detector Parts
    //
    //  the dimensions and materials: info from M. Kohl
    //
    gem_running_z_ = -gem_z/2.;
 
    // GAS SEAL
    gem_layer(Mylar,    1.0, 0.020*mm, (gem_label_+"_mylar0").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Al,       1.0, 0.001*mm, (gem_label_+"_al0").c_str(), gem_inner_volume_log, gem_colour);

    gem_layer(GEMGas,   1.0, 3.000*mm, (gem_label_+"_gas1").c_str(), gem_inner_volume_log, gas_colour);

    // HV (copper on one side only)
    gem_layer(Kapton,   0.7, 0.050*mm, (gem_label_+"_kapton1").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu1").c_str(), gem_inner_volume_log, gem_colour);
    gem_detector_log_ =
        gem_layer(GEMGas,   1.0, 3.000*mm, (gem_label_+"_gas2").c_str(), gem_inner_volume_log, gas_colour);

    // 1st-GEM
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu2a").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Kapton,   0.7, 0.050*mm, (gem_label_+"_kapton2").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu2b").c_str(), gem_inner_volume_log, gem_colour);

    gem_layer(GEMGas,   1.0, 2.000*mm, (gem_label_+"_gas3").c_str(), gem_inner_volume_log, gas_colour);

    // 2nd-GEM
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu3a").c_str(), gem_inner_volume_log, gem_colour);
	gem_layer(Kapton,   0.7, 0.050*mm, (gem_label_+"_kapton3").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu3b").c_str(), gem_inner_volume_log, gem_colour);

    gem_layer(GEMGas,   1.0, 2.000*mm, (gem_label_+"_gas4").c_str(), gem_inner_volume_log, gas_colour);

    // 3rd-GEM
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu4a").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Kapton,   0.7, 0.050*mm, (gem_label_+"_kapton4").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu4b").c_str(), gem_inner_volume_log, gem_colour);

    gem_layer(GEMGas,   1.0, 2.000*mm, (gem_label_+"_gas5").c_str(), gem_inner_volume_log, gas_colour);

    // Readout
    gem_layer(Copper,   0.7, 0.005*mm, (gem_label_+"_cu5a").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Kapton,   0.99, 0.050*mm, (gem_label_+"_kapton5").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.01, 0.05*mm, (gem_label_+"_cu5k").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Copper,   0.3, 0.005*mm, (gem_label_+"_cu5b").c_str(), gem_inner_volume_log, gem_colour);

    gem_layer(GEMGas,   1.0, 3.000*mm, (gem_label_+"_gas5").c_str(), gem_inner_volume_log, gas_colour);

    // GAS SEAL
    gem_layer(Al,       1.0, 0.001*mm, (gem_label_+"_al6").c_str(), gem_inner_volume_log, gem_colour);
    gem_layer(Mylar,    1.0, 0.020*mm, (gem_label_+"_mylar6").c_str(), gem_inner_volume_log, gem_colour);
}

void g4PSIGEM::SetSD(G4SDManager *SDman)
{
    G4String GEMSDname = "g4PSI/GEM/" + gem_label_;
    G4VSensitiveDetector* GEMSD = SDman->FindSensitiveDetector(GEMSDname);
    if (GEMSD == NULL) {
        GEMSD = SD_ = new g4PSITrackerSD(GEMSDname, gem_label_ + "_Collection" );
        SDman->AddNewDetector( GEMSD );
    };
    gem_detector_log_->SetSensitiveDetector( GEMSD );
}


void g4PSIGEM::Write()
{
    TVectorD v1(1);
    v1[0] = gem_r_;
    v1.Write(gem_label_);
}
void g4PSIGEM::InitTree(TTree *T)
{
    SD_->InitTree(T);
}

void g4PSIGEM::DeleteEventData()
{
    SD_->DeleteEventData();
}
