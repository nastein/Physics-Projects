#include "g4PSISTC.hh"
#include "g4PSITrackerSD.hh"
#include "g4PSIDetectorParts.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"

using namespace CLHEP;

g4PSISTC::g4PSISTC(G4String label, G4double angle, STTChamber part)
: g4PSIDetectorBase(label) {
    STT_Init(label, angle);
    stc_on_[part] = true;
}


g4PSISTC::g4PSISTC(G4String label, G4double angle)
: g4PSIDetectorBase(label) {
    STT_Init(label, angle);
    stc_on_[kFront] = true;
    stc_on_[kRear] = true;
}


void g4PSISTC::STT_Init(G4String label, G4double angle)
{
    stc_label_ = label;
    stc_angle_ = angle;
    
    // initialize geometry
    
    stc_r_[kFront] = 30.0*cm;
    stc_active_y_min_[kFront] = 0;
    stc_active_y_max_[kFront] = 55 * cm;
    stc_active_x_min_[kFront] = 0;
    stc_active_x_max_[kFront] = 60 * cm;
    NUMBER_OF_STRAWS_IN_YPLANE_[kFront] = 55;
    NUMBER_OF_STRAWS_IN_XPLANE_[kFront] = 60;
    straw_y_length_[kFront] = 60 * cm;
    straw_x_length_[kFront] = 55 * cm;
     
//    stc_active_y_min_[kFront] = 392.95 * mm;
//    stc_active_y_max_[kFront] = 1007.05 * mm;
//    stc_active_x_min_[kFront] = 39.5 * mm;
//    stc_active_x_max_[kFront] = 600.5 * mm;
    
    stc_r_[kRear] = 45.0*cm;
    stc_active_y_min_[kRear] = 0;
    stc_active_y_max_[kRear] = 90 * cm;
    stc_active_x_min_[kRear] = 0;
    stc_active_x_max_[kRear] = 90 * cm;
    NUMBER_OF_STRAWS_IN_YPLANE_[kRear] = 80;
    NUMBER_OF_STRAWS_IN_XPLANE_[kRear] = 90;
    straw_y_length_[kRear] = 90 * cm;
    straw_x_length_[kRear] = 80 * cm;

    
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        SD_[i] = NULL;
        SDF_[i] = NULL;
        stc_plane_log_[i] = NULL;
        stc_assembly_log_[i] = NULL;
        stc_frameX_log_[i] = NULL;
        stc_frameY_log_[i] = NULL;
        stc_active_y_[i] = stc_active_y_max_[i] - stc_active_y_min_[i];
        stc_active_x_[i] = stc_active_x_max_[i] - stc_active_x_min_[i];
        stc_tube_y_log_[NUMBER_OF_STT_CHAMBERS] = NULL;
        stc_tube_x_log_[NUMBER_OF_STT_CHAMBERS] = NULL;
    }
    
    t_mylar_ = pi * straw_wall_thickness_;  // <t> = (2 * pi * r * t) / (2 * r)
    t_arco2_ = pi/2 * (straw_spacing_/2);   // <t> = pi/2 * r
    stc_plane_z_ = t_mylar_ + t_arco2_;
    stc_full_chamber_width_ = (no_of_xplanes_ + no_of_yplanes_ - 2) * dz_planes_ + 2*dz_offset_ + stc_plane_z_;
}


g4PSISTC::~g4PSISTC() {
}


void g4PSISTC::Info() {
    InfoTitle("Straw Chamber");
    InfoParDouble("Angle", stc_angle_/deg, " deg");
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        InfoParInt("Chamber", i);
        InfoParDouble(Form("Distance from center r[%d]",i), stc_r_[i]/cm, " cm");
        InfoParDouble(Form("Active x[%d]", i), stc_active_x_[i]/cm, " cm");
        InfoParDouble(Form("Active y[%d]", i), stc_active_y_[i]/cm, " cm");
    }
    InfoParDouble("Effective plane thickness, mylar", t_mylar_/cm, " cm");
    InfoParDouble("Effective plane thickness, ArCo2", t_arco2_/cm, " cm");
    InfoParDouble("Effective plane thickness, total", stc_plane_z_/cm, " cm");
    InfoParDouble("Full chamber width", stc_full_chamber_width_/cm, " cm");
    InfoParDouble("Adjacent straw-plane offset", dz_planes_/cm, " cm");
    InfoParDouble("X-Y plane offset", dz_offset_/cm * 2, " cm");
    InfoParInt("No of x-planes", no_of_xplanes_);
    InfoParInt("No of y-planes", no_of_yplanes_);
    InfoParBool("Frame", !g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kWireChamber_noframe]);
}


void g4PSISTC::Placement() {
    
    // Material
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Argon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
    //  G4Material* Copper = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
    //  G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    G4Material* CO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    
    // STCGas at 2 atm & room temprature, 90% Argon and 10% CO2 gas
    // STCGas    density:  3.356 mg/cm3  RadL:  61.058 m
    G4Material* STCGas = G4Material::GetMaterial("STCGas", false);
    if(!STCGas) {
        G4double density = 2 * 1.678*mg/cm3;
        G4String name;
        G4int ncomponents;
        G4double fractionmass;
        STCGas = new G4Material(name="STCGas", density, ncomponents=2);
        STCGas->AddMaterial(Argon, fractionmass=0.9);
        STCGas->AddMaterial(CO2, fractionmass=0.1);
    }
    
    G4Colour stc_colour;
    G4Colour gas_colour;
    G4Colour rf_shield_colour;
    G4Colour::GetColour("stc_detector", stc_colour);
    G4Colour::GetColour("stc_rfshield", rf_shield_colour);
    G4Colour::GetColour("stc_detector_gas", gas_colour);
    
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateY(-stc_angle_);
    
    G4double stc_gap = stc_r_[kRear] - stc_r_[kFront] - stc_full_chamber_width_;
    
    // RF Shield
    // =========
    //
    // 2016      RF shield material is presently 0.085 mm copper.  This is an oversimplification
    //           and includes too much material.
    // 01/12/17  Polyester is almost identical to Mylar.  Use mylar as RF material for now.
    //           Use the manufacturers fabric weight information to determine the eff. thickness.
    //           Need to add nickle and copper to the material.
    // \todo create appropriate RF material
    
    G4Material* rfshield = Mylar;
    const G4double fabric_weight = 34 * g/m2;
    G4double t_shield = fabric_weight / rfshield->GetDensity(); // 0.085*mm;
    
    G4double rf_shield_front_dz = stc_full_chamber_width_ + stc_gap;
    
    G4Trd* rf_shield_front1 = new G4Trd("stc_rf_shield_front1", stc_active_x_[0]/2 + t_shield, stc_active_x_[1]/2 + t_shield, stc_active_y_[0]/2 + t_shield, stc_active_y_[1]/2 + t_shield, (rf_shield_front_dz + t_shield)/2);
    G4Trd* rf_shield_front2 = new G4Trd("stc_rf_shield_front2", stc_active_x_[0]/2, stc_active_x_[1]/2, stc_active_y_[0]/2, stc_active_y_[1]/2, rf_shield_front_dz/2);
    G4LogicalVolume* rf_shield_front1_log = new G4LogicalVolume(rf_shield_front1, rfshield, (stc_label_+ "_rf_shield_F1").c_str(), 0, 0, 0);
    G4LogicalVolume* rf_shield_front2_log = new G4LogicalVolume(rf_shield_front2, Air, (stc_label_+ "_rf_shield_F2").c_str(), 0, 0, 0);
    
    G4double r = stc_r_[kFront] + stc_gap / 2. - t_shield / 2.;

    new G4PVPlacement(zRot, G4ThreeVector(r * sin(stc_angle_), 0, r * cos(stc_angle_)),
                      rf_shield_front1_log, (stc_label_+"_F1_phys").c_str(), mother_volume_, false, 0, false);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, t_shield/2),
                      rf_shield_front2_log, (stc_label_+"_F2_phys").c_str(), rf_shield_front1_log, false, 0, false);
    
    G4Box* rf_shield_rear1 = new G4Box("stc_rf_shield_rear1", stc_active_x_[1]/2 + t_shield, stc_active_y_[1]/2 + t_shield, stc_full_chamber_width_/2 + t_shield/2);
    G4Box* rf_shield_rear2 = new G4Box("stc_rf_shield_rear2", stc_active_x_[1]/2, stc_active_y_[1]/2, stc_full_chamber_width_/2);
    G4LogicalVolume* rf_shield_rear1_log = new G4LogicalVolume(rf_shield_rear1, rfshield, (stc_label_ + "_rf_shield_R1").c_str(), 0, 0, 0);
    G4LogicalVolume* rf_shield_rear2_log = new G4LogicalVolume(rf_shield_rear2, Air, (stc_label_ + "_rf_shield_R2").c_str(), 0, 0, 0);
    
    r = stc_r_[kRear] + t_shield / 2.;

    new G4PVPlacement(zRot, G4ThreeVector(r * sin(stc_angle_), 0, r * cos(stc_angle_)),
                      rf_shield_rear1_log, (stc_label_+"_R1_phys").c_str(), mother_volume_, false, 0,
                      false);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, -t_shield/2),
                      rf_shield_rear2_log, (stc_label_+"_R2_phys").c_str(), rf_shield_rear1_log, false, 0, false);
    
    G4VisAttributes* stc_assembly_logVisAtt = new G4VisAttributes(rf_shield_colour);
    rf_shield_rear1_log->SetVisAttributes(stc_assembly_logVisAtt);
    rf_shield_front1_log->SetVisAttributes(stc_assembly_logVisAtt);
    rf_shield_rear2_log->SetVisAttributes(stc_assembly_logVisAtt);
    rf_shield_front2_log->SetVisAttributes(stc_assembly_logVisAtt);
    
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        
        //------------------------- Frame
        
        stc_frameX_log_[i] = NULL;
        stc_frameY_log_[i] = NULL;
        if (!g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kWireChamber_noframe]) {
            
            if (i == kFront) {
                G4double stc_frame_x = 1024.05*mm;  // outer dimension
                G4double stc_frame_y = 1470*mm;     // outer dimension
                G4double stc_frame_z = 100.0*mm;    // outer dimension
                G4double stc_frame_w =  20.0*mm;
                G4double sig = stc_angle_ > 0 ? 1.0 : -1.0;
                G4double th = fabs(stc_angle_);
                G4double stc_x0 = 56.51*mm + stc_frame_z/2 * sin(th) + stc_frame_w/2 * cos(th);
                G4double stc_z0 = 402.13*mm + stc_frame_z/2 * cos(th) - stc_frame_w/2 * sin(th);
                G4double dl = (stc_frame_x - stc_frame_w)/2.;
                G4double dz = stc_frame_y / 2. - (stc_active_y_min_[i] + stc_active_y_max_[i]) / 2.0;
                
                G4Box* stc_frameX = new G4Box("stc_frameX", stc_frame_x/2, stc_frame_w/2, stc_frame_z/2);
                G4Box* stc_frameY = new G4Box("stc_frameY", stc_frame_w/2, stc_frame_y/2 - stc_frame_w, stc_frame_z/2);
                stc_frameX_log_[i] = new G4LogicalVolume(stc_frameX, Al, "STC_frameX_log", 0, 0, 0);
                stc_frameY_log_[i] = new G4LogicalVolume(stc_frameY, Al, "STC_frameY_log", 0, 0, 0);
                new G4PVPlacement(zRot, G4ThreeVector(sig*(stc_x0 + dl*cos(th)), -stc_frame_y/2 + stc_frame_w/2 + dz, stc_z0 - dl*sin(th)), stc_frameX_log_[i], (GetLabel(i)+"_frameX1_phys").c_str(), mother_volume_, false, 0, false);
                new G4PVPlacement(zRot, G4ThreeVector(sig*(stc_x0 + dl*cos(th)), +stc_frame_y/2 - stc_frame_w/2 + dz, stc_z0 - dl*sin(th)), stc_frameX_log_[i], (GetLabel(i)+"_frameX2_phys").c_str(), mother_volume_, false, 0, false);
                new G4PVPlacement(zRot, G4ThreeVector(sig*stc_x0, dz, stc_z0),stc_frameY_log_[i], (GetLabel(i)+"_frameY1_phys").c_str(), mother_volume_, false, 0, false);
                new G4PVPlacement(zRot, G4ThreeVector(sig*(stc_x0 + 2*dl*cos(th)), dz, stc_z0 - 2*dl*sin(th)), stc_frameY_log_[i], (GetLabel(i)+"_frameY2_phys").c_str(), mother_volume_, false, 0, false);
                G4Colour stc_frame_colour;
                G4Colour::GetColour("stc_frame", stc_frame_colour);
                G4VisAttributes* stc_frame_logVisAtt = new G4VisAttributes(stc_frame_colour);
                stc_frameY_log_[i]->SetVisAttributes(stc_frame_logVisAtt);
                stc_frameX_log_[i]->SetVisAttributes(stc_frame_logVisAtt);
            }
        }
        
        
        // temporarily including SST support rod
        
//        G4Material* rod_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
//        double rod_diameter = 20*mm;
//        double rod_height = 60*cm;
//        G4Colour pmt_col;
//        G4Colour::GetColour("pmt", pmt_col);
//        G4Tubs* pmt = new G4Tubs((label_ + "_sst_rod_geom").c_str(), 0, rod_diameter/2, rod_height/2, 0, 360*deg);
//        G4LogicalVolume* rod_log = new G4LogicalVolume(pmt, rod_material, (label_ + "_sst_rod_log").c_str(), 0,0,0);
//        rod_log->SetVisAttributes(new G4VisAttributes(pmt_col));
//        G4RotationMatrix* Rot = new G4RotationMatrix;
//        Rot->rotateX(90*deg);
//        new G4PVPlacement(Rot, G4ThreeVector(91.31*mm,0,832.432*mm), rod_log, (GetLabel(i)+"_sst_rod").c_str(), mother_volume_, false, 0, false);
//        new G4PVPlacement(Rot, G4ThreeVector(86.475*mm,0,540.121*mm), rod_log, (GetLabel(i)+"_sst_rod").c_str(), mother_volume_, false, 0, false);
//        new G4PVPlacement(Rot, G4ThreeVector(-91.31*mm,0,832.432*mm), rod_log, (GetLabel(i)+"_sst_rod").c_str(), mother_volume_, false, 0, false);
//        new G4PVPlacement(Rot, G4ThreeVector(-86.475*mm,0,540.121*mm), rod_log, (GetLabel(i)+"_sst_rod").c_str(), mother_volume_, false, 0, false);
        
        //------------------------- Detector
        
        G4Box* stc_detector = new G4Box("stc_chamber", stc_active_x_[i]/2, stc_active_y_[i]/2, stc_full_chamber_width_/2);
        
        stc_assembly_log_[i] = new G4LogicalVolume(stc_detector, Air, (GetLabel(i)+"_log").c_str(), 0, 0, 0);
        if (i == kFront) {
            new G4PVPlacement(NULL, G4ThreeVector(0,0,-stc_gap/2), stc_assembly_log_[i], (GetLabel(i)+"_phys").c_str(), rf_shield_front2_log, false, 0, false);
        } else {
            new G4PVPlacement(NULL, G4ThreeVector(0,0,0), stc_assembly_log_[i], (GetLabel(i)+"_phys").c_str(), rf_shield_rear2_log, false, 0, false);
        }
//        G4VisAttributes* stc_assembly_logVisAtt = new G4VisAttributes(G4Colour::Red());
//        stc_assembly_logVisAtt->SetForceWireframe(true);
//        stc_assembly_log_[i]    ->SetVisAttributes(stc_assembly_logVisAtt);
        stc_assembly_log_[i]->SetVisAttributes(G4VisAttributes::Invisible);
        

        // wire planes

        G4double dz = stc_plane_z_;
        G4Box* stt_plane1 = new G4Box("stt_box1", stc_active_x_[i]/2, stc_active_y_[i]/2, dz/2);
        dz -= t_mylar_;
        G4Box* stt_plane2 = new G4Box("stt_box2", stc_active_x_[i]/2, stc_active_y_[i]/2, dz/2);
        
        G4LogicalVolume* stc_plane1_log = new G4LogicalVolume(stt_plane1, Mylar, (GetLabel(i) + "_mylar").c_str(), 0, 0, 0);
        G4LogicalVolume* stc_plane2_log = new G4LogicalVolume(stt_plane2, STCGas, (GetLabel(i) + "_gas").c_str(), 0, 0, 0);
        
        stc_plane1_log->SetVisAttributes(new G4VisAttributes(stc_colour));
        stc_plane2_log->SetVisAttributes(new G4VisAttributes(gas_colour));
        
        new G4PVPlacement(0, G4ThreeVector(0,0,0), stc_plane2_log, Form("%s_gas", GetLabel(i).c_str()), stc_plane1_log, false, 0);
        
        stc_plane_log_[i] = stc_plane2_log;
        
        int n = 0;
        for (int ip = 0; ip < no_of_xplanes_; ip++) {
            new G4PVPlacement(0, G4ThreeVector(0, 0, -dz_offset_ + (ip - no_of_xplanes_ + 1)*dz_planes_), stc_plane1_log, Form("%s_X_plane%d", GetLabel(i).c_str(), ip), stc_assembly_log_[i], false, n++);
        }
        for (int ip = 0; ip < no_of_yplanes_; ip++) {
            new G4PVPlacement(0, G4ThreeVector(0, 0, dz_offset_ + ip*dz_planes_), stc_plane1_log, Form("%s_Y_plane%d", GetLabel(i).c_str(), ip), stc_assembly_log_[i], false, n++);
        }
    }
}


void g4PSISTC::SetSD(G4SDManager *SDman) {
    
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        if (stc_plane_log_[i]) {
            G4String STCSDname = "g4PSI/STC/" + GetLabel(i);
            G4VSensitiveDetector* STCSD = SDman->FindSensitiveDetector(STCSDname);
            if (STCSD == NULL) {
                STCSD = SD_[i] = new g4PSITrackerSD(STCSDname, GetLabel(i) + "_Collection" );
                SD_[i]->SetDepth(1);
                SDman->AddNewDetector( STCSD );
            };
            stc_plane_log_[i]->SetSensitiveDetector( STCSD );
        }
        if (stc_frameY_log_[i]) {
            G4String STCFSDname = "g4PSI/STCF/" + GetLabel(i);
            G4VSensitiveDetector* STCFSD = SDman->FindSensitiveDetector(STCFSDname);
            if (STCFSD == NULL) {
                STCFSD = SDF_[i] = new g4PSITrackerSD(STCFSDname, GetLabel(i) + "_Collection" );
                SDman->AddNewDetector( STCFSD );
            };
            stc_frameY_log_[i]->SetSensitiveDetector( STCFSD );
        }
    }
}


void g4PSISTC::Write() {
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        TVectorD v1(4);
        v1[0] = stc_r_[i];
        v1[1] = stc_angle_;
        v1[2] = stc_active_x_[i];
        v1[3] = stc_active_y_[i];
        v1.Write(GetLabel(i));
    }
}

void g4PSISTC::InitTree(TTree *T) {
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        if (SD_[i]) SD_[i]->InitTree(T);
        if (SDF_[i]) SDF_[i]->InitTree(T, "F");
    }
}


void g4PSISTC::DeleteEventData() {
    for (int i = kFront; i < NUMBER_OF_STT_CHAMBERS; i++) {
        if (SD_[i]) SD_[i]->DeleteEventData();
        if (SDF_[i]) SDF_[i]->DeleteEventData();
    }
}


G4String g4PSISTC::GetLabel(int c) {
    std::ostringstream oss;
    oss << stc_label_ << c + 1;
    return oss.str();
}

