#include "g4PSIBeamMonitor.hh"
#include "g4PSIScintillatorSD.hh"
#include "g4PSIDetectorParts.hh"
#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TVectorD.h"
#include "G4GenericTrap.hh"
#include "G4Trd.hh"
#include "G4PVDivision.hh"

using namespace CLHEP;

g4PSIBeamMonitor::g4PSIBeamMonitor(G4String label, G4double z)
: g4PSIDetectorBase(label) {
    label_ = label;
    posz_ = z;
    gap_bar_ = 0.5 * mm;
    gap_plate_ = 0.5 * mm;
    m_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    N_ = 8;
    sizx_ = 12 * mm;
    sizy_ = 300 * mm;
    sizz_ = 3 * mm;
    sizw_ = 50 * mm;
    
    shiftX_ = (sizx_ + gap_plate_) / 2.0;
    shiftZ_ = 0.5*cm;
    
    nSC_vertical_bars_ = 2;
    nSC_horizontal_bars_ = 0;
    
    len_v_ = sizy_ + nSC_horizontal_bars_ * 2 * (sizw_ + gap_bar_);
    len_h_ = nSC_vertical_bars_ * 2 * (sizw_ + gap_bar_) + ((N_-1)*sizx_ + N_ * gap_plate_);
    
    SD_ = NULL;
}


g4PSIBeamMonitor::~g4PSIBeamMonitor() {
}

double g4PSIBeamMonitor::GetZPos() {return posz_;}

void g4PSIBeamMonitor::Info() {
    InfoTitle("Beamline Monitor");
    InfoParDouble("z", posz_/cm, "cm");
    InfoParInt("Number of sc plates front", N_);
    InfoParInt("Number of sc plates back", N_ - 1);
    InfoParDouble("Back plane offset X", shiftX_);
    InfoParDouble("Back plane offset Z", shiftZ_);
    InfoParDouble("Gap between plates", gap_plate_, "mm");
    InfoParDouble("Gap between bars", gap_bar_, "mm");
    InfoPar3Double("Plate dimensions", sizx_, sizy_, sizz_, "mm");
    InfoPar2Double("SC bar cross section", sizw_, sizw_, "mm");
    InfoParDouble("horizontal bar dimension", len_h_, "mm");
    InfoParDouble("vertical bar dimension", len_v_, "mm");
}


void g4PSIBeamMonitor::Placement() {
    
    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* silicon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

    double groupX = len_h_/2;
    double groupY = len_v_/2;
    double groupZ = (shiftZ_+ 2.*sizw_) / 2.;
    if (sipm_width_ > sizz_) groupZ += sipm_width_ - sizz_;
    if (nSC_horizontal_bars_ * sizw_ < sipm_height_) groupY += sipm_height_;
    
    double zoffset = -groupZ + sizz_/2.;
    
    G4Box* bcc_group = new G4Box("bm_group", groupX, groupY, groupZ);
    G4LogicalVolume *bcc_log = new G4LogicalVolume(bcc_group, air, (label_+"_log").c_str(), 0, 0, 0);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, posz_),
                      bcc_log, (label_+"_phys").c_str(), mother_volume_, false, 0);
    
    bcc_log->SetVisAttributes (G4VisAttributes::Invisible);
//    G4VisAttributes* logVisAtt = new G4VisAttributes();
//    logVisAtt->SetForceWireframe(true);
//    bcc_log->SetVisAttributes(logVisAtt);
    
    
    // Place single beam monitor plates
    // --------------------------------
    
    G4Box* bm_plate = new G4Box("bm_plate", sizx_/2, sizy_/2, sizz_/2);
    bm_log1_ = new G4LogicalVolume(bm_plate, m_, (label_+"_plate_log").c_str(), 0, 0, 0);
    G4Box* bcc_sipm = new G4Box("bcc_sipm", sizx_/2, sipm_height_/2, sipm_width_/2);
    bm_sipm_log_ = new G4LogicalVolume(bcc_sipm, silicon, (label_+"_sipm_log").c_str(), 0, 0, 0);
    G4Colour bcc_col;
    G4Colour::GetColour("pmt", bcc_col);
    bm_sipm_log_->SetVisAttributes(new G4VisAttributes(bcc_col));

    // CopyID:
    //
    // 0, 2, 4, ... for the upstream row of N_ paddles
    // 1, 3, 5, ... for the downstream row of N_ - 1 paddles
    // in total 2*N_1 - 1 scintillator paddles
    //
    //  -"- + nPaddles for the top SiPM
    //  -"- + 2*nPaddels for the bottom SiPM
    //
    char pName[60];
    G4int nPaddles = 2*N_ - 1;
    for (G4int i = 0; i < nPaddles; i++) {
        G4double x = (i+1-N_)*(sizx_+gap_plate_)/2. + 0*(i%2?0:shiftX_);
        G4double y = sizy_/2+sipm_height_/2;
        G4double z = zoffset+(i%2?shiftZ_:0);
        sprintf(pName, "%s_p%02d", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(x,0,z), bm_log1_, pName, bcc_log, false, i);
        sprintf(pName, "%s_sipm1_%02d", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(x,+y,z), bm_sipm_log_, pName, bcc_log, false, nPaddles + i);
        sprintf(pName, "%s_sipm2_%02d", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(x,-y,z), bm_sipm_log_, pName, bcc_log, false, 2*nPaddles + i);
    }
    
    // Place vertical SC bars
    // ----------------------

    G4int n = 3*nPaddles;
    double zoffsetSC_vertical = zoffset + shiftZ_ + sizw_/2. - sizz_/2.;
    G4Box* bm_sc_vertical = new G4Box("bm_sc", sizw_/2, len_v_/2, sizw_/2);
    bm_log2_ = new G4LogicalVolume(bm_sc_vertical, m_, (label_+"_sc_log").c_str(), 0, 0, 0);
    for (G4int i = 0; i < nSC_vertical_bars_; i++) {
        sprintf(pName, "%s_v_sc%02d", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(-(2*i+1)*(sizw_+ gap_bar_)/2.+(1-N_)*(sizx_+gap_plate_)/2.,0, zoffsetSC_vertical), bm_log2_, pName, bcc_log, false, n++);
        sprintf(pName, "%s_v_sc%02d", label_.c_str(), nSC_vertical_bars_ + i);
        new G4PVPlacement(0, G4ThreeVector(+(2*i+1)*(sizw_+ gap_bar_)/2.+(N_-1)*(sizx_+gap_plate_)/2.,0, zoffsetSC_vertical), bm_log2_, pName, bcc_log, false, n++);
    }

    // Place horizontal SC bars
    // ------------------------
    
    double zoffsetSC_horizontal = zoffset + shiftZ_ + 3*sizw_/2. - sizz_/2.;
    G4Box* bm_sc_horizontal = new G4Box("bm_sc", len_h_/2, sizw_/2, sizw_/2);
    bm_log3_ = new G4LogicalVolume(bm_sc_horizontal, m_, (label_+"_sc_log").c_str(), 0, 0, 0);
    for (G4int i = 0; i < nSC_horizontal_bars_; i++) {
        sprintf(pName, "%s_h_sc%02d", label_.c_str(), i);
        new G4PVPlacement(0, G4ThreeVector(0,-(2*i+1)*(sizw_+ gap_bar_)/2.-sizy_/2., zoffsetSC_horizontal), bm_log3_, pName, bcc_log, false, n++);
        sprintf(pName, "%s_h_sc%02d", label_.c_str(), nSC_vertical_bars_ + i);
        new G4PVPlacement(0, G4ThreeVector(0,+(2*i+1)*(sizw_+ gap_bar_)/2.+sizy_/2., zoffsetSC_horizontal), bm_log3_, pName, bcc_log, false, n++);
    }
    
    
    // Visualization
    // -------------
    
    G4Colour bm_col;
    G4Colour::GetColour("bm", bm_col);
    bm_log1_->SetVisAttributes(new G4VisAttributes(bm_col));
    bm_log2_->SetVisAttributes(new G4VisAttributes(bm_col));
    bm_log3_->SetVisAttributes(new G4VisAttributes(bm_col));
}


void g4PSIBeamMonitor::SetSD(G4SDManager *SDman) {
    if (bm_log1_ || bm_log2_ || bm_log3_ || bm_sipm_log_) {
        G4int nCells = 0;
        if (bm_log1_) nCells += 2*N_ - 1;
        if (bm_sipm_log_) nCells += 2*(2*N_ - 1);
        if (bm_log2_) nCells += 2*nSC_vertical_bars_;
        if (bm_log3_) nCells += 2*nSC_horizontal_bars_;
        
        G4String BeamMonitorSDname = "g4PSI/BeamMonitor/" + label_;
        G4VSensitiveDetector* bm_SD = SDman->FindSensitiveDetector(BeamMonitorSDname);
        if (bm_SD == NULL) {
            bm_SD = SD_ = new g4PSIScintillatorSD(BeamMonitorSDname, nCells, label_ + "_Collection" );
            SDman->AddNewDetector( bm_SD );
        };
        
        if (bm_log1_) bm_log1_->SetSensitiveDetector(bm_SD);
        if (bm_log2_) bm_log2_->SetSensitiveDetector(bm_SD);
        if (bm_log3_) bm_log3_->SetSensitiveDetector(bm_SD);
        if (bm_sipm_log_) bm_sipm_log_->SetSensitiveDetector(bm_SD);
    }
}


void g4PSIBeamMonitor::Write() {
    TVectorD v1(4);
    v1[0] = posz_;
    v1[1] = N_;
    v1[2] = nSC_vertical_bars_;
    v1[3] = nSC_horizontal_bars_;
    v1.Write(label_);
}

void g4PSIBeamMonitor::InitTree(TTree *T) {
    SD_->InitTree(T);
}

void g4PSIBeamMonitor::DeleteEventData() {
    SD_->DeleteEventData();
}
