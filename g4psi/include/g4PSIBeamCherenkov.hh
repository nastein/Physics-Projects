#ifndef g4PSIBeamCherenkov_h
#define g4PSIBeamCherenkov_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIBeamCherenkov
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

class g4PSIBeamCherenkov:public g4PSIDetectorBase {
    
public:
    g4PSIBeamCherenkov(G4String label, G4double z);
    g4PSIBeamCherenkov(G4String label, G4Material *mat, G4double z,
                       G4double angle, G4double len,
                       G4double sx, G4double sy);
    g4PSIBeamCherenkov(G4String label, G4Material *mat, G4double z,
                       G4double angle, G4double len,
                       G4int N1, G4double sx1, G4int N2, G4double sx2, G4double sy,
                       G4double roll, G4double gap);
    ~g4PSIBeamCherenkov();
    
public:
    enum BCCMode {
        kBCC_vMode,
        kBCC_SinglePaddle,
        kBCC_MultiplePaddles,
        kBCC_SISH
    };
    
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return label_;};
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
    const G4double sipm_height_ = 4*CLHEP::mm;
    const G4double sipm_width_ = 3*CLHEP::mm;
    
    BCCMode bcc_mode_;
    G4double bcc_z_;
    G4int bcc_N1_;
    G4int bcc_N2_;
    G4double bcc_len_;
    G4double bcc_paddle_width1_;
    G4double bcc_paddle_width2_;
    G4double bcc_paddle_thickness_;
    G4double bcc_angle_yaw_;
    G4Material *bcc_mat_;
    G4double bcc_gap_;
    G4double bcc_angle_roll_;
    G4LogicalVolume *bcc_log1_;
    G4LogicalVolume *bcc_log2_;
    G4LogicalVolume *bcc_sipm1_log_;
    G4LogicalVolume *bcc_sipm2_log_;
    g4PSIScintillatorSD* SD_;
    
    void place_paddle(G4LogicalVolume *paddle, G4LogicalVolume* sipm, G4LogicalVolume *mother, double &x, int &copyID, double step);
};

#endif
