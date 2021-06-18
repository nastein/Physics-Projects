#ifndef g4PSI_h
#define g4PSISTC_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSISTC
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITrackerSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"

class g4PSISTC : public g4PSIDetectorBase {

public:

    enum STTChamber {
        kFront,
        kRear,
        NUMBER_OF_STT_CHAMBERS
    };

    g4PSISTC(G4String label, G4double angle, STTChamber stc_part);
    g4PSISTC(G4String label, G4double angle);
    ~g4PSISTC();

public:
    void Placement();
//    G4LogicalVolume* GetLog() {return stc_detector_log_;};
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return stc_label_;};
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();

private:
    void STT_Init(G4String label, G4double angle);
    G4String GetLabel(int c);
    const double straw_spacing_ = 1.01*CLHEP::cm;
    const double straw_wall_thickness_ = 30.*CLHEP::um;
    const double dz_planes_ = 0.81*CLHEP::cm;
    const double dz_offset_ = 1.01*CLHEP::cm/2.0;
    G4double stc_full_chamber_width_;
    
    G4double t_mylar_;
    G4double t_arco2_;
    
    G4String stc_label_;
    G4double stc_angle_;
//    STTChamber stc_part_;
    
    int NUMBER_OF_STRAWS_IN_YPLANE_[NUMBER_OF_STT_CHAMBERS];
    int NUMBER_OF_STRAWS_IN_XPLANE_[NUMBER_OF_STT_CHAMBERS]; 
    G4double straw_x_length_[NUMBER_OF_STT_CHAMBERS];
    G4double straw_y_length_[NUMBER_OF_STT_CHAMBERS];

    static const int no_of_xplanes_ = 5;
    static const int no_of_yplanes_ = 5;

    G4bool stc_on_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_r_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_y_min_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_y_max_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_x_min_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_x_max_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_x_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_active_y_[NUMBER_OF_STT_CHAMBERS];
    G4double stc_plane_z_;
    
    G4LogicalVolume *stc_plane_log_[NUMBER_OF_STT_CHAMBERS];
    G4LogicalVolume *stc_assembly_log_[NUMBER_OF_STT_CHAMBERS];
    G4LogicalVolume *stc_frameX_log_[NUMBER_OF_STT_CHAMBERS];
    G4LogicalVolume *stc_frameY_log_[NUMBER_OF_STT_CHAMBERS];
    G4LogicalVolume *stc_tube_y_log_[NUMBER_OF_STT_CHAMBERS];
    G4LogicalVolume *stc_tube_x_log_[NUMBER_OF_STT_CHAMBERS];
    g4PSITrackerSD* SD_[NUMBER_OF_STT_CHAMBERS];
    g4PSITrackerSD* SDF_[NUMBER_OF_STT_CHAMBERS];
};

#endif
