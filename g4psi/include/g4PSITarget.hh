//
//  g4PSITarget.hh
//  g4PSI
//
//  Created by Steffen Strauch on 6/5/15.
//
//

#ifndef g4PSI_g4PSITarget_hh
#define g4PSI_g4PSITarget_hh

#define OPAQUE
//#define BAIL

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITargetSD.hh"
#include "G4PVPlacement.hh"

using namespace CLHEP;

class g4PSITarget: public g4PSIDetectorBase {

public:
    
    g4PSITarget(G4String label = "TARGET", G4double target_radius = 30*mm, G4double kapton_thickness = 120*um, G4double pipe_distance = 60*mm, G4double pipe_theta = 25*deg, G4double entrance_thickness = 200*um);

    void Placement();
    double GetZPos();
    void Write();
    void Info();
    void SetSD(G4SDManager *SDman);
    void InitTree(TTree *T);
    void DeleteEventData();
    
    double GetUpstreamPos();
    double GetTargetRadius();
    
private:
    void MakeWindow(G4String s,G4VSolid* &window, G4double w, G4double h, G4double t, G4double r);
    void MakeWindow(G4String s, G4LogicalVolume* &window_log, G4Material* mat, G4double w, G4double h, G4double t, G4double r);
    
    G4LogicalVolume* mother_for_target_log_;
    G4UserLimits* StepLimit_;            // pointer to user step limits
    g4PSITargetSD* SD_target_;
    
    G4double min_p_;
    G4double min_theta_;

    G4LogicalVolume* target_log_;
    G4LogicalVolume* target_cell_film_log_;
    G4LogicalVolume* target_si_log_;
    G4LogicalVolume* target_base_top_log_;
    G4LogicalVolume* target_base_bottom_log_;
    G4LogicalVolume* target_nipple_log_;
    G4LogicalVolume* target_pipe_v_log_;
    G4LogicalVolume* target_pipe_h_log_;
    G4LogicalVolume* target_probe_t_log_;
    G4LogicalVolume* target_probe_p_log_;
    G4LogicalVolume* chamber_full_log_;
    G4LogicalVolume* chamber_full_vac_log_;
    G4LogicalVolume* chamber_entrance_window_log_;
    G4LogicalVolume* chamber_exit_window_log_;
    G4LogicalVolume* window_post_log_;
//    G4LogicalVolume* scintillator_bar_A_log_;
//    G4LogicalVolume* scintillator_bar_B_log_;
    G4LogicalVolume* chamber_side_exit_window_log_;
    G4LogicalVolume* chamber_back_exit_window_log_;
    
    // geometry - scattering chamber
    
    G4double chamber_outer_radius_;
    G4double chamber_lid_radius_;
    G4double chamber_lid_height_;
    G4double chamber_inner_radius_;
    G4double chamber_height_;
    G4double top_part_chamber_height_;
    
    G4double bottom_flange_inner_radius;
    G4double bottom_flange_height_;
    
    G4double exit_window_max_angle_;
    G4double exit_window_height_;
    G4double chamber_exit_window_thickness_;
    G4double side_exit_window_width_;
    G4double back_exit_window_width_;
    //  G4double side_wall_width_;
    //  G4double back_wall_width_;
    G4double side_wall_center_to_front_;
    G4double back_wall_center_to_front_;
    
    G4double window_offset_;
    
    G4double entrance_window_zpos_;
    G4double entrance_window_ypos_;
    G4double entrance_flange_vac_closer_height_;
    G4double entrance_flange_vac_farther_height_;
    G4double entrance_tube_inner_radius_;
    G4double entrance_tube_outer_radius_;
    G4double entrance_window_thickness_;
    G4double exit_window_thickness_;
    G4double kapton_sheet_thickness_;
    G4double entrance_tube_height_;
    
    G4double height_above_exit_window = 0;
    G4double height_below_exit_window = 0;
    
    // geometry - target
    
    G4double target_dx_;
    G4double target_dy_;
    G4double target_dz_;
    
    G4double target_entrance_cap_thickness_;
    G4double target_exit_cap_thickness_;
    G4double target_wall_thickness_;
    
    G4double target_mylar_thickness_;
    G4double target_flask_mylar_gap_;
    
    G4double chamberPos_y_offset;

    //UM added class variables 
//    G4LogicalVolume* bottomcap_log_;
//    G4LogicalVolume* topcap_log_;
    G4LogicalVolume* rightframe_log_;
    G4LogicalVolume* leftframe_log_;
    G4LogicalVolume* Tpipe_log_;
    G4LogicalVolume* Tconnect_log_;
    G4LogicalVolume* Bpipe_log_;
    G4LogicalVolume* Bconnect_log_;
    G4LogicalVolume* mylar1_log_;
    G4LogicalVolume* mylar2_log_;
    G4LogicalVolume* mylar3_log_;

    G4double pipe_dist;
    G4double pipe_angle;
    G4double target_film_thickness;
    G4double radius;
    G4double base_oneOutterR;
    G4double base_twoOutterR;
    G4double base_one_thickness;
    G4double cellHeight;
    G4double capHeight;
    G4double cell_distance;



};

#endif
