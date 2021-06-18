//
//  g4PSITargets.cc
//  g4PSI
//
//  Created by Steffen Strauch on 6/5/15.
//
//

#include "g4PSITarget.hh"
#include "g4PSIDetectorParts.hh"
#include "g4PSIBOptrMultiParticleChangeCrossSection.hh"

#include "G4PVParameterised.hh"
#include "G4NistManager.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"
#include "g4PSIScintillatorSD.hh"

#include "TVectorD.h"

using namespace CLHEP;

#define in (2.54*cm)
#define mil (25.4*um)

g4PSITarget::g4PSITarget(G4String label, G4double target_radius, G4double kapton_thickness, G4double pipe_distance, G4double pipe_theta, G4double entrance_thickness) : g4PSIDetectorBase(label) {
    
    mother_volume_ = NULL;
    StepLimit_ = new G4UserLimits(0.5*mm);
    target_log_ = NULL;
    target_cell_film_log_ = NULL;
    target_base_top_log_ = NULL;
    target_base_bottom_log_ = NULL;
    target_nipple_log_ = NULL;
    target_pipe_v_log_ = NULL;
    target_pipe_h_log_ = NULL;
    target_probe_t_log_ = NULL;
    target_probe_p_log_ = NULL;
    target_si_log_ = NULL;
    
    rightframe_log_ = NULL;
    leftframe_log_ = NULL;
    Tpipe_log_ = NULL;
    Tconnect_log_ = NULL;
    Bpipe_log_ = NULL;
    Bconnect_log_ = NULL;
    window_post_log_ = NULL;
    
    chamber_full_log_ = NULL;
    chamber_full_vac_log_ = NULL;
    chamber_entrance_window_log_ = NULL;
    chamber_exit_window_log_ = NULL;
    chamber_side_exit_window_log_ = NULL;
    chamber_back_exit_window_log_ = NULL;
    
    SD_target_ = NULL;
    min_p_ = 10 * CLHEP::MeV;
    min_theta_ = 10 * CLHEP::deg;
    
    // geometry
    
    top_part_chamber_height_ = 0;
    chamber_outer_radius_ = 0;
    chamber_inner_radius_ = 0;
    chamber_lid_radius_ = 0;
    chamber_lid_height_ = 0;
    chamber_height_ = 0;
    exit_window_max_angle_ = 0;
    exit_window_height_ = 0;
    chamber_exit_window_thickness_ = 0;
    side_exit_window_width_ = 0;
    back_exit_window_width_ = 0;
    
    entrance_window_zpos_ = 0;
    entrance_window_ypos_ = 0;
    entrance_flange_vac_closer_height_ = 0;
    entrance_flange_vac_farther_height_ = 0;
    entrance_tube_inner_radius_ = 0;
    entrance_tube_outer_radius_ = 0;
    kapton_sheet_thickness_ = 0;
    entrance_window_thickness_ = 0;
    window_offset_ = 0;
    
    bottom_flange_inner_radius = 0;
    bottom_flange_height_ = 0;
    
    height_above_exit_window = 0;
    height_below_exit_window = 0;

    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    
    if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type1]) {
        chamber_outer_radius_ = 20 * cm;
        chamber_inner_radius_ = 18 * cm;
        chamber_lid_radius_ = chamber_outer_radius_ + 1 * cm;
        chamber_lid_height_ = 2 * cm;
        chamber_height_ = 60 * cm;
        
        exit_window_max_angle_ = 110.*deg;
        exit_window_height_ = 46*cm;
        chamber_exit_window_thickness_ = 200*um;
        
        entrance_window_zpos_ = -25*cm;
        entrance_tube_inner_radius_ = 3*cm;
        entrance_tube_outer_radius_ = 4*cm;
        entrance_window_thickness_ = 200*um;
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type2]) {
        chamber_outer_radius_ = 7 * cm;
        chamber_inner_radius_ = 6 * cm;
        chamber_lid_radius_ = chamber_outer_radius_ + 1 * cm;
        chamber_lid_height_ = 1 * cm;
        chamber_height_ = 60 * cm;
        
        exit_window_max_angle_ = 110.*deg;
        exit_window_height_ = 15*cm;
        chamber_exit_window_thickness_ = 200*um;
        
        entrance_window_zpos_ = -7*cm;
        entrance_tube_inner_radius_ = 3*cm;
        entrance_tube_outer_radius_ = 4*cm;
        entrance_window_thickness_ = 200*um;
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type4]) {
        chamber_outer_radius_ = 3 * in; // for target chamber, not vessel?
        chamber_inner_radius_ = 2.5 * in;
        chamber_lid_radius_ = chamber_outer_radius_ + 1 * cm; // temporary
        chamber_lid_height_ = 2 * cm; // temporary
        chamber_height_ = 60 * cm;
        
        exit_window_max_angle_ = 110.*deg;
        exit_window_height_ = 15*cm;
        chamber_exit_window_thickness_ = 200*um;
        
        entrance_window_zpos_ = -8.265*cm;
        entrance_tube_inner_radius_ = 2 * in;
        entrance_tube_outer_radius_ = 2.19 * in;
        entrance_window_thickness_ = 200*um; // adjust        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeUM]) {

    chamber_outer_radius_ = 350.004 / 2. * mm; // (*)
    chamber_inner_radius_ = 340.004 / 2. * mm; // (*)
    exit_window_height_ = 365.5 * mm;          // (*)
        
    side_wall_center_to_front_ = 112.099 * mm + 5.000 * mm;
    side_exit_window_width_ = 220.161 * mm;  // (*)
        
    back_wall_center_to_front_ = 160.586 * mm; // (*)
    back_exit_window_width_ = 72 * mm; // (*)
    top_part_chamber_height_ = 197.25 * mm; // (*)
    height_above_exit_window = 20.000 * mm; // (*)
    height_below_exit_window = 84.500 * mm; // (*)
        
    entrance_tube_inner_radius_ = 70.200 / 2 * mm;
    entrance_tube_outer_radius_ = 117.094 / 2 * mm;
    entrance_window_thickness_ = 75 * um; // e-mail Dany Horovitz 3/13/16
    entrance_tube_height_ = 200.000 * mm - chamber_outer_radius_;
    exit_window_thickness_ = 300 * um; // e-mail Dany Horovitz 3/13/16
    
    chamber_height_ = exit_window_height_ + height_above_exit_window + height_below_exit_window + top_part_chamber_height_; // okay
        
    entrance_window_ypos_ = 267.25 * mm - (chamber_height_ - top_part_chamber_height_)/2.0;
    entrance_window_zpos_ = -(chamber_outer_radius_ + entrance_tube_height_ + entrance_window_thickness_ / 2);
        
    // make sure you at least set this parameter; it is used elsewhere:
    //entrance_window_zpos_ = -200 * mm;   // upstream, usually z < 0, position
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru] ||
             DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru2]) {
        
        // parameter marked with "(*)" are taken from
        // file: Camber_ Assembly-Coffee Can.EASM  (March 13, 2016)
        
        G4double chamber_size_reduce_ = 15 *mm; // Parameter to change the size of the scattering chamber
        
        //Nominal outer radius is 350.004/2mm and inner radius is 340.004/2*mm
        chamber_outer_radius_ = 350.004 / 2. * mm - chamber_size_reduce_; // (*)
        chamber_inner_radius_ = 340.004 / 2. * mm - chamber_size_reduce_; // (*)
        
        exit_window_height_ = 2*chamber_outer_radius_;          //was 365.6*mm (*)
        
        side_wall_center_to_front_ = 112.099 * mm + 5.000 * mm;
        side_exit_window_width_ = 220.161 * mm;  // (*)
        
        back_wall_center_to_front_ = 160.586 * mm; // (*)
        back_exit_window_width_ = 72 * mm; // (*)
        top_part_chamber_height_ = 197.25 * mm; // (*)
        height_above_exit_window = 20.000 * mm; // (*)
        height_below_exit_window = 84.500 * mm; // (*)
        
        entrance_tube_inner_radius_ = 70.200 / 2 * mm;
        entrance_tube_outer_radius_ = 117.094 / 2 * mm;
        entrance_window_thickness_ = entrance_thickness > 0 ? entrance_thickness : 75 * um; // e-mail Dany Horovitz 3/13/16
        entrance_tube_height_ = 200.000 * mm - chamber_outer_radius_ - chamber_size_reduce_;
        exit_window_thickness_ = 200*um; //was 200 *um
        
        G4double chamber_addition = 365.5 *mm;
        
        chamber_height_ = chamber_addition + height_above_exit_window + height_below_exit_window + top_part_chamber_height_; // okay
        
        entrance_window_ypos_ = 267.25 * mm - (chamber_height_ - top_part_chamber_height_)/2.0;
        entrance_window_zpos_ = -(chamber_outer_radius_ + entrance_tube_height_ + entrance_window_thickness_ / 2);
    }
    
        else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeTrapezoid]) {
        
        // parameter marked with "(*)" are taken from
        // file: Camber_ Assembly-Coffee Can.EASM  (March 13, 2016)

        G4double chamber_size_reduce_ = 15 *mm; // Parameter to change the size of the scattering chamber
        
        //Nominal outer radius is 350.004/2mm and inner radius is 340.004/2*mm
        chamber_outer_radius_ = 350.004 / 2. * mm - chamber_size_reduce_; // (*)
        chamber_inner_radius_ = 340.004 / 2. * mm - chamber_size_reduce_; // (*)
        
        exit_window_height_ = 2*chamber_outer_radius_;          //was 365.6*mm (*)
        
        side_wall_center_to_front_ = 112.099 * mm + 5.000 * mm;
        side_exit_window_width_ = 220.161 * mm;  // (*)
        
        back_wall_center_to_front_ = 160.586 * mm; // (*)
        back_exit_window_width_ = 72 * mm; // (*)
        top_part_chamber_height_ = 197.25 * mm; // (*)
        height_above_exit_window = 20.000 * mm; // (*)
        height_below_exit_window = 84.500 * mm; // (*)
        
        entrance_tube_inner_radius_ = 70.200 / 2 * mm;
        entrance_tube_outer_radius_ = 117.094 / 2 * mm;
        entrance_window_thickness_ = entrance_thickness > 0 ? entrance_thickness : 75 * um; // e-mail Dany Horovitz 3/13/16
        entrance_tube_height_ = 200.000 * mm - chamber_outer_radius_ - chamber_size_reduce_;
        exit_window_thickness_ = 300*um; //was 200 *um

        G4double chamber_addition = 365.5 *mm;
        
        chamber_height_ = chamber_addition + height_above_exit_window + height_below_exit_window + top_part_chamber_height_; // okay
        
        entrance_window_ypos_ = 267.25 * mm - (chamber_height_ - top_part_chamber_height_)/2.0;
        entrance_window_zpos_ = -(chamber_outer_radius_ + entrance_tube_height_ + entrance_window_thickness_ / 2);
            
            entrance_window_ypos_ = -161.33969999999999;
    }
    
    target_dx_ = 0;  // radius
    target_dy_ = 0;  // radius
    target_dz_ = 0;  // full length
    
    target_entrance_cap_thickness_ = 0;
    target_exit_cap_thickness_ = 0;
    target_wall_thickness_ = 0;
    
    target_mylar_thickness_ = 0.;
    target_flask_mylar_gap_ = 0.;

    pipe_dist = 0;
    pipe_angle = 0;
    target_film_thickness = 0;
    radius = 0;
    base_oneOutterR = 0;
    base_twoOutterR = 0;
    base_one_thickness = 0;
    cellHeight = 0;
    capHeight = 0 ;
    cell_distance = 0;
    
    if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type1] ||
        DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type1]) {
        
        target_dx_ = 2.*cm;
        target_dy_ = 2.*cm;
        target_dz_ = 4.*cm;
        
        target_entrance_cap_thickness_ = 0.0625 * mm;
        target_exit_cap_thickness_ = 0.0625 * mm;
        target_wall_thickness_ = 0.0625 * mm;  // 125/2 um for flask
        
        target_mylar_thickness_ = 0.100 * mm;  // dummy
        target_flask_mylar_gap_ = 2 * cm; // dummy
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type2] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type2]) {
        
        target_dx_ = 3.*cm;
        target_dy_ = 3.*cm;
        target_dz_ = 4.*cm;
        
        target_entrance_cap_thickness_ = 0.0625 * mm;
        target_exit_cap_thickness_ = 0.0625 * mm;
        target_wall_thickness_ = 0.0625 * mm;  // 125/2 um for flask
        
        target_mylar_thickness_ = 0.050 * mm;  // dummy, updated from TypeIA
        target_flask_mylar_gap_ = 1 * cm; // dummy, updated from TypeIA
    }
    
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type3] ||
                  DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type3]) {
        
        target_dx_ = 3.*cm;
        target_dy_ = 3.*cm;
        target_dz_ = 10.*cm;
        
        target_entrance_cap_thickness_ = 0.0625 * mm;
        target_exit_cap_thickness_ = 0.0625 * mm;
        target_wall_thickness_ = 0.0625 * mm;  // 125/2 um for flask
        
        target_mylar_thickness_ = 0.050 * mm;
        target_flask_mylar_gap_ = 1 * cm;
    }
    
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type4] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type4]) {
        
        target_dx_ = 3.*cm;
        target_dy_ = 3.*cm;
        target_dz_ = 6.*cm;
        
        target_entrance_cap_thickness_ = 0.0625 * mm;
        target_exit_cap_thickness_ = 0.0625 * mm;
        target_wall_thickness_ = 0.0625 * mm;  // 125/2 um for flask
        
        target_mylar_thickness_ = 0.050 * mm;  // dummy, updated from TypeIA
        target_flask_mylar_gap_ = 1 * cm; // dummy, updated from TypeIA
    }
    
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeUM]) {

        target_dx_ = target_radius; 
        target_dy_ = target_radius;
        target_dz_ = 100 *mm;

        pipe_angle = pipe_theta;
        pipe_dist = pipe_distance;
        radius = target_radius;
        target_film_thickness = kapton_thickness;
        
    }

    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM2] || 
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeUM2]) {

        pipe_angle = pipe_theta;
        pipe_dist = pipe_distance;
        radius = target_radius;
        target_film_thickness = kapton_thickness;

        cell_distance = 16 *cm;
        base_oneOutterR = 61.5*mm/2;
        base_twoOutterR = 63.5*mm/2;
        capHeight = 10/2 *mm;
        base_one_thickness = 1.7*mm/2;
        cellHeight = 100/2 *mm;

        target_dx_ = target_radius; 
        target_dy_ = target_radius;
        target_dz_ = 100 *mm;

        

    }
    
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru] || DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeJeru]) {
        target_dx_ = 3.*cm;
        target_dy_ = 3.*cm;
        target_dz_ = 4.*cm;
        
        target_entrance_cap_thickness_ = 0.0625 * mm;
        target_exit_cap_thickness_ = 0.0625 * mm;
        target_wall_thickness_ = 0.0625 * mm;  // 125/2 um for flask
        
        target_mylar_thickness_ = 0.050 * mm;  // dummy
        target_flask_mylar_gap_ = 1 * cm; // dummy
    }
    
    //------- Placement of chamber and target in world space. Note: in Indiana configuration, chamber slides up or down depending on which entrance window being used and target remains stationary -------
    
    chamberPos_y_offset = entrance_window_ypos_; //temporary; so that window and target in front of beam

}

double g4PSITarget::GetZPos() {return 0.0;}

double g4PSITarget::GetTargetRadius() {
    
    return target_dx_ > target_dy_ ? target_dx_ : target_dy_;
}

double g4PSITarget::GetUpstreamPos() {
    
    // get most upstream part of the scattering chamber or target
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    
    double zt = 0.0;
    double zs = 0.0;
    if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type1] ||
        DetectorParts->build[g4PSIDetectorParts::kTarget_Type2] ||
        DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type1] ||
        DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type2] ||
        DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru] ||
        DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeJeru]) {
        zt = -(target_dz_/2.0 + target_entrance_cap_thickness_ + target_flask_mylar_gap_ + target_mylar_thickness_);
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeUM]) {
        // Most upstream (z < 0) part of the target. This is not needed if the
        // target is put into a scattering chamber. For testing purposes it might be
        // useful.  Set it to, e.g., -(target radius + wall thickness)
        zt = 0.;
    }
    if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type1] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type2] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type4] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeUM] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeTrapezoid] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru2]) {
        zs = entrance_window_zpos_;
    }
    return zt < zs ? zt : zs;  // it is usually zs < zt < 0
}

void g4PSITarget::Placement() {
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    mother_for_target_log_ = mother_volume_;
    
    // ===============================================================
    //
    //   Build scattering chamber
    //
    // ===============================================================
    
    if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type1] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type2] ||
        DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type4]) {
        
        G4double entrance_tube_height = -entrance_window_zpos_ - sqrt(pow(chamber_inner_radius_,2) - pow(entrance_tube_inner_radius_, 2));
        
        G4Material* Aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
        G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        
        G4Tubs* chamber_lid = new G4Tubs("scat_lid", 0.0, chamber_lid_radius_, chamber_lid_height_/2.0, 0, twopi);
        
        // Aluminum volume, including aluminum chamber walls and entrance flange
        G4RotationMatrix* RotUnion = new G4RotationMatrix;
        RotUnion->rotateY(90*deg);
        G4Tubs* chamber_wall = new G4Tubs("scat_wall", 0.0, chamber_outer_radius_, chamber_height_/2.0, 0, twopi);
        G4Tubs* entrance_wall = new G4Tubs("scat_et_wall", 0.0, entrance_tube_outer_radius_, entrance_tube_height/2.0, 0, twopi);
        G4VSolid* full_al_chamber = new G4UnionSolid("chamber_wall", chamber_wall, entrance_wall,
                                                     RotUnion, G4ThreeVector(entrance_window_zpos_ + entrance_tube_height/2,0,0));
        
        // Vacuum volume for chamber and flange
        G4Tubs* chamber_inner_volume = new G4Tubs("scat_vacuum", 0.0, chamber_inner_radius_, chamber_height_/2.0, 0, twopi);
        G4Tubs* entrance_inner_volume = new G4Tubs("scat_et_vacuum", 0.0, entrance_tube_inner_radius_, entrance_tube_height/2.0, 0, twopi);
        G4Tubs* chamber_exitgap = new G4Tubs("scat_exit_vacuum", chamber_inner_radius_, chamber_outer_radius_,
                                             exit_window_height_/2., -exit_window_max_angle_, 2.*exit_window_max_angle_);
        G4VSolid* chamber_plus_gap = new G4UnionSolid("chamber_plus_gap_vacuum", chamber_inner_volume, chamber_exitgap);
        
        G4VSolid* full_vacuum_chamber = new G4UnionSolid("chamber_vacuum",
                                                         chamber_plus_gap, entrance_inner_volume,
                                                         RotUnion, G4ThreeVector(entrance_window_zpos_ + entrance_tube_height/2,0,0));
        
        // Chamber entrance and exit windows
        G4Tubs* chamber_exit_window
        = new G4Tubs("scat_exit_window",
                     chamber_outer_radius_,
                     chamber_outer_radius_ + chamber_exit_window_thickness_,
                     exit_window_height_/2., -exit_window_max_angle_, 2.*exit_window_max_angle_);
        G4Tubs* chamber_entrance_window
        = new G4Tubs("scat_entrance_window",
                     0, (entrance_tube_inner_radius_+entrance_tube_outer_radius_)/2,
                     entrance_window_thickness_/2., 0, twopi);
        
        G4LogicalVolume *full_al_chamber_log = new G4LogicalVolume(full_al_chamber, Aluminum, "chamber_wall.log");
        G4LogicalVolume *full_vacuum_chamber_log = new G4LogicalVolume(full_vacuum_chamber, Vacuum, "chamber_vacuum.log");
        G4LogicalVolume *chamber_lid_log = new G4LogicalVolume(chamber_lid, Aluminum, "chamber_lid.log");
        G4LogicalVolume *chamber_exit_window_log = new G4LogicalVolume(chamber_exit_window, Kapton, "chamber_exit_window.log");
        G4LogicalVolume *chamber_entrance_window_log = new G4LogicalVolume(chamber_entrance_window, Kapton, "chamber_entrance_window.log");
        
        DetectorParts->AttachBiasingOperator(full_al_chamber_log);
        DetectorParts->AttachBiasingOperator(chamber_entrance_window_log);
        DetectorParts->AttachBiasingOperator(chamber_exit_window_log);
        
        G4RotationMatrix* Rot = new G4RotationMatrix;
        Rot->rotateX(90*deg);
        Rot->rotateZ(90*deg);
        
        new G4PVPlacement(Rot, G4ThreeVector(0,0,0),
                          full_al_chamber_log, "full_al_chamber.vol", mother_volume_, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0,0,0),
                          full_vacuum_chamber_log, "full_vacuum_chamber.vol", full_al_chamber_log, false, 0);
        new G4PVPlacement(Rot, G4ThreeVector(0,0,0),
                          chamber_exit_window_log, "scat_exit_window.vol", mother_volume_, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0,0, entrance_window_zpos_ - entrance_window_thickness_/2),
                          chamber_entrance_window_log, "scat_entrance_window.vol", mother_volume_, false, 0);
        new G4PVPlacement(Rot, G4ThreeVector(0, chamber_height_/2 + chamber_lid_height_/2,0),
                          chamber_lid_log, "scat_top_lid.vol", mother_volume_, false, 0);
        new G4PVPlacement(Rot, G4ThreeVector(0, -chamber_height_/2 - chamber_lid_height_/2,0),
                          chamber_lid_log, "scat_bottom_lid.vol", mother_volume_, false, 0);
        
        G4Colour colour_chamber;
        G4Colour colour_chamber_window;
        G4Colour colour_chamber_air;
        G4Colour::GetColour("scattering_chamber", colour_chamber);
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4Colour::GetColour("scattering_chamber_airgap", colour_chamber_air);
        G4VisAttributes* scatVisAtt = new G4VisAttributes(colour_chamber);
        G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
        G4VisAttributes* scatGVisAtt = new G4VisAttributes(colour_chamber_air);
        
        scatVisAtt->SetVisibility(true);
        scatVisAtt->SetForceLineSegmentsPerCircle(36);
        scatWVisAtt->SetVisibility(true);
        scatWVisAtt->SetForceLineSegmentsPerCircle(50);
        scatGVisAtt->SetVisibility(true);
        
        full_al_chamber_log->SetVisAttributes(scatVisAtt);
        full_vacuum_chamber_log->SetVisAttributes(scatGVisAtt);
        chamber_lid_log->SetVisAttributes(scatVisAtt);
        chamber_exit_window_log->SetVisAttributes(scatWVisAtt);
        chamber_entrance_window_log->SetVisAttributes(scatWVisAtt);
        
        // limit step size within the scattering chamber to 1 mm
        // does not propagate to daughter, shall we use regions?
        
        full_vacuum_chamber_log->SetUserLimits(StepLimit_);
        mother_for_target_log_ = full_vacuum_chamber_log;
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type3]) {
        
        /// \todo rewrite dummy chamber in terms of data member variables.
        
        // dummy scattering chamber
        
        G4double scat_r = 12.5*cm; // previously 3 * inch;
        G4double scat_dr = 0.05 * mm; // width of the wall
        G4double scat_height = 50 * cm; // total height
        G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        
        G4RotationMatrix* Rot = new G4RotationMatrix;
        Rot->rotateX(90*deg);
        Rot->rotateZ(90*deg);
        
        G4Tubs* tub = new G4Tubs("scat_tub", scat_r-scat_dr, scat_r, scat_height/2.0, 0, twopi);
        G4Tubs* innerVol = new G4Tubs("scat_inner_vol", 0.0, scat_r-scat_dr, scat_height/2.0, 0, twopi);
        G4Tubs* lid = new G4Tubs("scat_lid", 0., scat_r, scat_dr/2.0, 0, twopi);
        
        G4VSolid* logo1 = new G4UnionSolid("scat_3", tub, lid, 0, G4ThreeVector(0, 0, (scat_height + scat_dr)/2.));
        G4VSolid* logo2 = new G4UnionSolid("scat_4", logo1, lid, 0, G4ThreeVector(0, 0, -(scat_height + scat_dr)/2.));
        G4LogicalVolume* scat_innerVol_log = new G4LogicalVolume(innerVol, Vacuum, "scat_innerVol_log", 0,0,0);
        G4LogicalVolume* scat_log = new G4LogicalVolume(logo2, Kapton, "scat_chamber_log", 0,0,0);
        
        new G4PVPlacement(Rot, G4ThreeVector(0,0,0),
                          scat_innerVol_log,"scat_innerVol",mother_volume_,false,0);
        new G4PVPlacement(Rot, G4ThreeVector(0,0,0),
                          scat_log,"scat_chamber",mother_volume_,false,0);
        
        G4Colour colour_chamber_window;
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
        //	scatWVisAtt->SetForceWireframe(true);
        scat_innerVol_log->SetVisAttributes(scatWVisAtt);
        scat_log->SetVisAttributes(scatWVisAtt);
        
        // limit step size within the scattering chamber to 1 mm
        // does not propagate to daughter, shall we use regions?
        
        scat_innerVol_log->SetUserLimits(StepLimit_);
        
        mother_for_target_log_ = scat_innerVol_log;
    }

    else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru]) {
        
        //  Materials for chamber and vacuum inside chamber
        //
        G4NistManager* man = G4NistManager::Instance();
        
        G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic"); //  temporary; Is this right, galactic?
        G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
        
        G4double density = 8.06*g/cm3;
        G4int ncomponents = 6;
        G4double fractionmass = 0.;
        G4Material* StainlessSteel = new G4Material("StainlessSteel", density, ncomponents);
        G4Element* C  = man->FindOrBuildElement("C");
        G4Element* Si = man->FindOrBuildElement("Si");
        G4Element* Cr = man->FindOrBuildElement("Cr");
        G4Element* Mn = man->FindOrBuildElement("Mn");
        G4Element* Fe = man->FindOrBuildElement("Fe");
        G4Element* Ni = man->FindOrBuildElement("Ni");
        StainlessSteel->AddElement(C, fractionmass=0.001);
        StainlessSteel->AddElement(Si, fractionmass=0.007);
        StainlessSteel->AddElement(Cr, fractionmass=0.18);
        StainlessSteel->AddElement(Mn, fractionmass=0.01);
        StainlessSteel->AddElement(Fe, fractionmass=0.712);
        StainlessSteel->AddElement(Ni, fractionmass=0.09);
        
        
        /// ------------- Stainless Steel Chamber -------------
        
        G4double bottom_chamber_height = chamber_height_ - top_part_chamber_height_;
        G4double chamber_wall_thickness = chamber_outer_radius_ - chamber_inner_radius_;
        G4double theta_max_outer = acos(side_wall_center_to_front_ / chamber_outer_radius_) + 60 * deg;
        G4double theta_max_inner = acos((side_wall_center_to_front_ - chamber_wall_thickness) / chamber_inner_radius_) + 60 * deg;
       

        G4double back_partial_chamber_trap_zlength = chamber_outer_radius_ * cos(180*deg - theta_max_outer) + back_wall_center_to_front_;
        
        G4double back_partial_chamber_trap_ylength_in = 2 * chamber_outer_radius_ * sin(theta_max_outer);
        G4double back_partial_chamber_trap_ylength_out = 2 * (chamber_outer_radius_ * sin(180*deg - theta_max_outer) -  back_partial_chamber_trap_zlength * tan(30*deg));
        G4double back_partial_chamber_trap_xdisp = -back_partial_chamber_trap_zlength / 2 - chamber_outer_radius_ * cos(theta_max_outer);
        
        const G4double delta = 1 * mm;
        G4Tubs *top_chamber_cyl = new G4Tubs("top_chamber_cylinder", 0, chamber_outer_radius_, 0.5 * top_part_chamber_height_ + delta, 0, 2 * pi);
        G4Tubs *front_partial_chamber_cyl = new G4Tubs("front_partial_chamber_cylinder", 0, chamber_outer_radius_, 0.5 * bottom_chamber_height, -(180 * deg - theta_max_outer), 360 * deg - 2 * theta_max_outer);
        G4Trd *back_partial_chamber_trap =
        new G4Trd("back_partial_chamber_trap",
                  0.5 * bottom_chamber_height, 0.5 * bottom_chamber_height,
                  0.5 * back_partial_chamber_trap_ylength_in, 0.5 * back_partial_chamber_trap_ylength_out,
                  0.5 * back_partial_chamber_trap_zlength);
        
        double entrance_tube_dh = 5 * cm;
        G4Tubs *entrance_tube = new G4Tubs("entrance_tube", 0, entrance_tube_outer_radius_, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);
        
        G4VSolid *full_chamber = new G4UnionSolid("front_partial_chamber_cylinder + top_chamber_cylinder", front_partial_chamber_cyl, top_chamber_cyl, NULL, G4ThreeVector(0, 0, 0.5 * bottom_chamber_height + 0.5 * top_part_chamber_height_ - 0.5 * delta));
        
        G4RotationMatrix *rot_trap = new G4RotationMatrix;
        rot_trap->rotateY(90 * deg);
        full_chamber = new G4UnionSolid("full_chamber + back_partial_chamber_trap", full_chamber, back_partial_chamber_trap, rot_trap, G4ThreeVector(back_partial_chamber_trap_xdisp, 0, 0));
        
        G4RotationMatrix *rot_entrance_window = new G4RotationMatrix;
        rot_entrance_window->rotateY(90 * deg);
        full_chamber = new G4UnionSolid("full_chamber_vac + entrance_tube", full_chamber, entrance_tube, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_));
        
        chamber_full_log_ = new G4LogicalVolume(full_chamber, StainlessSteel, "full_chamber.log");
        
        /// ------------- Vacuum --------------
        
        G4double back_partial_chamber_trap_vac_zlength = chamber_inner_radius_ * cos(180*deg - theta_max_inner) + back_wall_center_to_front_ - chamber_wall_thickness;
        G4double back_partial_chamber_trap_vac_ylength_in = 2 * chamber_inner_radius_ * sin(theta_max_inner);
        G4double back_partial_chamber_trap_vac_ylength_out = 2 * (chamber_inner_radius_ * sin(180*deg - theta_max_inner) -  back_partial_chamber_trap_vac_zlength * tan(30*deg));
        G4double back_partial_chamber_trap_vac_xdisp = -back_partial_chamber_trap_vac_zlength / 2 - chamber_inner_radius_ * cos(theta_max_inner);
        
        G4Tubs *top_chamber_cyl_vac = new G4Tubs("top_chamber_cylinder_vac", 0, chamber_inner_radius_, (top_part_chamber_height_ - 2 *chamber_wall_thickness) / 2, 0, 2 * pi);
        
        G4Tubs *front_partial_chamber_cyl_vac = new G4Tubs("front_partial_chamber_cylinder_vac", 0, chamber_inner_radius_, 0.5 * bottom_chamber_height, -(180 * deg - theta_max_inner), 2 * (180 * deg - theta_max_inner));
        G4Trd *back_partial_chamber_trap_vac =
        new G4Trd("back_partial_chamber_trap_vac",
                  0.5 * bottom_chamber_height, 0.5 * bottom_chamber_height,
                  0.5 * back_partial_chamber_trap_vac_ylength_in, 0.5 * back_partial_chamber_trap_vac_ylength_out,
                  0.5 * back_partial_chamber_trap_vac_zlength);
        G4Box *side_exit_window_vac = new G4Box("side_exit_window_vac", chamber_wall_thickness / 2, side_exit_window_width_ / 2, exit_window_height_ / 2);
        G4Box *back_exit_window_vac = new G4Box("side_exit_window_vac", chamber_wall_thickness / 2, back_exit_window_width_ / 2, exit_window_height_ / 2);
        G4Tubs *entrance_tube_vac = new G4Tubs("entrance_tube_vac", 0, entrance_tube_inner_radius_, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);
        
        G4VSolid *full_chamber_vac = new G4UnionSolid("top_chamber_cylinder_vac + front_partial_chamber_cylinder_vac", front_partial_chamber_cyl_vac, top_chamber_cyl_vac, NULL, G4ThreeVector(0, 0, bottom_chamber_height / 2 + top_part_chamber_height_ / 2 - chamber_wall_thickness));
        full_chamber_vac = new G4UnionSolid("full_chamber_vac + back_partial_chamber_trap_vac", full_chamber_vac, back_partial_chamber_trap_vac, rot_trap, G4ThreeVector(back_partial_chamber_trap_vac_xdisp, 0, 0));
        full_chamber_vac = new G4UnionSolid("full_chamber_vac + entrance window", full_chamber_vac, entrance_tube_vac, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_ - chamber_wall_thickness));
        
        G4double side_window_vac_xdisp = - (side_wall_center_to_front_ - 0.5 * chamber_wall_thickness) * cos(60 * deg);
        G4double side_window_vac_ydisp = - (side_wall_center_to_front_ - 0.5 * chamber_wall_thickness) * sin(60 * deg);
        G4double side_window_vac_zdisp = -bottom_chamber_height/2 + height_below_exit_window + exit_window_height_/2 - chamber_wall_thickness;
        
        G4RotationMatrix *rot_side_exit_window = new G4RotationMatrix;
        rot_side_exit_window->rotateZ(-60*deg);
        
        full_chamber_vac = new G4UnionSolid("full_chamber_vac + side exit window1", full_chamber_vac, side_exit_window_vac, rot_side_exit_window, G4ThreeVector(side_window_vac_xdisp, side_window_vac_ydisp, side_window_vac_zdisp + entrance_window_ypos_));
        
        rot_side_exit_window = new G4RotationMatrix;
        rot_side_exit_window->rotateZ(60*deg);
        side_window_vac_ydisp *= -1;
        
        full_chamber_vac = new G4UnionSolid("full_chamber_vac + side exit window2", full_chamber_vac, side_exit_window_vac, rot_side_exit_window, G4ThreeVector(side_window_vac_xdisp, side_window_vac_ydisp, side_window_vac_zdisp + entrance_window_ypos_));
        full_chamber_vac = new G4UnionSolid("full_chamber_vac + back_exit_window", full_chamber_vac, back_exit_window_vac, NULL, G4ThreeVector(-back_wall_center_to_front_ + chamber_wall_thickness / 2, 0, side_window_vac_zdisp + entrance_window_ypos_));
        
        chamber_full_vac_log_ = new G4LogicalVolume(full_chamber_vac, Vacuum, "full_chamber_vac.log");
        
        /// ------------ Kapton windows --------------
        
        G4Tubs *entrance_window = new G4Tubs("entrance_window", 0, entrance_tube_inner_radius_, entrance_window_thickness_ / 2, 0, 2 * pi);
        G4Box *side_exit_window = new G4Box("side_exit_window", side_exit_window_width_ / 2, exit_window_height_ / 2, exit_window_thickness_ / 2);
        G4Box *back_exit_window = new G4Box("back_exit_window", back_exit_window_width_ / 2, exit_window_height_ / 2, exit_window_thickness_ / 2);
        
        chamber_entrance_window_log_ = new G4LogicalVolume(entrance_window, Kapton, "entrance_window.log");
        chamber_side_exit_window_log_ = new G4LogicalVolume(side_exit_window, Kapton, "side_exit_window.log");
        chamber_back_exit_window_log_ = new G4LogicalVolume(back_exit_window, Kapton, "back_exit_window.log");
        
        // ------ placements
        
        G4RotationMatrix *rot_chamber = new G4RotationMatrix;
        rot_chamber->rotateY(90 * deg);
        rot_chamber->rotateX(90 * deg);
        rot_chamber->rotateZ(180 * deg);
        

        new G4PVPlacement(rot_chamber, G4ThreeVector(0,-entrance_window_ypos_,0), chamber_full_log_, "chamber_full", mother_volume_, false, 0);
        
        new G4PVPlacement(NULL, G4ThreeVector(0,0,chamber_wall_thickness), chamber_full_vac_log_, "chamber_full_vac", chamber_full_log_, false, 1);
        new G4PVPlacement(NULL, G4ThreeVector(0, 0, entrance_window_zpos_), chamber_entrance_window_log_, "chamber_entrance_window", mother_volume_, false, 2);
        
        rot_side_exit_window = new G4RotationMatrix;
        rot_side_exit_window->rotateY(-60 * deg);
        G4double side_window_zdisp = (side_wall_center_to_front_ + 0.5 * exit_window_thickness_) * cos(60 * deg);
        G4double side_window_xdisp = (side_wall_center_to_front_ + 0.5 * exit_window_thickness_) * sin(60 * deg);
        
        new G4PVPlacement(rot_side_exit_window, G4ThreeVector(side_window_xdisp, 0, side_window_zdisp), chamber_side_exit_window_log_, "chamber_side_exit_window", mother_volume_, false, 3);
        
        rot_side_exit_window = new G4RotationMatrix;
        rot_side_exit_window->rotateY(60 * deg);
        side_window_xdisp *= -1;
        
        new G4PVPlacement(rot_side_exit_window, G4ThreeVector(side_window_xdisp, 0, side_window_zdisp), chamber_side_exit_window_log_, "chamber_side_exit_window", mother_volume_, false, 4);
        
        new G4PVPlacement(NULL, G4ThreeVector(0, 0, back_wall_center_to_front_ + exit_window_thickness_ / 2), chamber_back_exit_window_log_, "chamber_back_exit_window", mother_volume_, false, 5);
        
        // ------ vis
        
        G4Colour colour_chamber_window;
        G4Colour colour_chamber;
        G4Colour colour_chamber_vacuum;
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4Colour::GetColour("scattering_chamber", colour_chamber);
        G4Colour::GetColour("scattering_chamber_airgap", colour_chamber_vacuum);
        G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
        G4VisAttributes* scatVisAtt = new G4VisAttributes(colour_chamber);
        G4VisAttributes* scatGVisAtt = new G4VisAttributes(colour_chamber_vacuum);
        
        // Set to false so I could see the inside of the target
        scatGVisAtt->SetVisibility(true);
        scatVisAtt->SetVisibility(true);
        scatWVisAtt->SetVisibility(true);
        
        chamber_full_log_->SetVisAttributes(scatVisAtt);
        chamber_full_vac_log_->SetVisAttributes(scatGVisAtt);
        chamber_side_exit_window_log_->SetVisAttributes(scatWVisAtt);
        chamber_back_exit_window_log_->SetVisAttributes(scatWVisAtt);
        chamber_entrance_window_log_->SetVisAttributes(scatWVisAtt);
        
        DetectorParts->AttachBiasingOperator(chamber_full_log_);
        DetectorParts->AttachBiasingOperator(chamber_side_exit_window_log_);
        DetectorParts->AttachBiasingOperator(chamber_back_exit_window_log_);
        DetectorParts->AttachBiasingOperator(chamber_entrance_window_log_);
        
        chamber_full_log_->SetUserLimits(StepLimit_);
        mother_for_target_log_ = chamber_full_vac_log_;
    } else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeJeru2]) {
        
        //  kScatteringChamber_TypeJeru2 is a cylindrical version of the trapezoidal TypeJeru chamber
        
        //  Materials for chamber and vacuum inside chamber
        G4NistManager* man = G4NistManager::Instance();

        G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic"); 
        G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
        G4Material* SiO2   = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        G4Material* Sc_material = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        G4Element* C  = man->FindOrBuildElement("C");
        G4Element* Si = man->FindOrBuildElement("Si");
        G4Element* Cr = man->FindOrBuildElement("Cr");
        G4Element* Mn = man->FindOrBuildElement("Mn");
        G4Element* Fe = man->FindOrBuildElement("Fe");
        G4Element* Ni = man->FindOrBuildElement("Ni");
        G4Element* H  = man ->FindOrBuildElement("H");

        //FR4 Glass + Epoxy 
        G4int ncomp, natoms;
        G4double densityEpox = 1.2*g/cm3;
        G4double densityFR4 = 1.86*g/cm3;
        G4double frac;
        G4Material* Epoxy = new G4Material("Epoxy" , densityEpox, ncomp=2);
        Epoxy->AddElement(H, natoms=2);
        Epoxy->AddElement(C, natoms=2);
        G4Material* G10 = new G4Material("G10"  , densityFR4, ncomp=2);
        G10->AddMaterial(SiO2, frac=0.528);
        G10->AddMaterial(Epoxy, frac=0.472);

        //Stainless steel
        G4double density = 8.06*g/cm3;
        G4int ncomponents = 6;
        G4double fractionmass = 0.;
        G4Material* StainlessSteel = new G4Material("StainlessSteel", density, ncomponents);
        StainlessSteel->AddElement(C, fractionmass=0.001);
        StainlessSteel->AddElement(Si, fractionmass=0.007);
        StainlessSteel->AddElement(Cr, fractionmass=0.18);
        StainlessSteel->AddElement(Mn, fractionmass=0.01);
        StainlessSteel->AddElement(Fe, fractionmass=0.712);
        StainlessSteel->AddElement(Ni, fractionmass=0.09);
        
        /// ------------- Stainless Steel Chamber -------------
        
        G4double bottom_chamber_height = chamber_height_ - top_part_chamber_height_;
        G4double chamber_wall_thickness = chamber_outer_radius_ - chamber_inner_radius_;
        G4double theta_max_outer = 110*deg/*acos(side_wall_center_to_front_ / chamber_outer_radius_) + 60 * deg*/;

        G4cout << "Chamber outer radius is: " << chamber_outer_radius_/mm << G4endl;
        G4cout << "Entrance window thickness is: " << entrance_window_thickness_/um << G4endl;
        G4cout << "Exit window thickness is: " << exit_window_thickness_/um << G4endl;
        
        
        G4Tubs *chamber_cyl = new G4Tubs("chamber_cylinder", 0, chamber_outer_radius_, 0.5 * top_part_chamber_height_ + 0.5 * bottom_chamber_height, 0, 2 * pi);
        
        double entrance_tube_dh = 5 * cm;
        G4Tubs *entrance_tube = new G4Tubs("entrance_tube", 0, entrance_tube_outer_radius_, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);
        
        G4RotationMatrix *rot_entrance_window = new G4RotationMatrix;
        rot_entrance_window->rotateY(90 * deg);
        G4VSolid *full_chamber = new G4UnionSolid("chamber_cyl + entrance_tube", chamber_cyl, entrance_tube, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_));
        
        chamber_full_log_ = new G4LogicalVolume(full_chamber, StainlessSteel, "full_chamber.log");
        
        /// ------------- Vacuum --------------
        

        G4Tubs *chamber_cyl_vac = new G4Tubs("top_chamber_cylinder_vac", 0, chamber_inner_radius_, (top_part_chamber_height_ + bottom_chamber_height - 2 *chamber_wall_thickness) / 2, 0, 2 * pi);
        G4Tubs *side_cyl_vac = new G4Tubs("side_cylinder_vac", 0, chamber_outer_radius_, exit_window_height_ / 2, 180*deg - theta_max_outer, 2*theta_max_outer);
        G4Tubs *entrance_tube_vac = new G4Tubs("entrance_tube_vac", 0, entrance_tube_outer_radius_ - chamber_wall_thickness, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);

        G4VSolid *full_chamber_vac = new G4UnionSolid("full_chamber_vac + entrance window", chamber_cyl_vac, entrance_tube_vac, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_));

        full_chamber_vac = new G4UnionSolid("full_chamber_vac + entrance window", full_chamber_vac, side_cyl_vac, NULL, G4ThreeVector(0, 0, entrance_window_ypos_));
        
        chamber_full_vac_log_ = new G4LogicalVolume(full_chamber_vac, Vacuum, "full_chamber_vac.log");
        
        /// ------------ Kapton windows --------------
        
        G4Tubs *entrance_window = new G4Tubs("entrance_window", 0, entrance_tube_inner_radius_, entrance_window_thickness_ / 2, 0, 2 * pi);
        G4Tubs *exit_window = new G4Tubs("side_cylinder_vac", chamber_outer_radius_, chamber_outer_radius_ + exit_window_thickness_, exit_window_height_ / 2, 180*deg - 110*deg/*theta_max_outer*/, 2*110*deg/*theta_max_outer*/);

        chamber_entrance_window_log_ = new G4LogicalVolume(entrance_window, Kapton, "entrance_window.log");
        chamber_exit_window_log_ = new G4LogicalVolume(exit_window, Kapton, "exit_window.log");
        
        /// ------------ Window Reinforcement --------

        G4double post_angle_ = asin(3.995*cm/(chamber_outer_radius_ - .25*cm))*(180/pi)*deg;
        G4cout << "Post angle is: " << (post_angle_) / deg << G4endl;

        G4Tubs *window_post = new G4Tubs("post", 0, .25*cm, exit_window_height_ /2, 0, 2*pi);
        window_post_log_ = new G4LogicalVolume(window_post, StainlessSteel, "window_post.log");

        /// ----------- Scintillator bar -----------
       
        /*G4double scintillator_bar_thickness_ = .5*cm;
        G4Tubs *scintillator_bar_A = new G4Tubs("scintillator_bars", chamber_outer_radius_ + exit_window_thickness_, chamber_outer_radius_ + exit_window_thickness_ + scintillator_bar_thickness_, exit_window_height_/2, 163.3*deg, 4*deg);
        G4Tubs *scintillator_bar_B = new G4Tubs("scintillator_bars", chamber_outer_radius_ + exit_window_thickness_, chamber_outer_radius_ + exit_window_thickness_ + scintillator_bar_thickness_, exit_window_height_/2, 192.7*deg, 4*deg);

        scintillator_bar_A_log_ = new G4LogicalVolume(scintillator_bar_A, Sc_material, "scintillator_bar_A.log");
        scintillator_bar_B_log_ = new G4LogicalVolume(scintillator_bar_B, Sc_material, "scintillator_bar_B.log");
        */

        // ------ placements
        
        G4RotationMatrix *rot_chamber = new G4RotationMatrix;
        rot_chamber->rotateY(90 * deg);
        rot_chamber->rotateX(90 * deg);
        rot_chamber->rotateZ(180 * deg);
        
        new G4PVPlacement(rot_chamber, G4ThreeVector(0,-entrance_window_ypos_,0), chamber_full_log_, "chamber_full", mother_volume_, false, 0);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,0), chamber_full_vac_log_, "chamber_full_vac", chamber_full_log_, false, 1);
        new G4PVPlacement(NULL, G4ThreeVector(0, 0, entrance_window_zpos_), chamber_entrance_window_log_, "chamber_entrance_window", mother_volume_, false, 2);
        new G4PVPlacement(rot_chamber, G4ThreeVector(0,0,0), chamber_exit_window_log_, "chamber_exit_window", mother_volume_, false, 3);
        new G4PVPlacement(NULL, G4ThreeVector(-(chamber_outer_radius_ - .25*cm)*sin(pi/2 - post_angle_),-(chamber_outer_radius_ - .25*cm)*cos(pi/2 - post_angle_), + entrance_window_ypos_), window_post_log_, "window_post", chamber_full_vac_log_, false, 27);
        new G4PVPlacement(NULL, G4ThreeVector(-(chamber_outer_radius_ - .25*cm)*sin(pi/2 - post_angle_),(chamber_outer_radius_ - .25*cm)*cos(pi/2 - post_angle_), + entrance_window_ypos_), window_post_log_, "window_post", chamber_full_vac_log_, false, 28);
        //new G4PVPlacement(rot_chamber, G4ThreeVector(0, 0, 0), scintillator_bar_A_log_, "scintillator_bar_A", mother_volume_, false, 35);
        //new G4PVPlacement(rot_chamber, G4ThreeVector(0, 0, 0), scintillator_bar_B_log_, "scintillator_bar_B", mother_volume_, false, 36);


        // ------ vis
        
        G4Colour colour_chamber_window;
        G4Colour colour_chamber;
        G4Colour colour_chamber_vacuum;
        //G4Colour colour_scintillator;
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4Colour::GetColour("scattering_chamber", colour_chamber);
        G4Colour::GetColour("scattering_chamber_airgap", colour_chamber_vacuum);
        //G4Colour::GetColour("scwall", colour_scintillator);
        G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
        G4VisAttributes* scatVisAtt = new G4VisAttributes(colour_chamber);
        G4VisAttributes* scatGVisAtt = new G4VisAttributes(colour_chamber_vacuum);
        //G4VisAttributes* scatSVisAtt = new G4VisAttributes(colour_scintillator);
        
        scatGVisAtt->SetVisibility(true);
        scatVisAtt->SetVisibility(true);
        scatWVisAtt->SetVisibility(true);
        //scatSVisAtt->SetVisibility(true);
        
        chamber_full_log_->SetVisAttributes(scatVisAtt);
        chamber_full_vac_log_->SetVisAttributes(scatGVisAtt);
        chamber_exit_window_log_->SetVisAttributes(scatWVisAtt);
        chamber_entrance_window_log_->SetVisAttributes(scatWVisAtt);
        window_post_log_->SetVisAttributes(scatVisAtt);
        //scintillator_bar_A_log_->SetVisAttributes(scatVisAtt);
        //scintillator_bar_B_log_->SetVisAttributes(scatVisAtt);

        DetectorParts->AttachBiasingOperator(chamber_full_log_);
        DetectorParts->AttachBiasingOperator(chamber_exit_window_log_);
        DetectorParts->AttachBiasingOperator(chamber_entrance_window_log_);
        DetectorParts->AttachBiasingOperator(window_post_log_);
        //DetectorParts->AttachBiasingOperator(scintillator_bar_A_log_);
        //DetectorParts->AttachBiasingOperator(scintillator_bar_B_log_);
        
        chamber_full_log_->SetUserLimits(StepLimit_);
        mother_for_target_log_ = chamber_full_vac_log_;
        
        
    } else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeTrapezoid]) {

        //  Materials for chamber and vacuum inside chamber
        G4NistManager* man = G4NistManager::Instance();
        
        G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");
        G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
        G4Material* SiO2   = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        G4Material* Sc_material = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        G4Element* C  = man->FindOrBuildElement("C");
        G4Element* Si = man->FindOrBuildElement("Si");
        G4Element* Cr = man->FindOrBuildElement("Cr");
        G4Element* Mn = man->FindOrBuildElement("Mn");
        G4Element* Fe = man->FindOrBuildElement("Fe");
        G4Element* Ni = man->FindOrBuildElement("Ni");
        G4Element* H  = man ->FindOrBuildElement("H");
        
        //FR4 Glass + Epoxy
        G4int ncomp, natoms;
        G4double densityEpox = 1.2*g/cm3;
        G4double densityFR4 = 1.86*g/cm3;
        G4double frac;
        G4Material* Epoxy = new G4Material("Epoxy" , densityEpox, ncomp=2);
        Epoxy->AddElement(H, natoms=2);
        Epoxy->AddElement(C, natoms=2);
        G4Material* G10 = new G4Material("G10"  , densityFR4, ncomp=2);
        G10->AddMaterial(SiO2, frac=0.528);
        G10->AddMaterial(Epoxy, frac=0.472);
        
        //Stainless steel
        G4double density = 8.06*g/cm3;
        G4int ncomponents = 6;
        G4double fractionmass = 0.;
        G4Material* StainlessSteel = new G4Material("StainlessSteel", density, ncomponents);
        StainlessSteel->AddElement(C, fractionmass=0.001);
        StainlessSteel->AddElement(Si, fractionmass=0.007);
        StainlessSteel->AddElement(Cr, fractionmass=0.18);
        StainlessSteel->AddElement(Mn, fractionmass=0.01);
        StainlessSteel->AddElement(Fe, fractionmass=0.712);
        StainlessSteel->AddElement(Ni, fractionmass=0.09);

        /// ------------- Stainless Steel Chamber -------------
        
        const G4double inches = 2.54 * cm;
        
        G4double theta = 65 * deg;
        G4double thetaR = 90 * deg - theta;
       
        G4double sc_wall_thickness = 0.375 * inches;
        
        G4double sc_half_height = 26 * inches / 2.0;
        G4double sc_half_y = 14.254 * inches / 2.0 * cos(thetaR) + sc_wall_thickness;
        G4double sc_half_x_small = 4.094 * inches / 2.0;
        G4double sc_half_x_large = sc_half_x_small + 2 * sc_half_y / tan(theta);
        
        G4Trap *chamber_ss_trap = new G4Trap("chamber_ss_trap", sc_half_height, 0*deg, 0*deg, sc_half_y, sc_half_x_small, sc_half_x_large, 0*deg,sc_half_y, sc_half_x_small, sc_half_x_large, 0*deg);
        
        G4double entrance_window_height = 5*cm;
        G4double entrance_window_outer_radius = 6 * inches / 2.0;
        G4double entrance_window_inner_radius = 7.5 * cm / 2.0;
        G4double entrance_window_sc_position = 16 * inches;
        G4double entrance_window_thickness = 50 * um;
        
        G4double downstream_exit_window_w = 2.756 * inches;
        G4double downstream_exit_window_h = 14 * inches;
        G4double downstream_exit_window_r = 1.181 * inches;
        G4double downstream_exit_window_thickness = 200*um;
        
        G4double side_exit_window_w = 11.867 * inches;
        G4double side_exit_window_h = downstream_exit_window_h;
        G4double side_exit_window_r = 1.969 * inches;
        G4double side_exit_window_thickness = 300*um;

        G4double chamber_vac_for_target_r = 10.0 *cm;
        G4double top_ss_height = 5*cm;
        
        
        G4Tubs *entrance_ss_tube = new G4Tubs("entrance_ss_tube", 0, entrance_window_outer_radius, entrance_window_height / 2.0, 0, 2 * pi);
        
        G4RotationMatrix *rot_entrance_window = new G4RotationMatrix;
        rot_entrance_window->rotateX(90 * deg);
        G4VSolid *full_chamber = new G4UnionSolid("chamber_ss", chamber_ss_trap, entrance_ss_tube, rot_entrance_window, G4ThreeVector(0, sc_half_y + entrance_window_height / 2, entrance_window_sc_position - sc_half_height));

        G4double sc_z_shift = (8.386 * inches + sc_wall_thickness) / 2.0 - sc_half_y;

        G4Tubs *top_cyl = new G4Tubs("top_ss_tube", 0, chamber_vac_for_target_r + sc_wall_thickness, top_ss_height, 0, 2*pi);
        
        full_chamber = new G4UnionSolid("full_chamber_vac + entrance window", full_chamber, top_cyl, NULL, G4ThreeVector(0, -sc_z_shift, sc_half_height + top_ss_height/2.0));

        chamber_full_log_ = new G4LogicalVolume(full_chamber, StainlessSteel, "full_chamber.log");

        // ------------- Vacuum --------------
        
        G4double dx_s = sc_wall_thickness * (-1. * tan(theta/2));
        G4double dx_l = sc_wall_thickness * (-1. / tan(theta/2));

        G4Trap *chamber_vac_trap = new G4Trap("chamber_vac_trap", sc_half_height - sc_wall_thickness, 0*deg, 0*deg, sc_half_y - sc_wall_thickness, sc_half_x_small + dx_s, sc_half_x_large + dx_l, 0*deg, sc_half_y - sc_wall_thickness, sc_half_x_small + dx_s, sc_half_x_large + dx_l, 0*deg);
        G4Tubs *entrance_vac_tube = new G4Tubs("entrance_vac_tube", 0, entrance_window_inner_radius, (entrance_window_height + sc_wall_thickness) / 2.0, 0, 2 * pi);
        G4Tubs *top_vac_cyl = new G4Tubs("top_vac_tube", 0, chamber_vac_for_target_r, top_ss_height, 0, 2*pi);

        G4VSolid *full_chamber_vac = new G4UnionSolid("chamber_vac", chamber_vac_trap, entrance_vac_tube, rot_entrance_window, G4ThreeVector(0, sc_half_y + (entrance_window_height - sc_wall_thickness) / 2, entrance_window_sc_position - sc_half_height));
        full_chamber_vac = new G4UnionSolid("full_chamber_vac", full_chamber_vac, top_vac_cyl, NULL, G4ThreeVector(0, -sc_z_shift, sc_half_height - sc_wall_thickness + top_ss_height/2.0));
        
        // openings for exit windows
        
        G4VSolid *downstream_exit_window = NULL;
        
        MakeWindow("downstream_exit_window", downstream_exit_window, downstream_exit_window_w, downstream_exit_window_h, sc_wall_thickness, downstream_exit_window_r);
        if (downstream_exit_window) {
            full_chamber_vac = new G4UnionSolid("full_chamber_vac", full_chamber_vac, downstream_exit_window, rot_entrance_window, G4ThreeVector(0,-sc_half_y + sc_wall_thickness / 2.0, entrance_window_sc_position - sc_half_height));
        }
        
        G4VSolid *side_exit_window = NULL;
        G4double shift = sc_half_y / sin(theta) - side_exit_window_w/2. - sc_wall_thickness / 2 / tan(theta) - 0.504 * inches;
        G4double xwin = (sc_half_x_small + dx_s + sc_half_x_large + dx_l) / 2. + 0.5 * sc_wall_thickness / cos(thetaR) - shift * cos(theta);
        G4double ywin = - shift * sin(theta);
        G4double zwin = entrance_window_sc_position - sc_half_height;
        G4RotationMatrix *rot_side_window_1 = new G4RotationMatrix;
        rot_side_window_1->rotateX(90 * deg);
        rot_side_window_1->rotateY(65 * deg);
        G4RotationMatrix *rot_side_window_2 = new G4RotationMatrix;
        rot_side_window_2->rotateX(90 * deg);
        rot_side_window_2->rotateY(-65 * deg);
        
        MakeWindow("side_exit_window", side_exit_window, side_exit_window_w, side_exit_window_h, sc_wall_thickness, side_exit_window_r);
        if (side_exit_window) {
            full_chamber_vac = new G4UnionSolid("full_chamber_vac", full_chamber_vac, side_exit_window, rot_side_window_1, G4ThreeVector( xwin, ywin, zwin));
            full_chamber_vac = new G4UnionSolid("full_chamber_vac", full_chamber_vac, side_exit_window, rot_side_window_2, G4ThreeVector(-xwin, ywin, zwin));
        }
        
        chamber_full_vac_log_ = new G4LogicalVolume(full_chamber_vac, Vacuum, "full_chamber_vac.log");
        
        // volume for target placement (we should not have to do that in the future)

        G4Tubs *chamber_vac_for_target = new G4Tubs("chamber_vac_for_target", 0, chamber_vac_for_target_r, sc_half_height, 0, 2 * pi);
        G4LogicalVolume *chamber_vac_for_target_log = new G4LogicalVolume(chamber_vac_for_target, Vacuum, "chamber_vac_for_target.log");

        // ------------- Window foils --------------

        G4Tubs *entrance_window = new G4Tubs("entrance_window", 0, entrance_window_inner_radius, entrance_window_thickness/2., 0, 2 * pi);
        chamber_entrance_window_log_ = new G4LogicalVolume(entrance_window, Kapton, "entrance_window.log");
        
        chamber_side_exit_window_log_ = NULL;
        chamber_exit_window_log_ = NULL;
        MakeWindow("side_exit_window", chamber_side_exit_window_log_, Kapton, side_exit_window_w, side_exit_window_h, side_exit_window_thickness, side_exit_window_r);
        MakeWindow("downstream_exit_window", chamber_exit_window_log_, Kapton, downstream_exit_window_w, downstream_exit_window_h,downstream_exit_window_thickness, downstream_exit_window_r);
        
        // ------ placements
        
        G4RotationMatrix *rot_chamber = new G4RotationMatrix;
        rot_chamber->rotateX(90 * deg);
        
        G4double sc_y_shift = -(entrance_window_sc_position - sc_half_height);
        new G4PVPlacement(rot_chamber, G4ThreeVector(0, sc_y_shift, -sc_z_shift), chamber_full_log_, "chamber_full", mother_volume_, false, 0);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,0), chamber_full_vac_log_, "chamber_full_vac", chamber_full_log_, false, 1);
        new G4PVPlacement(rot_chamber, G4ThreeVector(0, sc_half_y + entrance_window_height - entrance_window_thickness / 2, entrance_window_sc_position - sc_half_height), chamber_entrance_window_log_, "entrance_window", chamber_full_vac_log_, false, 3);
        new G4PVPlacement(rot_entrance_window, G4ThreeVector(0,-sc_half_y + downstream_exit_window_thickness / 2.0 , entrance_window_sc_position - sc_half_height), chamber_exit_window_log_, "downstream_exit_window", chamber_full_vac_log_, false, 4);

        G4double wshift = (sc_wall_thickness - side_exit_window_thickness) / 2.0;
        new G4PVPlacement(rot_side_window_1, G4ThreeVector(xwin + wshift * cos(thetaR), ywin - wshift * sin(thetaR), zwin), chamber_side_exit_window_log_, "side_exit_window1", chamber_full_vac_log_, false, 5);
        new G4PVPlacement(rot_side_window_2, G4ThreeVector(-xwin - wshift * cos(thetaR), ywin - wshift * sin(thetaR), zwin), chamber_side_exit_window_log_, "side_exit_window2", chamber_full_vac_log_, false, 6);
        
        G4RotationMatrix *rot_for_target = new G4RotationMatrix;
        rot_for_target->rotateZ(-90 * deg);
        G4double zpos = 49 * mm; // this is a manually adjusted term to bring the target center into the beamline
        new G4PVPlacement(rot_for_target, G4ThreeVector(0,-sc_z_shift, zpos), chamber_vac_for_target_log, "chamber_vac_for_target", chamber_full_vac_log_, false, 2);

        // ------ vis
        
        G4Colour colour_chamber_window;
        G4Colour colour_chamber;
        G4Colour colour_chamber_vacuum;
        //G4Colour colour_scintillator;
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4Colour::GetColour("scattering_chamber", colour_chamber);
        G4Colour::GetColour("scattering_chamber_vacuum", colour_chamber_vacuum);
        //G4Colour::GetColour("scwall", colour_scintillator);
        G4VisAttributes* scatWindowVisAtt = new G4VisAttributes(colour_chamber_window);
        G4VisAttributes* scatSteelAtt = new G4VisAttributes(colour_chamber);
        G4VisAttributes* scatVacuumVisAtt = new G4VisAttributes(colour_chamber_vacuum);
        //G4VisAttributes* scatSVisAtt = new G4VisAttributes(colour_scintillator);
        
        scatVacuumVisAtt->SetVisibility(true);
        scatSteelAtt->SetVisibility(true);
        scatWindowVisAtt->SetVisibility(true);
        //scatSVisAtt->SetVisibility(true);

        chamber_vac_for_target_log->SetVisAttributes(G4VisAttributes::Invisible);

//        test->SetVisAttributes(scatWVisAtt);
        
        chamber_full_log_->SetVisAttributes(scatSteelAtt);
        chamber_full_vac_log_->SetVisAttributes(scatVacuumVisAtt);
        if (chamber_exit_window_log_) chamber_exit_window_log_->SetVisAttributes(scatWindowVisAtt);
        if (chamber_entrance_window_log_) chamber_entrance_window_log_->SetVisAttributes(scatWindowVisAtt);
        if (chamber_side_exit_window_log_) chamber_side_exit_window_log_->SetVisAttributes(scatWindowVisAtt);
//        window_post_log_->SetVisAttributes(scatVisAtt);
        //scintillator_bar_A_log_->SetVisAttributes(scatVisAtt);
        //scintillator_bar_B_log_->SetVisAttributes(scatVisAtt);
        
        DetectorParts->AttachBiasingOperator(chamber_full_log_);
//        DetectorParts->AttachBiasingOperator(chamber_exit_window_log_);
//        DetectorParts->AttachBiasingOperator(chamber_entrance_window_log_);
//        DetectorParts->AttachBiasingOperator(window_post_log_);
        //DetectorParts->AttachBiasingOperator(scintillator_bar_A_log_);
        //DetectorParts->AttachBiasingOperator(scintillator_bar_B_log_);
        
        chamber_full_log_->SetUserLimits(StepLimit_);
        mother_for_target_log_ = chamber_vac_for_target_log;
        
        
    } else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeUM]) {

        //  Materials for chamber and vacuum inside chamber
        //
        G4NistManager* man = G4NistManager::Instance();

        G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic"); 
        G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Al = man->FindOrBuildMaterial("G4_Al");
        G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");

        G4double density = 8.06*g/cm3;
        G4int ncomponents = 6;
        G4double fractionmass = 0.;
        G4Material* StainlessSteel = new G4Material("StainlessSteel", density, ncomponents);
        G4Element* C  = man->FindOrBuildElement("C");
        G4Element* Si = man->FindOrBuildElement("Si");
        G4Element* Cr = man->FindOrBuildElement("Cr");
        G4Element* Mn = man->FindOrBuildElement("Mn");
        G4Element* Fe = man->FindOrBuildElement("Fe");
        G4Element* Ni = man->FindOrBuildElement("Ni");
        StainlessSteel->AddElement(C, fractionmass=0.001);
        StainlessSteel->AddElement(Si, fractionmass=0.007);
        StainlessSteel->AddElement(Cr, fractionmass=0.18);
        StainlessSteel->AddElement(Mn, fractionmass=0.01);
        StainlessSteel->AddElement(Fe, fractionmass=0.712);
        StainlessSteel->AddElement(Ni, fractionmass=0.09);

        G4Material* alMylar = new G4Material("AluminizedMylar", 1.762 * g/cm3, 2);
        alMylar->AddMaterial(Mylar, 0.7217);
        alMylar->AddMaterial(Al, 0.2783);
        

        G4double bottom_chamber_height = chamber_height_ - top_part_chamber_height_;
        G4double chamber_wall_thickness = chamber_outer_radius_ - chamber_inner_radius_;
        G4double theta_max_outer = acos(side_wall_center_to_front_ / chamber_outer_radius_) + 60 * deg;
        
        
        G4cout << "Chamber height is " << .5*bottom_chamber_height + .5*top_part_chamber_height_ << G4endl;

        /// --------- Stainless Steel ------------
        G4Tubs *chamber_cyl = new G4Tubs("chamber_cylinder", 0, chamber_outer_radius_, 0.5 * top_part_chamber_height_ + 0.5 * bottom_chamber_height, 0, 2 * pi);
        
        double entrance_tube_dh = 5 * cm;
        G4Tubs *entrance_tube = new G4Tubs("entrance_tube", 0, entrance_tube_outer_radius_, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);
        
        G4RotationMatrix *rot_entrance_window = new G4RotationMatrix;
        rot_entrance_window->rotateY(90 * deg);
        G4VSolid *full_chamber = new G4UnionSolid("chamber_cyl + entrance_tube", chamber_cyl, entrance_tube, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_));
        
        chamber_full_log_ = new G4LogicalVolume(full_chamber, StainlessSteel, "full_chamber.log");

        /// ------------- Vacuum --------------
        G4Tubs *chamber_cyl_vac = new G4Tubs("top_chamber_cylinder_vac", 0, chamber_inner_radius_, (top_part_chamber_height_ + bottom_chamber_height - 2 *chamber_wall_thickness) / 2, 0, 2 * pi);
        G4Tubs *side_cyl_vac = new G4Tubs("side_cylinder_vac", 0, chamber_outer_radius_, exit_window_height_ / 2, 180*deg - theta_max_outer, 2*theta_max_outer);
        G4Tubs *entrance_tube_vac = new G4Tubs("entrance_tube_vac", 0, entrance_tube_inner_radius_, 0.5 * (entrance_tube_height_ + entrance_tube_dh), 0, 2 * pi);

        G4VSolid *full_chamber_vac = new G4UnionSolid("full_chamber_vac + entrance window", chamber_cyl_vac, entrance_tube_vac, rot_entrance_window, G4ThreeVector(chamber_outer_radius_ + (entrance_tube_height_ - entrance_tube_dh) / 2, 0, entrance_window_ypos_ - chamber_wall_thickness));

        full_chamber_vac = new G4UnionSolid("full_chamber_vac + entrance window", full_chamber_vac, side_cyl_vac);
        
        chamber_full_vac_log_ = new G4LogicalVolume(full_chamber_vac, Vacuum, "full_chamber_vac.log");
  
        /// ------------- Entrance Window ------------
        G4Tubs *entrance_window = new G4Tubs("entrance_window", 0, entrance_tube_inner_radius_, entrance_window_thickness_ / 2, 0, 2 * pi);
        chamber_entrance_window_log_ = new G4LogicalVolume(entrance_window, Kapton, "entrance_window.log");

        // -------------- Exit Window ----------------
        G4Tubs *exit_window = new G4Tubs("exit_window", chamber_outer_radius_, chamber_outer_radius_ + exit_window_thickness_, exit_window_height_ / 2, 180*deg - theta_max_outer, 2*theta_max_outer);
        chamber_exit_window_log_ = new G4LogicalVolume(exit_window, Kapton, "exit_window.log");
        

        // --------- Placement ----------
        G4RotationMatrix *rot_chamber = new G4RotationMatrix;
        rot_chamber->rotateY(90 * deg);
        rot_chamber->rotateX(90 * deg);
        rot_chamber->rotateZ(180 * deg);
        
        new G4PVPlacement(rot_chamber, G4ThreeVector(0,-entrance_window_ypos_,0), chamber_full_log_, "chamber_full", mother_volume_, false, 0);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,chamber_wall_thickness), chamber_full_vac_log_, "chamber_full_vac", chamber_full_log_, false, 1);
        new G4PVPlacement(NULL, G4ThreeVector(0, 0, entrance_window_zpos_), chamber_entrance_window_log_, "chamber_entrance_window", mother_volume_, false, 2);
        new G4PVPlacement(rot_chamber, G4ThreeVector(0,0,0), chamber_exit_window_log_, "chamber_exit_window", mother_volume_, false, 3);
        
        // ------ vis
        
        G4Colour colour_chamber_window;
        G4Colour colour_chamber;
        G4Colour colour_chamber_vacuum;
        G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
        G4Colour::GetColour("scattering_chamber", colour_chamber);
        G4Colour::GetColour("scattering_chamber_airgap", colour_chamber_vacuum);
        G4VisAttributes* scatWVisAtt = new G4VisAttributes(colour_chamber_window);
        G4VisAttributes* scatVisAtt = new G4VisAttributes(colour_chamber);
        G4VisAttributes* scatGVisAtt = new G4VisAttributes(colour_chamber_vacuum);
        
        
        scatGVisAtt->SetVisibility(true);
        scatVisAtt->SetVisibility(true);
        scatWVisAtt->SetVisibility(true);
        

        /*
        scatGVisAtt->SetVisibility(false);
        scatVisAtt->SetVisibility(false);
        scatWVisAtt->SetVisibility(false);
        */

        chamber_full_log_->SetVisAttributes(scatVisAtt);
        chamber_full_vac_log_->SetVisAttributes(scatGVisAtt);
        chamber_entrance_window_log_->SetVisAttributes(scatWVisAtt);
        chamber_exit_window_log_->SetVisAttributes(scatWVisAtt);

        DetectorParts->AttachBiasingOperator(chamber_full_log_);
        DetectorParts->AttachBiasingOperator(chamber_exit_window_log_);
        DetectorParts->AttachBiasingOperator(chamber_entrance_window_log_);
        
        chamber_full_log_->SetUserLimits(StepLimit_);
        mother_for_target_log_ = chamber_full_vac_log_;

        // take TypeJeru2 as an example
        // - be mindful about your copy IDs we have typically:
        //     0 for the chamber wall,
        //     1 for the chamber vacuum (nested)
        //     2 for the entrance window
        //     3 for the exit window
        //  (>= 10 for the target, see below)
        // - resuse whenever possible the G4LogicalVolume pointers that are
        //   already defined as member variables in the class; we later attach
        //   the TGT sensitive detector to all kinds of chamber + target
        //   logical volumes.  This would be automatically done for your version
        //   if you reuse those variables
        // - limite the step size for the (most outer) logical volume
        // - set the mother_for_target_log_ to your vacuum logical volume.  The
        //   target will be placed into that volume and usually not into the world
        //   volume.
        // - attach the biasing operator to the logical volumes of the scattering
        //   chamber
        // ...
    }
    
    // ===============================================================
    //
    //   Build target
    //
    // ===============================================================
    
    if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type0]) {
        
        target_dz_ = 1 * CLHEP::cm;
        target_dx_ = 0.1 * CLHEP::cm;
        target_dy_ = 0.1 * CLHEP::cm;
        
        int target_segments = 1;
        min_p_ = 0;
        min_theta_ = 0;
        
        double dl = target_dz_ / target_segments;
        
        G4Tubs* solidtarget = new G4Tubs("target", 0., target_dx_, dl/2, 0, 360*CLHEP::deg);
        
        G4double Z, a, density;
        density = 0.0708*g/cm3;
        a = 1.00794*g/mole;
        
        G4Material* target_mat = new G4Material("liquidH2", Z=1., a, density);
        
        G4LogicalVolume* target_log = new G4LogicalVolume(solidtarget, target_mat, "target");
        
        for (int i = 0; i < target_segments; i++) {
            double offset = 0*CLHEP::mm;
            double z = -target_dz_/2 + (i + 0.5) * dl + offset;
            if (i % 2 == 0) {
                new G4PVPlacement (0, G4ThreeVector(0, 0, z), target_log,
                                   "target", mother_volume_, false, 0);
            }
        }
        
        G4Colour target_col;
        G4Colour::GetColour("target", target_col);
        target_log->SetVisAttributes(new G4VisAttributes(target_col));
        
        target_log->SetUserLimits(new G4UserLimits(0.5*CLHEP::mm));
        
        target_log_ = target_log;
        
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type1] ||
             DetectorParts->build[g4PSIDetectorParts::kTarget_Type2] ||
             DetectorParts->build[g4PSIDetectorParts::kTarget_Type3] ||
             DetectorParts->build[g4PSIDetectorParts::kTarget_Type4] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type1] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type2] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type3] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type4]) {
        
        // --- material
        
        G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
        G4Material* target_material = NULL;
        if (g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kEmptyTarget_Type1] ||
            g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kEmptyTarget_Type2] ||
            g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kEmptyTarget_Type3] ||
            g4PSIDetectorParts::getInstance()->build[g4PSIDetectorParts::kEmptyTarget_Type4]) {
            std::cout << "empty\n";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
        } else {
            std::cout << "full\n";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
        };
        std::cout << "Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;
        
        
        // --- geometry
        
        G4EllipticalTube* target_si_geom = new G4EllipticalTube("target_si.geom",
                                                                target_dx_ + target_wall_thickness_ + target_flask_mylar_gap_ + target_mylar_thickness_,
                                                                target_dy_ + target_wall_thickness_ + target_flask_mylar_gap_ + target_mylar_thickness_,
                                                                (target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/2 +
                                                                target_flask_mylar_gap_ + target_mylar_thickness_);
        G4EllipticalTube* target_vacuum_geom = new G4EllipticalTube("target_vacuum.geom",
                                                                    target_dx_ + target_wall_thickness_ + target_flask_mylar_gap_,
                                                                    target_dy_ + target_wall_thickness_ + target_flask_mylar_gap_,
                                                                    (target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/2 +
                                                                    target_flask_mylar_gap_);
        G4EllipticalTube* target_wall_geom = new G4EllipticalTube("target_wall_geom",
                                                                  target_dx_ + target_wall_thickness_,
                                                                  target_dy_ + target_wall_thickness_,
                                                                  (target_dz_ + target_entrance_cap_thickness_ + target_exit_cap_thickness_)/ 2.0);
        G4EllipticalTube* target_cell_geom = new G4EllipticalTube("target_geom",
                                                                  target_dx_,
                                                                  target_dy_,
                                                                  target_dz_ / 2.0);
        double test_plane_dz = 1*mm;
        G4EllipticalTube* target_out_test_geom = new G4EllipticalTube("target_out_test",
                                                                      target_dx_,
                                                                      target_dy_,
                                                                      test_plane_dz / 2.0);
        
        G4LogicalVolume* target_log = new G4LogicalVolume(target_cell_geom, target_material, "target_log", 0,0,0);
        G4LogicalVolume* target_wall_log = new G4LogicalVolume(target_wall_geom, Kapton, "target_wall_log", 0,0,0);
        G4LogicalVolume* target_vacuum_log = new G4LogicalVolume(target_vacuum_geom, Vacuum, "target_vacuum_log", 0,0,0);
        G4LogicalVolume* target_out_test_log = new G4LogicalVolume(target_out_test_geom, Vacuum, "target_out_test_log", 0,0,0);
        G4LogicalVolume* target_si_log = new G4LogicalVolume(target_si_geom, Mylar, "target_si_log", 0,0,0);
        
        G4Colour target_color;
        G4Colour::GetColour("target", target_color);
        G4Colour target_si_color;
        G4Colour::GetColour("target_si", target_si_color);
        G4Colour test_color;
        G4Colour::GetColour("white", test_color);
        target_log->SetVisAttributes(new G4VisAttributes(target_color));
        target_wall_log->SetVisAttributes(new G4VisAttributes(target_color));
        target_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);
        target_si_log->SetVisAttributes(new G4VisAttributes(target_si_color));
        target_out_test_log->SetVisAttributes(new G4VisAttributes(test_color));
        
        G4RotationMatrix* Rot = NULL;
        G4double targetPos_x = 0.0*m;
        G4double targetPos_y = 0.0*m;
        G4double targetPos_z = 0.0*m;
        
        if (mother_for_target_log_ != mother_volume_) {
            // the mother volume is the scattering chamber
            if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type3] ||
                DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type3]) {
                Rot = NULL;
            } else {
                Rot = new G4RotationMatrix;
                Rot->rotateZ(-90*deg);
                Rot->rotateX(-90*deg);
            }
        } else {
            if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type3] ||
                DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_Type3]) {
                Rot = new G4RotationMatrix;
                Rot->rotateX(-90*deg);
            } else {
                Rot = NULL;
            }
        }
        
        new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                          target_log, "target", target_wall_log, false, 0);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                          target_wall_log, "target_wall", target_vacuum_log, false, 0);
        new G4PVPlacement(NULL, G4ThreeVector(0,0,0),
                          target_vacuum_log, "target_vacuum", target_si_log, false, 0);
//        new G4PVPlacement(NULL, G4ThreeVector(0,0,target_dz_/2 + target_exit_cap_thickness_/2. + test_plane_dz/2. + 1*mm),
//                          target_out_test_log, "target_out_test", target_vacuum_log, false, 0);
        new G4PVPlacement(Rot, G4ThreeVector(targetPos_x,targetPos_y,targetPos_z),
                          target_si_log, "target_si", mother_for_target_log_, false, 0);
        
        
        target_log->SetUserLimits(StepLimit_);
        
        target_log_ = target_log;
        
    }
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeJeru]) {
        
        G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* target_material = NULL;
        if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru]) {
            G4cout << "Full Target: Filling target with liquid hydrogen\n" << G4endl;
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
        } else {
            G4cout << "Empty Target: Filling target with hydryogen gas\n" << G4endl;
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
        };
        std::cout << "Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;
        
        
        // --- geometry ---
        
        G4double target_radius = 30*mm;
        G4double target_film_thickness = 75*um;  // Dany Horovitz e-mail, 03/13/16
        G4double base_inner_radius = 43.8*mm / 2;
        
        // target volume
        
        G4double target_inner[] = {0, 0, 0, 0, 0, 0};
        G4double target_outer[] = {base_inner_radius, base_inner_radius, target_radius, target_radius, base_inner_radius, base_inner_radius};
        G4double target_z[] = {-54*mm, -40*mm, -40*mm, 40*mm, 40*mm, 54*mm};
        G4Polycone* target_volume = new G4Polycone("target_cell_volume", 0, 2*pi, 6 , target_z, target_inner, target_outer);
        
        target_log_ = new G4LogicalVolume(target_volume, target_material, "target_cell_volume.log", 0, 0, 0);

        // base
        
        G4Tubs* base1 = new G4Tubs("target_cell_base_1", base_inner_radius, 60*mm/2, 5*mm/2, 0., twopi);
        G4Tubs* base2 = new G4Tubs("target_cell_base_2", base_inner_radius, 61*mm/2, 9*mm/2, 0., twopi);
        G4Tubs* base3 = new G4Tubs("target_cell_base_3",    0*mm/2, 61*mm/2, 3*mm/2, 0., twopi);
        G4UnionSolid* base = new G4UnionSolid("target_cell_base12", base1, base2, NULL, G4ThreeVector(0., 0., (5+9)*mm/2));
        base = new G4UnionSolid("target_cell_base_top", base, base3, NULL, G4ThreeVector(0., 0., (5+3)*mm/2 + 9*mm));
        
        G4SubtractionSolid* target_base_top = (G4SubtractionSolid*) base;
        G4SubtractionSolid* target_base_bottom = (G4SubtractionSolid*) base;
    
        // for now do not include holes in the base, assume solid base

//        G4Tubs* base_hole = new G4Tubs("target_base_hole", 0, 6.8*mm/2, 8.6*mm/2, 0., twopi);
//        for (G4int i = 0; i < 2; i++) {
//            G4double r = base_inner_radius + 8.6*mm/2;
//            G4RotationMatrix* Rot1 = new G4RotationMatrix;
//            Rot1->rotateY(90*deg);
//            Rot1->rotateX(90*deg + angle1);
//            G4RotationMatrix* Rot2 = new G4RotationMatrix;
//            Rot2->rotateY(90*deg);
//            Rot2->rotateX(90*deg + angle2);
//            target_base_top = new G4SubtractionSolid("target_cell_base.log", target_base_top, base_hole, Rot1, G4ThreeVector(r*sin(angle1),r*cos(angle1), 8.5*mm));
//            target_base_bottom = new G4SubtractionSolid("target_cell_base.log", target_base_bottom, base_hole, Rot2, G4ThreeVector(r*sin(angle2),r*cos(angle2),8.5*mm));
//            angle1 *= -1;
//            angle2 *= -1;
//        }
        
        target_base_top_log_ = new G4LogicalVolume(target_base_top, Kapton, "target_cell_base.log", 0, 0, 0);
        target_base_bottom_log_ = new G4LogicalVolume(target_base_bottom, Kapton, "target_cell_base.log", 0, 0, 0);

        // film
        
        G4Tubs* target_cell_film = new G4Tubs("target_cell_film", target_radius, target_radius + target_film_thickness, 40*mm, 0, 2*pi);
        target_cell_film_log_ = new G4LogicalVolume(target_cell_film, Kapton, "target_cell_film.log", 0, 0, 0);

        // nipple
        
        G4Tubs* nipple = new G4Tubs("target_nipple", 0., 25*mm/2, 20*mm/2, 0, 2*pi);
        target_nipple_log_ = new G4LogicalVolume(nipple, Kapton, "target_nipple.log", 0, 0, 0);
        
        // pipes
        
        G4double pipe_length_h = 70.5*mm - 61*mm/2 - 25*mm/2;
        G4Tubs* pipe_h = new G4Tubs("target_pipe_h", 0., 10*mm/2, pipe_length_h/2, 0, 2*pi);
        target_pipe_h_log_ = new G4LogicalVolume(pipe_h, Kapton, "target_pipe.log", 0, 0, 0);
        G4double pipe_length_v1 = 104*mm;
        G4Tubs* pipe_v1 = new G4Tubs("target_pipe_v1", 0., 10*mm/2, pipe_length_v1/2, 0, 2*pi);
        target_pipe_v_log_ = new G4LogicalVolume(pipe_v1, Kapton, "target_pipe.log", 0, 0, 0);
        
        // pressure probe
        
        const G4double probe_p_length1 = 8.5 * mm;
        const G4double probe_p_length2 = 20 * mm;
        const G4double probe_p_r1 = 10 * mm / 2;
        const G4double probe_p_r2 = 15 * mm / 2;
        
        const G4double probe_p_inner[] = {0, 0, 0, 0};
        const G4double probe_p_outer[] = {probe_p_r1, probe_p_r1, probe_p_r2, probe_p_r2};
        const G4double probe_p_z[] = {0, probe_p_length1, probe_p_length1, probe_p_length2};
        G4Polycone* probe_p_volume = new G4Polycone("probe_p_cell_volume", 0, 2*pi, 4, probe_p_z, probe_p_inner, probe_p_outer);
        
        target_probe_p_log_ = new G4LogicalVolume(probe_p_volume, Kapton, "probe_p_cell_volume.log", 0, 0, 0);
        
        // temperature probe
        
        const G4double probe_t_length1 = 8.5 * mm;
        const G4double probe_t_length2 = 30 * mm;
        const G4double probe_t_r1 = 10 * mm / 2;
        const G4double probe_t_r2 = 15 * mm / 2;
        
        const G4double probe_t_inner[] = {0, 0, 0, 0};
        const G4double probe_t_outer[] = {probe_t_r1, probe_t_r1, probe_t_r2, probe_t_r2};
        const G4double probe_t_z[] = {0, probe_t_length1, probe_t_length1, probe_t_length2};
        G4Polycone* probe_t_volume = new G4Polycone("probe_t_cell_volume", 0, 2*pi, 4, probe_t_z, probe_t_inner, probe_t_outer);
        
        target_probe_t_log_ = new G4LogicalVolume(probe_t_volume, Kapton, "probe_t_cell_volume.log", 0, 0, 0);
        
        // placements
        // ----------
        
        const G4int NCELL = 2;
        G4double angle1 = 55*deg + 90*deg;
        G4double angle2 = 35*deg - 90*deg;
        G4RotationMatrix* Rot = new G4RotationMatrix;
        Rot->rotateX(180*deg);
        Rot->rotateZ(180*deg);
        G4RotationMatrix* Rot1A = new G4RotationMatrix;
        Rot1A->rotateY(90*deg);
        Rot1A->rotateX(90*deg + angle1);
        G4RotationMatrix* Rot1B = new G4RotationMatrix;
        Rot1B->rotateY(90*deg);
        Rot1B->rotateX(-90*deg - angle1);
        G4RotationMatrix* Rot2A = new G4RotationMatrix;
        Rot2A->rotateY(90*deg);
        Rot2A->rotateX(90*deg - angle2);
        G4RotationMatrix* Rot2B = new G4RotationMatrix;
        Rot2B->rotateY(90*deg);
        Rot2B->rotateX(-90*deg + angle2);
        for (G4int icell = 0; icell < NCELL; icell++) {
            G4double dz = icell * 124 * mm + 27.25 * mm;
            G4int id = (icell + 1) * 10; // target cell copy IDs are 10 and 20
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_log_,
                              Form("target_cell%d_volume", icell), mother_for_target_log_, false, id);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_cell_film_log_,
                              Form("target_cell%d_film", icell), mother_for_target_log_, false, id+1);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz+40*mm+5*mm/2), target_base_top_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);
            new G4PVPlacement(Rot, G4ThreeVector(0,0, dz-40*mm-5*mm/2), target_base_bottom_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);
            
            // place nipples
            
            G4double r = 70.5 * mm;

            new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm), target_nipple_log_, Form("target_nipple%d", 2*icell), mother_for_target_log_, false, id+3);
            new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm), target_nipple_log_, Form("target_nipple%d", 2*icell+1), mother_for_target_log_, false, id+3);
            
            // place pipes
            
            new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), mother_for_target_log_, false, id+4);
            new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), mother_for_target_log_, false, id+4);

            r = 61*mm/2 + pipe_length_h / 2;
            new G4PVPlacement(Rot1A, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+11*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell), mother_for_target_log_, false, id+5);
            new G4PVPlacement(Rot2A, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-11*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell+1), mother_for_target_log_, false, id+5);
            
            // place probes
            
            r = 61*mm/2;
            new G4PVPlacement(Rot1B, G4ThreeVector(r*sin(angle1), -r*cos(angle1), dz+40*mm+11*mm), target_probe_t_log_, Form("target_probe_t%d", 2*icell), mother_for_target_log_, false, id+6);
            new G4PVPlacement(Rot2B, G4ThreeVector(-r*sin(angle2), -r*cos(angle2), dz-40*mm-11*mm), target_probe_p_log_, Form("target_probe_p%d", 2*icell), mother_for_target_log_, false, id+6);

        };
        
        //
        G4Colour target_color;
        G4Colour::GetColour("target", target_color);
        G4Colour target_si_color;
        G4Colour::GetColour("target_si", target_si_color);
        G4Colour target_base_color;
        G4Colour::GetColour("target_base", target_base_color);
        G4Colour target_film_color;
        G4Colour::GetColour("target_film", target_film_color);
        G4Colour target_pipe_color;
        G4Colour::GetColour("target_pipe", target_pipe_color);
        target_log_->SetVisAttributes(new G4VisAttributes(target_color));
        target_base_top_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_base_bottom_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_cell_film_log_->SetVisAttributes(new G4VisAttributes(target_film_color));
        target_nipple_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_pipe_h_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_pipe_v_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_probe_p_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_probe_t_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        
        target_log_->SetUserLimits(StepLimit_);
        
    }
    
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeUM]) {
        
        G4NistManager* man = G4NistManager::Instance();

        G4String target_setting;

        G4Material* target_material = NULL;
        if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM]) {
            target_setting = "Full Target";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
        } else {
            target_setting = "Empty Target";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
        };

        G4cout << "Target setting: " << target_setting << "; Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << G4endl;

        G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Al = man->FindOrBuildMaterial("G4_Al");
        G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");

        G4double density = 8.06*g/cm3;
        G4int ncomponents = 6;
        G4double fractionmass = 0.;
        G4Material* StainlessSteel = new G4Material("StainlessSteel", density, ncomponents);
        G4Element* C  = man->FindOrBuildElement("C");
        G4Element* Si = man->FindOrBuildElement("Si");
        G4Element* Cr = man->FindOrBuildElement("Cr");
        G4Element* Mn = man->FindOrBuildElement("Mn");
        G4Element* Fe = man->FindOrBuildElement("Fe");
        G4Element* Ni = man->FindOrBuildElement("Ni");
        StainlessSteel->AddElement(C, fractionmass=0.001);
        StainlessSteel->AddElement(Si, fractionmass=0.007);
        StainlessSteel->AddElement(Cr, fractionmass=0.18);
        StainlessSteel->AddElement(Mn, fractionmass=0.01);
        StainlessSteel->AddElement(Fe, fractionmass=0.712);
        StainlessSteel->AddElement(Ni, fractionmass=0.09);

        G4Material* alMylar = new G4Material("AluminizedMylar", 1.762 * g/cm3, 2);
        alMylar->AddMaterial(Mylar, 0.7217);
        alMylar->AddMaterial(Al, 0.2783);

        //Target cell 
        G4double cell_distance = 13.4546 *cm;
        G4double capInnerR = 0;
        G4double capOuterR = 3.0 *cm;
        G4double capHeightT = .49085 *cm;
        G4double capHeightB = .5296 *cm;
        G4double kaptonInnerR = 2.9880 *cm;
        G4double kaptonOuterR = 3 *cm;
        G4double cellHeight = 9.0221/2 *cm;
        G4double lH2InnerR = 0 *cm;
        G4double lH2OuterR = 2.988 *cm;

        // tube/frame 
        G4double tubeInnerR = (.5 - .01) *cm;
        G4double tubeOutterR = .5 *cm;
        G4double outtertubeheight = 54.88525/2 *cm;
        G4double uppertubelength = 7.2244/2 *cm; 
        G4double connectortubelength = 1.3135/2 *cm;
        G4double cuttingcube_side = pow(2, .5) * tubeOutterR; // sqrt(2) * diameter of tube/2
        G4double connect_cap_dist = 1.2756 *cm;
        G4double frame_distance = 7 *cm;
        G4double frame_angle = 25 *deg; 
        G4double pipe_angle = 30.27 *deg; 
        G4double bottom_frame_thickness = .5312/2 *cm;
        G4double bottom_frame_length = 17.6326/4 *cm; 
        G4double bottom_frame_distance = .2586 *cm;
        G4double bottom_frame_width = 1.5/2 *cm;
        G4double frame_to_cell_distance = 20.59695 *cm;

        // Alumized mylar layers
        G4double mylarheight = 16 *cm;
        G4double mylarthickness = .001 *cm;
        G4double mylar1R = 9 *cm;
        G4double mylar2R = 10 *cm;
        G4double mylar3R = 11 *cm;

        //Build target
        G4VSolid* target_base_top = new G4Tubs("topcap", capInnerR, capOuterR, capHeightT, 0, twopi);
        G4VSolid* target_base_bottom = new G4Tubs("bottomcap", capInnerR, capOuterR, capHeightB, 0, twopi);
        G4VSolid* target_cell_film = new G4Tubs("kaptoncell", kaptonInnerR, kaptonOuterR, cellHeight, 0, twopi);
        G4VSolid* target_volume = new G4Tubs("hydrogen", lH2InnerR, lH2OuterR, cellHeight, 0, twopi);

        //Build heat shields
        G4VSolid* mylar1 = new G4Tubs("mylar1", mylar1R - mylarthickness, mylar1R, mylarheight, 0, twopi);
        G4VSolid* mylar2 = new G4Tubs("mylar2", mylar2R - mylarthickness, mylar2R, mylarheight, 0, twopi);
        G4VSolid* mylar3 = new G4Tubs("mylar3", mylar3R - mylarthickness, mylar3R, mylarheight, 0, twopi);


         // Shape for tube/frame
        G4VSolid* frame = new G4Tubs("frame", tubeInnerR, tubeOutterR, outtertubeheight - 4 *cm, 0, twopi);
        G4VSolid* uppertube = new G4Tubs("uppertube", tubeInnerR, tubeOutterR, uppertubelength, 0, twopi);
        G4VSolid* connectortube = new G4Tubs("connectortube", tubeInnerR, tubeOutterR, connectortubelength, 0, twopi);
        G4VSolid* cuttingcube = new G4Box("cuttingcube", cuttingcube_side, cuttingcube_side, cuttingcube_side);
        G4RotationMatrix* tRot = new G4RotationMatrix();
        tRot->rotateZ(0);
        tRot->rotateX(-45*deg);
        tRot->rotateY(0);
        G4VSolid* pipesub = new G4SubtractionSolid("uppertube - cuttingcube", uppertube, cuttingcube, tRot, G4ThreeVector(0, -tubeOutterR , uppertubelength));
        G4VSolid* connectsub = new G4SubtractionSolid("connectortube - cuttingcube", connectortube, cuttingcube, tRot, G4ThreeVector(0, -tubeOutterR, connectortubelength));

        G4RotationMatrix* fRot = new G4RotationMatrix();
        fRot->rotateZ(130*deg);
        G4VSolid* bottomframepart = new G4Box("bottomframepart", bottom_frame_length, bottom_frame_width, bottom_frame_thickness);
        G4VSolid* bottomframe = new G4UnionSolid("bottomframe", bottomframepart, bottomframepart,fRot, G4ThreeVector(-(bottom_frame_length)*(1 + sin(40*deg)), -(bottom_frame_length)*cos(40*deg), 0));

        target_log_ = new G4LogicalVolume(target_volume, target_material, "target_cell_volume.log", 0, 0, 0);
        target_base_top_log_ = new G4LogicalVolume(target_base_top, StainlessSteel, "target_cell_base.log", 0, 0, 0);
        target_base_bottom_log_ = new G4LogicalVolume(target_base_bottom, StainlessSteel, "target_cell_base.log", 0, 0, 0);
        target_cell_film_log_ = new G4LogicalVolume(target_cell_film, Kapton, "target_cell_film.log");

        // build logical volume for frame
        rightframe_log_ = new G4LogicalVolume(frame, StainlessSteel, "rightframe.log");
        leftframe_log_ = new G4LogicalVolume(frame, StainlessSteel, "leftframe.log");
        Tpipe_log_ = new G4LogicalVolume(pipesub, StainlessSteel, "Tpipe.log");
        Tconnect_log_ = new G4LogicalVolume(connectsub, StainlessSteel, "Tconnect.log");
        Bpipe_log_ = new G4LogicalVolume(pipesub, StainlessSteel, "Bpipe.log");
        Bconnect_log_ = new G4LogicalVolume(connectsub, StainlessSteel, "Bconnect.log");

        // build logical volume for mylar
        mylar1_log_ = new G4LogicalVolume(mylar1, alMylar, "mylar1Log");
        mylar2_log_ = new G4LogicalVolume(mylar2, alMylar, "mylar2Log");
        mylar3_log_ = new G4LogicalVolume(mylar3, alMylar, "mylar3Log");


        const G4int NCELL = 2;
        G4double angle1 = 55*deg + 90*deg;
        G4double angle2 = 35*deg - 90*deg;
        G4RotationMatrix* pipRot = new G4RotationMatrix();
        pipRot->rotateX(-90*deg);
        pipRot->rotateY(180*deg - pipe_angle);
        pipRot->rotateZ(0);
        G4RotationMatrix* conRot = new G4RotationMatrix();
        conRot->rotateX(0);
        conRot->rotateY(0);
        conRot->rotateZ(-pipe_angle);
        G4RotationMatrix* pipRot2 = new G4RotationMatrix();
        pipRot2->rotateX(-90*deg);
        pipRot2->rotateY(+pipe_angle );
        pipRot2->rotateZ(180*deg);
        G4RotationMatrix* conRot2 = new G4RotationMatrix();
        conRot2->rotateX(0);
        conRot2->rotateY(180*deg);
        conRot2->rotateZ(180*deg - pipe_angle);

        for (G4int icell = 0; icell < NCELL; icell++) {
            G4double dz = icell * -cell_distance + 27.25 * mm;
            G4int id = (icell + 1) * 10; // target cell copy IDs are 10 and 20
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_log_,
                              Form("target_cell%d_volume", icell), mother_for_target_log_, false, id);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_cell_film_log_,
                              Form("target_cell%d_film", icell), mother_for_target_log_, false, id+1);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz+cellHeight + capHeightT), target_base_top_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);
            new G4PVPlacement(NULL, G4ThreeVector(0,0, dz-cellHeight - capHeightT), target_base_bottom_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);

            //Place pipes and connectors
            new G4PVPlacement(pipRot, G4ThreeVector( (uppertubelength - tubeOutterR)*sin(pipe_angle), -(uppertubelength  - tubeOutterR + connect_cap_dist  - (uppertubelength - tubeOutterR)*(1-cos(pipe_angle))),(cellHeight + 2*capHeightT + 2*connectortubelength - tubeOutterR) + dz), Tpipe_log_, 
                              Form("Tpipe_%d", icell),mother_for_target_log_, false, id+5);
            new G4PVPlacement(conRot, G4ThreeVector(0, -connect_cap_dist, cellHeight + 2*capHeightT + connectortubelength + dz), Tconnect_log_, 
                              Form("Tconnect_%d", icell), mother_for_target_log_, false, id+6);
            new G4PVPlacement(pipRot2, G4ThreeVector( (uppertubelength - tubeOutterR)*sin(pipe_angle), +(uppertubelength  - tubeOutterR + connect_cap_dist  - (uppertubelength - tubeOutterR)*(1-cos(pipe_angle))), (-cellHeight - 2*capHeightT  - 2*connectortubelength + tubeOutterR) + dz), Bpipe_log_, 
                              Form("Bpipe_%d", icell),mother_for_target_log_, false, id+7);
            new G4PVPlacement(conRot2, G4ThreeVector(0, connect_cap_dist, (-cellHeight - 2*capHeightT - connectortubelength + dz)), Bconnect_log_, 
                              Form("Bconnect_%d", icell), mother_for_target_log_, false, id+8);
        

        }

        // Now we must place the frame 
        new G4PVPlacement(NULL, G4ThreeVector(frame_distance*sin(frame_angle), frame_distance*cos(frame_angle), outtertubeheight - frame_to_cell_distance + 27.25 *mm - 2.5*cm), rightframe_log_, "frame", mother_for_target_log_, false, 13);
        new G4PVPlacement(NULL, G4ThreeVector(frame_distance*sin(frame_angle), -frame_distance*cos(frame_angle), outtertubeheight - frame_to_cell_distance + 27.25 *mm - 2.5*cm), leftframe_log_, "frame", mother_for_target_log_, false, 14);


        G4Colour target_color;
        G4Colour::GetColour("target", target_color);
        G4Colour target_si_color;
        G4Colour::GetColour("target_si", target_si_color);
        G4Colour target_base_color;
        G4Colour::GetColour("target_base", target_base_color);
        G4Colour target_film_color;
        G4Colour::GetColour("target_film", target_film_color);
        G4Colour target_pipe_color;
        G4Colour::GetColour("target_pipe", target_pipe_color);
        target_log_->SetVisAttributes(new G4VisAttributes(target_color));
        target_base_top_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_base_bottom_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_cell_film_log_->SetVisAttributes(new G4VisAttributes(target_film_color));
        Bpipe_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        Tpipe_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        Bconnect_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        Tconnect_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        rightframe_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        leftframe_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));


        target_log_->SetUserLimits(StepLimit_);


        // take the TypeJeru target as example
        // - be mindful about the copy IDs
        //     10 : cell with 'target_material' (LH2 or H2 gas)
        //     11 : cell walls
        //     12 : ...
        //     ...
        //     20 : for 2nd cell
        //     21 : ...
        // - make sureto communicate your logical volumes via
        //      g4PSIDetectorParts::getInstance()->SetTargetVolume( ... )
        //   this is needed to record 'true' scattering events in the Stepping Action
        
    }

    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM2] ||
             DetectorParts->build[g4PSIDetectorParts::kEmptyTarget_TypeUM2]) {

        
        const G4double base_inner_radius = base_oneOutterR - base_one_thickness; // measurements from Priay 2/24/17 Old: 43.8*mm / 2;
        

        G4cout << "Pipe angle = " << pipe_angle/deg << " deg" << G4endl;
        G4cout << "Pipe distance = " << pipe_dist/mm << " mm" << G4endl;
        G4cout << "Target radius = " << radius/mm << " mm" << G4endl;
        G4cout << "Kapton thickness = " << target_film_thickness/um << " um" << G4endl;
  
        G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
        G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
        G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* Copper = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
        G4Material* AlMylar = new G4Material("AluminizedMylar", 1.762 * g/cm3, 2);
        AlMylar->AddMaterial(Mylar, 0.7217);
        AlMylar->AddMaterial(Al, 0.2783);

        G4Material* target_material = NULL;
        if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM2]) {
            std::cout << "Full Target\n";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
        } else {
            std::cout << "Empty Target\n";
            target_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
            //target_material = Vacuum;
        };
        std::cout << "Target density: " << target_material->GetDensity() / (g/cm3) << " g/cm^3" << std::endl;

        
        // target volume
        
        G4double target_inner[] = {0, 0, 0, 0, 0, 0};
        G4double target_outer[] = {base_inner_radius, base_inner_radius, radius, radius, base_inner_radius, base_inner_radius};
        G4double target_z[] = {-cellHeight, -cellHeight + capHeight, -cellHeight + capHeight, cellHeight - capHeight, cellHeight - capHeight, cellHeight};
        G4Polycone* target_volume = new G4Polycone("target_cell_volume", 0, 2*pi, 6 , target_z, target_inner, target_outer);
        
        target_log_ = new G4LogicalVolume(target_volume, target_material, "target_cell_volume.log", 0, 0, 0);

        // base
        
        G4Tubs* base1 = new G4Tubs("target_cell_base_1", base_inner_radius, base_oneOutterR, capHeight, 0., twopi);
        G4Tubs* base2 = new G4Tubs("target_cell_base_3",    0*mm/2, base_twoOutterR, capHeight, 0., twopi);
        G4UnionSolid* base = new G4UnionSolid("target_cell_base_top", base1, base2, NULL, G4ThreeVector(0., 0., 2*capHeight));
        
        G4SubtractionSolid* target_base_top = (G4SubtractionSolid*) base;
        G4SubtractionSolid* target_base_bottom = (G4SubtractionSolid*) base;
    
        // for now do not include holes in the base, assume solid base

//        G4Tubs* base_hole = new G4Tubs("target_base_hole", 0, 6.8*mm/2, 8.6*mm/2, 0., twopi);
//        for (G4int i = 0; i < 2; i++) {
//            G4double r = base_inner_radius + 8.6*mm/2;
//            G4RotationMatrix* Rot1 = new G4RotationMatrix;
//            Rot1->rotateY(90*deg);
//            Rot1->rotateX(90*deg + angle1);
//            G4RotationMatrix* Rot2 = new G4RotationMatrix;
//            Rot2->rotateY(90*deg);
//            Rot2->rotateX(90*deg + angle2);
//            target_base_top = new G4SubtractionSolid("target_cell_base.log", target_base_top, base_hole, Rot1, G4ThreeVector(r*sin(angle1),r*cos(angle1), 8.5*mm));
//            target_base_bottom = new G4SubtractionSolid("target_cell_base.log", target_base_bottom, base_hole, Rot2, G4ThreeVector(r*sin(angle2),r*cos(angle2),8.5*mm));
//            angle1 *= -1;
//            angle2 *= -1;
//        }
        
        target_base_top_log_ = new G4LogicalVolume(target_base_top, Copper, "target_cell_base.log", 0, 0, 0);
        target_base_bottom_log_ = new G4LogicalVolume(target_base_bottom, Copper, "target_cell_base.log", 0, 0, 0);

        // film
        
        G4Tubs* target_cell_film = new G4Tubs("target_cell_film", base_oneOutterR, base_oneOutterR + target_film_thickness, cellHeight, 0, 2*pi);
        target_cell_film_log_ = new G4LogicalVolume(target_cell_film, Kapton, "target_cell_film.log", 0, 0, 0);

        // nipple
        
        G4Tubs* nipple = new G4Tubs("target_nipple", 0., 7*mm/2, 20*mm/2, 0, 2*pi);// raidus was 25/2*mm
        target_nipple_log_ = new G4LogicalVolume(nipple, Copper, "target_nipple.log", 0, 0, 0);
        

        // pipes
        
        G4double pipe_length_h = 60*mm - 61*mm/2 - 7*mm/2; // first was 70.5 , second was 61/2*mm, third was 25/2*mm
        G4Tubs* pipe_h = new G4Tubs("target_pipe_h", 5.35*mm/2, 6.35*mm/2, pipe_length_h/2, 0, 2*pi);
        target_pipe_h_log_ = new G4LogicalVolume(pipe_h, Copper, "target_pipe.log", 0, 0, 0);
        G4double pipe_length_v1 = cell_distance; // was 104
        G4Tubs* pipe_v1 = new G4Tubs("target_pipe_v1", 5.35*mm/2, 6.35*mm/2, pipe_length_v1/2, 0, 2*pi);
        target_pipe_v_log_ = new G4LogicalVolume(pipe_v1, Copper, "target_pipe.log", 0, 0, 0);

        // aluminzed mylar superinsulation
        G4double insulation_thickness = .03125*mm;              //~5 layers of superinsulation
        G4Trd *Trd1 = new G4Trd("Trd", 2*cm, 7*cm, 8*cm, 8*cm, 3.4*cm);
        G4Trd *Trd2 = new G4Trd("Trd", 2*cm - insulation_thickness, 7*cm - insulation_thickness, 9*cm, 9*cm, 3.4*cm - insulation_thickness);
        G4SubtractionSolid* Trd = new G4SubtractionSolid("Trd", Trd1, Trd2, NULL, G4ThreeVector(0., 0., 0.));
        
        target_si_log_ = new G4LogicalVolume(Trd, AlMylar, "target_si", 0,0,0);

        
        // pressure probe
        /*
        const G4double probe_p_length1 = 8.5 * mm;
        const G4double probe_p_length2 = 20 * mm;
        const G4double probe_p_r1 = 10 * mm / 2;
        const G4double probe_p_r2 = 15 * mm / 2;
        
        const G4double probe_p_inner[] = {0, 0, 0, 0};
        const G4double probe_p_outer[] = {probe_p_r1, probe_p_r1, probe_p_r2, probe_p_r2};
        const G4double probe_p_z[] = {0, probe_p_length1, probe_p_length1, probe_p_length2};
        G4Polycone* probe_p_volume = new G4Polycone("probe_p_cell_volume", 0, 2*pi, 4, probe_p_z, probe_p_inner, probe_p_outer);
        
        target_probe_p_log_ = new G4LogicalVolume(probe_p_volume, Kapton, "probe_p_cell_volume.log", 0, 0, 0);
        
        // temperature probe
        
        const G4double probe_t_length1 = 8.5 * mm;
        const G4double probe_t_length2 = 30 * mm;
        const G4double probe_t_r1 = 10 * mm / 2;
        const G4double probe_t_r2 = 15 * mm / 2;
        
        const G4double probe_t_inner[] = {0, 0, 0, 0};
        const G4double probe_t_outer[] = {probe_t_r1, probe_t_r1, probe_t_r2, probe_t_r2};
        const G4double probe_t_z[] = {0, probe_t_length1, probe_t_length1, probe_t_length2};
        G4Polycone* probe_t_volume = new G4Polycone("probe_t_cell_volume", 0, 2*pi, 4, probe_t_z, probe_t_inner, probe_t_outer);
        
        target_probe_t_log_ = new G4LogicalVolume(probe_t_volume, Kapton, "probe_t_cell_volume.log", 0, 0, 0);
        */
        // placements
        // ----------

        const G4double target_offset = entrance_window_ypos_; //was 27.25*mm 
        
        const G4int NCELL = 2;
        G4double angle1 = 180*deg - pipe_angle;
        G4double angle2 = -pipe_angle;
        G4RotationMatrix* Rot = new G4RotationMatrix;
        Rot->rotateX(180*deg);
        Rot->rotateZ(180*deg);
        G4RotationMatrix* Rot1A = new G4RotationMatrix;
        Rot1A->rotateY(90*deg);
        Rot1A->rotateX(90*deg + angle1);
        G4RotationMatrix* Rot1B = new G4RotationMatrix;
        Rot1B->rotateY(90*deg);
        Rot1B->rotateX(-90*deg - angle1);
        G4RotationMatrix* Rot2A = new G4RotationMatrix;
        Rot2A->rotateY(90*deg);
        Rot2A->rotateX(90*deg - angle2);
        G4RotationMatrix* Rot2B = new G4RotationMatrix;
        Rot2B->rotateY(90*deg);
        Rot2B->rotateX(-90*deg + angle2);
        for (G4int icell = 0; icell < NCELL; icell++) {
            G4double dz = icell * -cell_distance + target_offset;
            G4int id = (icell + 1) * 10; // target cell copy IDs are 10 and 20
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_log_,
                              Form("target_cell%d_volume", icell), mother_for_target_log_, false, id);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz), target_cell_film_log_,
                              Form("target_cell%d_film", icell), mother_for_target_log_, false, id+1);
            new G4PVPlacement(NULL, G4ThreeVector(0,0,dz + cellHeight - 10/2*mm), target_base_top_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);
            new G4PVPlacement(Rot, G4ThreeVector(0,0, dz- cellHeight + 10/2*mm), target_base_bottom_log_,
                              Form("target_cell%d_base", icell), mother_for_target_log_, false, id+2);
            
            G4double r = pipe_dist;
         
            // place nipples
            
            new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+14.5*mm), target_nipple_log_, Form("target_nipple%d", 2*icell), mother_for_target_log_, false, id+3);
            new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-14.5*mm), target_nipple_log_, Form("target_nipple%d", 2*icell+1), mother_for_target_log_, false, id+3);
            
            // place pipes
            
            new G4PVPlacement(NULL, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+14.5*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), mother_for_target_log_, false, id+4);
            new G4PVPlacement(NULL, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-14.5*mm+10*mm + pipe_length_v1/2), target_pipe_v_log_, Form("target_pipe_v1%d", 2*icell), mother_for_target_log_, false, id+4);

            r = base_twoOutterR + pipe_length_h / 2;

            new G4PVPlacement(Rot1A, G4ThreeVector(r*sin(angle1), r*cos(angle1), dz+40*mm+14.5*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell), mother_for_target_log_, false, id+5);
            new G4PVPlacement(Rot2A, G4ThreeVector(-r*sin(angle2), r*cos(angle2), dz-40*mm-14.5*mm), target_pipe_h_log_, Form("target_pipe_h%d", 2*icell+1), mother_for_target_log_, false, id+5);
            /*
            // place probes
            
            r = 61*mm/2;
            new G4PVPlacement(Rot1B, G4ThreeVector(r*sin(angle1), -r*cos(angle1), dz+40*mm+11*mm), target_probe_t_log_, Form("target_probe_t%d", 2*icell), mother_for_target_log_, false, id+6);
            new G4PVPlacement(Rot2B, G4ThreeVector(-r*sin(angle2), -r*cos(angle2), dz-40*mm-11*mm), target_probe_p_log_, Form("target_probe_p%d", 2*icell), mother_for_target_log_, false, id+6);
            */
        };

        //place insulation
        G4RotationMatrix* RotAl = new G4RotationMatrix;
        RotAl->rotateX(-90*deg);
        RotAl->rotateY(-90*deg);
        RotAl->rotateZ(0*deg);

        new G4PVPlacement(RotAl, G4ThreeVector(0,0, target_offset), target_si_log_, "Trapez", mother_for_target_log_, false, 26);

        
        //
        G4Colour target_color;
        G4Colour::GetColour("target", target_color);
        G4Colour target_si_color;
        G4Colour::GetColour("target_si", target_si_color);
        G4Colour target_base_color;
        G4Colour::GetColour("target_base", target_base_color);
        G4Colour target_film_color;
        G4Colour::GetColour("target_film", target_film_color);
        G4Colour target_pipe_color;
        G4Colour::GetColour("target_pipe", target_pipe_color);
        target_log_->SetVisAttributes(new G4VisAttributes(target_color));
        target_base_top_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_base_bottom_log_->SetVisAttributes(new G4VisAttributes(target_base_color));
        target_cell_film_log_->SetVisAttributes(new G4VisAttributes(target_film_color));
        target_nipple_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_pipe_h_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_pipe_v_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        target_si_log_->SetVisAttributes(new G4VisAttributes(target_si_color));
        //target_probe_p_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        //target_probe_t_log_->SetVisAttributes(new G4VisAttributes(target_pipe_color));
        
        target_log_->SetUserLimits(StepLimit_);

    }

}

void g4PSITarget::Info() {
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    
    InfoTitle("Target and Scattering Chamber");
    if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type0])
        InfoParString("Target type", "kTarget_Type0");
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type1])
        InfoParString("Target type", "kTarget_Type1");
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type2])
        InfoParString("Target type", "kTarget_Type2");
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru])
        InfoParString("Target type", "kTarget_TypeJeru");
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM])
            InfoParString("Target type", "kTarget_TypeUM");
    else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUM2])
            InfoParString("Target type", "kTarget_TypeUM2");
    else
        InfoParString("Target type", "unknown");
    
    //    InfoPar2Double("Target cross section (rx,ry)", target_dx_/cm, target_dy_/cm, "cm");
    InfoParDouble("Target/chamber entrance", GetUpstreamPos()/cm, "cm");
}


void g4PSITarget::SetSD(G4SDManager *SDman) {
    
    G4String TargetSDname = "g4PSI/" + label_;
    G4VSensitiveDetector* TargetSD = SDman->FindSensitiveDetector(TargetSDname);
    if (TargetSD == NULL) {
        TargetSD = SD_target_ = new g4PSITargetSD(TargetSDname, label_ + "_Collection", min_p_, min_theta_);
        SDman->AddNewDetector( TargetSD );
    };

    /* G4String SDname = "g4PSI/Scint" + label_;
    G4VSensitiveDetector* ScintSD_ = SDman->FindSensitiveDetector(SDname);
    if (ScintSD_ == NULL) {
        ScintSD_ = new g4PSIScintillatorSD(SDname, 1, label_ + "_Collection" );
        SDman->AddNewDetector( ScintSD_ );
    };

    if (scintillator_bar_A_log_) scintillator_bar_A_log_->SetSensitiveDetector( ScintSD_ );
    if (scintillator_bar_B_log_) scintillator_bar_B_log_->SetSensitiveDetector( ScintSD_ );
    */

    // attach Target sensitive detector to all relevant logical volumes
    
    if (target_log_) target_log_->SetSensitiveDetector( TargetSD );
    if (target_cell_film_log_) target_cell_film_log_->SetSensitiveDetector( TargetSD );
    if (target_si_log_) target_si_log_->SetSensitiveDetector( TargetSD );
    if (target_base_top_log_) target_base_top_log_->SetSensitiveDetector( TargetSD );
    if (target_base_bottom_log_) target_base_bottom_log_->SetSensitiveDetector( TargetSD );
    if (target_nipple_log_) target_nipple_log_->SetSensitiveDetector( TargetSD );
    if (target_pipe_v_log_) target_pipe_v_log_->SetSensitiveDetector( TargetSD );
    if (target_pipe_h_log_) target_pipe_h_log_->SetSensitiveDetector( TargetSD );
    if (target_probe_t_log_) target_probe_t_log_->SetSensitiveDetector( TargetSD );
    if (target_probe_p_log_) target_probe_p_log_->SetSensitiveDetector( TargetSD );
    if (rightframe_log_) rightframe_log_->SetSensitiveDetector( TargetSD );
    if (leftframe_log_) leftframe_log_->SetSensitiveDetector( TargetSD );
    if (Tpipe_log_) Tpipe_log_->SetSensitiveDetector ( TargetSD );
    if (Tconnect_log_) Tconnect_log_->SetSensitiveDetector( TargetSD );
    if (Bpipe_log_) Bpipe_log_->SetSensitiveDetector( TargetSD );
    if (Bconnect_log_) Bconnect_log_->SetSensitiveDetector( TargetSD );


    if (chamber_full_log_) chamber_full_log_->SetSensitiveDetector( TargetSD );
    //    if (chamber_full_vac_log_) chamber_full_vac_log_->SetSensitiveDetector( TargetSD );
    if (chamber_exit_window_log_) chamber_exit_window_log_->SetSensitiveDetector( TargetSD);
    if (chamber_entrance_window_log_) chamber_entrance_window_log_->SetSensitiveDetector( TargetSD );
    if (window_post_log_) window_post_log_->SetSensitiveDetector( TargetSD );
    if (chamber_side_exit_window_log_) chamber_side_exit_window_log_->SetSensitiveDetector( TargetSD );
    if (chamber_back_exit_window_log_) chamber_back_exit_window_log_->SetSensitiveDetector( TargetSD );

}


void g4PSITarget::Write() {
    TVectorD v1(6);
    v1[0] = target_dz_;
    v1[1] = target_dx_;
    v1[2] = target_dy_;
    v1[3] = chamber_outer_radius_;
    v1[4] = chamber_inner_radius_;
    v1[5] = entrance_window_zpos_;
    v1.Write(label_);
}


void g4PSITarget::InitTree(TTree *T) {
    if (SD_target_) SD_target_->InitTree(T);
}


void g4PSITarget::DeleteEventData() {
    if (SD_target_) SD_target_->DeleteEventData();
}


void g4PSITarget::MakeWindow(G4String name, G4VSolid *&window, G4double dx, G4double dy, G4double t, G4double r) {
    
    window = NULL;
    if (t > 0) {
        G4Box *w1 = new G4Box(name + "_w1", dx/2, dy/2 - r, t/2);
        G4Box *w2 = new G4Box(name + "_w2", dx/2 - r, dy/2, t/2);
        window = new G4UnionSolid(name + "_c0", w1, w2);
        
        if (r > 0) {

            
            //  without the increase of r by a small amount (dr)
            //  we get the error message below: (why?)
            //
            //  ERROR: G4VSceneHandler::RequestPrimitives
            //  Polyhedron not available for exit_window_c4.
            //  This means it cannot be visualized on most systems.
            //  Contact the Visualization Coordinator.
            
            G4double tx = dx/2 - r;
            G4double ty = dy/2 - r;
            
            G4Tubs *t1 = new G4Tubs(name + "_tub_c1", 0, r, t/2, 0, pi/2);
            window = new G4UnionSolid(name + "_uni_c1", window, t1, NULL, G4ThreeVector(tx, ty, 0.));
            G4Tubs *t2 = new G4Tubs(name + "_tub_c2", 0, r, t/2, pi/2, pi/2);
            window = new G4UnionSolid(name + "_uni_c2", window, t2, NULL, G4ThreeVector(-tx, ty, 0.));
            G4Tubs *t3 = new G4Tubs(name + "_tub_c3", 0, r, t/2, pi, pi/2);
            window = new G4UnionSolid(name + "_uni_c3", window, t3, NULL, G4ThreeVector(-tx, -ty, 0.));
            G4Tubs *t4 = new G4Tubs(name + "_tub_c4", 0, r, t/2, 3*pi/2, pi/2);
            window = new G4UnionSolid(name + "_uni_c4", window, t4, NULL, G4ThreeVector(tx, -ty, 0.));
            
//            G4ThreeVector trans1(dx/2 - r, dy/2 - r, 0);
//            G4ThreeVector trans2(-(dx/2 - r), dy/2 - r, 0);
//            G4ThreeVector trans3(-(dx/2 - r), -(dy/2 - r), 0);
//            G4ThreeVector trans4((dx/2 - r), -(dy/2 - r), 0);
//            window = new G4UnionSolid(name + "_c1", window, t1, NULL, trans1);
//            window = new G4UnionSolid(name + "_c2", window, t1, NULL, trans2);
//            window = new G4UnionSolid(name + "_c3", window, t1, NULL, trans3);
//            window = new G4UnionSolid(name + "_c4", window, t1, NULL, trans4);
            
        } else {
            window = w1;
        }
    }
}


void g4PSITarget::MakeWindow(G4String name, G4LogicalVolume* &window_log, G4Material* mat, G4double dx, G4double dy, G4double t, G4double r) {
    
    G4VSolid *window = NULL;
    MakeWindow(name, window, dx, dy, t, r);
    window_log = new G4LogicalVolume(window, mat, name + "_log");
}
