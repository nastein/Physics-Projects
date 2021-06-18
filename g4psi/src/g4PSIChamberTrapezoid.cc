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

#include "g4PSIDetectorParts.hh"
#include "g4PSIChamberTrapezoid.hh"
#include "TVectorD.h"


g4PSIChamberTrapezoid::g4PSIChamberTrapezoid(G4String label) : g4PSIChamberBase(label)
{
    chamber_sd_ = NULL;
    
    chamber_full_log_ = NULL;
    chamber_full_vac_log_ = NULL;
    chamber_conflat_log_ = NULL;
    chamber_entrance_window_log_ = NULL;
    chamber_exit_window_kapton_log_ = NULL;
    chamber_side_exit_window_mylar_log_ = NULL;
    chamber_side_exit_window_kevlar_log_ = NULL;
    chamber_side_exit_window_poly_log_ = NULL;
    chamber_veto_log_ = NULL;
    chamber_side_veto_log_ = NULL;
    chamber_downstream_window_reinforcement_log_ = NULL;
    chamber_downstream_window_frame_log_ = NULL;
    chamber_downstream_window_vacuum_log_ = NULL;
    chamber_window_cover_log_ = NULL;
    chamber_left_side_exit_frame_log_ = NULL;
    chamber_right_side_exit_frame_log_ = NULL;
    chamber_side_exit_frame_vacuum_log_ = NULL;
}


void g4PSIChamberTrapezoid::Info() {
    InfoTitle("Chamber (trapezoid)");
}


void g4PSIChamberTrapezoid::Placement() {
    
    using namespace CLHEP;

    g4PSIDetectorParts* detector_parts = g4PSIDetectorParts::getInstance();
    
    //  Materials for chamber and vacuum inside chamber
    G4NistManager* man = G4NistManager::Instance();
    
    G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");
    G4Material* Kapton = man->FindOrBuildMaterial("G4_KAPTON");
    G4Material* SiO2   = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* Polycarb = man->FindOrBuildMaterial("G4_POLYCARBONATE");
    G4Material* Polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
    G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");
    G4Element* C  = man->FindOrBuildElement("C");
    G4Element* Si = man->FindOrBuildElement("Si");
    G4Element* Cr = man->FindOrBuildElement("Cr");
    G4Element* Mn = man->FindOrBuildElement("Mn");
    G4Element* Fe = man->FindOrBuildElement("Fe");
    G4Element* Ni = man->FindOrBuildElement("Ni");
    G4Element* H  = man ->FindOrBuildElement("H");
    G4Element* O = man->FindOrBuildElement("O");
    G4Element* N = man->FindOrBuildElement("N");
    
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

    // Kevlar
    // (-NH-C6H4-NH-CO-C6H4-CO-)*n
    G4String name = "Kevlar";
    double densityKev = 1.44 *g/cm3;
    int ncomponentsKev = 4;
    G4Material* Kevlar = new G4Material(name, densityKev, ncomponentsKev);
    Kevlar -> AddElement(H, natoms=10 );
    Kevlar -> AddElement(C, natoms=14);
    Kevlar -> AddElement(O, natoms= 2);
    Kevlar -> AddElement(N, natoms= 2);
    
    /// ------------- Stainless Steel Chamber -------------
    
    const G4double inches = 2.54 * cm;
    
    G4double theta = 65 * deg;
    G4double thetaR = 90 * deg - theta;
    
    G4double sc_wall_thickness = .375 * inches;
    
    G4double sc_half_height = 26 * inches / 2.0;
    G4double sc_half_y = 14.254 * inches / 2.0 * cos(thetaR) + sc_wall_thickness;
    G4double sc_half_x_small = 4.094 * inches / 2.0;
    G4double sc_half_x_large = sc_half_x_small + 2 * sc_half_y / tan(theta);
    
    G4Trap *chamber_ss_trap = new G4Trap("chamber_ss_trap", sc_half_height, 0*deg, 0*deg, sc_half_y, sc_half_x_small, sc_half_x_large, 0*deg,sc_half_y, sc_half_x_small, sc_half_x_large, 0*deg);
    
    G4double input_stub_width = 25.4 * mm;
    G4double input_stub_outer_radius = 50.8 * mm;
    G4double input_stub_inner_radius =  48.692 * mm;
    G4double stub_conflat_overlap = 8.712 * mm;
    G4double conflat_one_width = 19.812 * mm;
    G4double conflat_two_width = 18.796 * mm;
    G4double conflat_one_outer_radius = 75.819 * mm;
    G4double conflat_one_inner_radius = 75.0 * mm / 2.0;
    G4double conflat_two_outer_radius = 75.819 * mm;
    G4double conflat_two_inner_radius = input_stub_inner_radius;
    G4double entrance_window_thickness = 50 * um;
    G4double entrance_window_sc_position = 16 * inches;

    G4double downstream_exit_window_w = 3.071 * inches;
    G4double downstream_exit_window_h = 14 * inches;
    G4double downstream_exit_window_r = 1.181 * inches;
    G4double downstream_exit_window_thickness = 130*um;

    G4double downstream_frame_r = 1.685 * inches;
    G4double downstream_frame_w = 103.985 * mm;
    G4double downstream_frame_h = 381.2 * mm;
    G4double downstream_frame_thickness = .188 * inches;
    G4double downstream_block_w = 43.762 * mm;
    G4double downstream_block_h = 12.802 * mm;
    G4double downstream_block_thickness = 11.708 * mm;

    
    G4double side_exit_window_w = 11.867 * inches;
    G4double side_exit_window_h = downstream_exit_window_h;
    G4double side_exit_window_r = 1.969 * inches;
    G4double side_exit_window_thickness = 120*um;
    G4double mylar_window_thickness = 125*um;
    G4double kevlar_window_thickness = 175*um;
    G4double poly_window_thickness = .25 * mm;

    G4double side_window_frame_bar_thickness = .188 * inches;
    G4double side_window_frame_bar_width = .505 * inches;
    G4double side_window_frame_w = side_exit_window_w + 2*side_window_frame_bar_width;
    G4double side_window_frame_h = side_exit_window_h + 2*side_window_frame_bar_width;
    G4double side_window_frame_r = 2.473 * inches;

    G4double frame_w = sc_half_x_small - downstream_exit_window_w/2;
    
    G4double chamber_extension_r = 10.0 *cm;
    G4double chamber_extension_h = 5*cm;
    
    // ---------------------------------------------------------------------------------------------------

    G4Tubs *entrance_conflat_one = new G4Tubs("entrance_conflat_one", conflat_one_inner_radius, conflat_one_outer_radius, conflat_one_width / 2.0, 0, 2 * pi);
    G4Tubs *entrance_conflat_two = new G4Tubs("entrance_conflat_two", conflat_two_inner_radius, conflat_two_outer_radius, conflat_two_width / 2.0, 0, 2 * pi);
    G4VSolid* entrance_conflat = new G4UnionSolid("entrance_conflat", entrance_conflat_one, entrance_conflat_two, NULL, G4ThreeVector(0,0, + (conflat_one_width + conflat_two_width) / 2.0));

    G4Tubs *conflat_one_vacuum = new G4Tubs("conflat_one_vacuum", 0, conflat_one_inner_radius, conflat_one_width / 2.0, 0, 2* pi);
    G4Tubs *conflat_two_vacuum = new G4Tubs("conflat_two_vacuum", 0, conflat_two_inner_radius, (conflat_two_width - stub_conflat_overlap) / 2.0, 0, 2* pi);
    G4VSolid* conflat_vacuum = new G4UnionSolid("conflat_vacuum", conflat_one_vacuum, conflat_two_vacuum, NULL, G4ThreeVector(0,0, + (conflat_one_width + conflat_two_width - stub_conflat_overlap) / 2.0));

    G4Tubs *entrance_ss_tube = new G4Tubs("entrance_ss_tube", 0, input_stub_outer_radius, input_stub_width / 2.0, 0, 2 * pi);
    G4RotationMatrix *rot_entrance_window = new G4RotationMatrix;
    rot_entrance_window->rotateX(90 * deg);
    G4VSolid *full_chamber = new G4UnionSolid("chamber_ss", chamber_ss_trap, entrance_ss_tube, rot_entrance_window, G4ThreeVector(0, sc_half_y + input_stub_width / 2, entrance_window_sc_position - sc_half_height));

    G4double sc_z_shift = (8.386 * inches + sc_wall_thickness) / 2.0 - sc_half_y;
    
    G4Tubs *top_cyl = new G4Tubs("top_ss_tube", 0, chamber_extension_r + sc_wall_thickness, chamber_extension_h, 0, 2*pi);
    
    full_chamber = new G4UnionSolid("full_chamber_vac + entrance window", full_chamber, top_cyl, NULL, G4ThreeVector(0, -sc_z_shift, sc_half_height + chamber_extension_h/2.0));
    
    chamber_full_log_ = new G4LogicalVolume(full_chamber, StainlessSteel, "full_chamber.log");


    // ------------- Downstream exit window frame --------
    G4VSolid* chamber_downstream_window_frame = NULL;
    G4VSolid* chamber_downstream_window_vacuum = NULL;

    MakeWindow("downstream_exit_window_frame", chamber_downstream_window_frame, downstream_frame_w, downstream_frame_h, downstream_frame_thickness, downstream_frame_r);
    MakeWindow("downstream_exit_window_vacuum", chamber_downstream_window_vacuum, downstream_exit_window_w, downstream_exit_window_h, downstream_frame_thickness, downstream_exit_window_r );

    chamber_downstream_window_frame_log_ = new G4LogicalVolume(chamber_downstream_window_frame, StainlessSteel, "chamber_downstream_window_frame_log");
    chamber_downstream_window_vacuum_log_ = new G4LogicalVolume(chamber_downstream_window_vacuum, Vacuum, "chamber_downstream_window_vacuum_log");   

    G4Box* chamber_downstream_block = new G4Box("downstream_block", downstream_block_h/2, downstream_block_w/2, downstream_block_thickness/2);
    G4Trap* chamber_downstream_trap = new G4Trap("downstream_trap", downstream_block_h, 42.8 * mm, downstream_block_thickness, 3.827 * mm);

    G4double middle_trap_side = (12.689*(11.708 - 3.827) + 3.827*42.8)/42.8;
    G4Trap* chamber_downstream_long_trap = new G4Trap("downstream_long_trap", 381.203*mm, 12.689*mm, middle_trap_side * mm, 3.827 *mm);

    G4RotationMatrix *traprot = new G4RotationMatrix;
    traprot->rotateY(-90*deg);
    G4RotationMatrix *traprotreverse = new G4RotationMatrix;
    traprotreverse->rotateY(-90*deg);
    traprotreverse->rotateX(180*deg);

    G4double block_shift_y = .52 * mm;
    G4double block_shift_z = 1.98 * mm;

    G4VSolid* reinforced_frame_top = new G4UnionSolid("reinforced_downstream_frame_top", chamber_downstream_block, chamber_downstream_trap,traprot, G4ThreeVector(0, +downstream_block_w - block_shift_y, + block_shift_z));
    reinforced_frame_top = new G4UnionSolid("reinforced_downstream_frame_top", reinforced_frame_top, chamber_downstream_trap,traprotreverse, G4ThreeVector(0, -downstream_block_w + block_shift_y, + block_shift_z));

    G4VSolid* reinforced_frame = new G4UnionSolid("reinforced_downstream_frame", reinforced_frame_top, chamber_downstream_long_trap, traprot, G4ThreeVector(downstream_frame_h/2 - downstream_block_h/2, +129.362*mm/2 - 12.689*mm/2, + 3.827*mm - .5*mm));

    reinforced_frame = new G4UnionSolid("reinforced_downstream_frame", reinforced_frame, chamber_downstream_long_trap, traprotreverse, G4ThreeVector(downstream_frame_h/2 - downstream_block_h/2, -129.362*mm/2 + 12.689*mm/2, +3.827*mm - .5 *mm));

    reinforced_frame = new G4UnionSolid("reinforced_downstream_frame", reinforced_frame_top, reinforced_frame, NULL, G4ThreeVector(-downstream_frame_h,0,0));

    chamber_downstream_window_reinforcement_log_ = new G4LogicalVolume(reinforced_frame, StainlessSteel, "reinforced_frame_log");

   
    // ------------ Window covers ------------

    G4double window_cover_width = 4*cm;
    G4double window_cover_t = 3*mm;
        
    G4Box *window_cover = new G4Box("window_cover", downstream_exit_window_h/2, window_cover_width/2, window_cover_t/2);
    chamber_window_cover_log_ = new G4LogicalVolume(window_cover, Polycarb, "window_cover.log");
    
    // ------------- Vacuum --------------
    
    G4double dx_s = sc_wall_thickness * (-1. * tan(theta/2));
    G4double dx_l = sc_wall_thickness * (-1. / tan(theta/2));
    
    G4Trap *chamber_vac_trap = new G4Trap("chamber_vac_trap", sc_half_height - sc_wall_thickness, 0*deg, 0*deg, sc_half_y - sc_wall_thickness, sc_half_x_small + dx_s, sc_half_x_large + dx_l, 0*deg, sc_half_y - sc_wall_thickness, sc_half_x_small + dx_s, sc_half_x_large + dx_l, 0*deg);
    G4Tubs *entrance_vac_tube = new G4Tubs("entrance_vac_tube", 0, input_stub_inner_radius, (input_stub_width + sc_wall_thickness) / 2.0, 0, 2 * pi);
    G4Tubs *top_vac_cyl = new G4Tubs("top_vac_tube", 0, chamber_extension_r, chamber_extension_h, 0, 2*pi);
    
    G4VSolid *full_chamber_vac = new G4UnionSolid("chamber_vac", chamber_vac_trap, entrance_vac_tube, rot_entrance_window, G4ThreeVector(0, sc_half_y + (input_stub_width - sc_wall_thickness) / 2, entrance_window_sc_position - sc_half_height));
    full_chamber_vac = new G4UnionSolid("full_chamber_vac", full_chamber_vac, top_vac_cyl, NULL, G4ThreeVector(0, -sc_z_shift, sc_half_height - sc_wall_thickness + chamber_extension_h/2.0));
    
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
    
    // ------------- Window foils --------------
    
    G4Tubs *entrance_window = new G4Tubs("entrance_window", 0, conflat_one_inner_radius, entrance_window_thickness/2., 0, 2 * pi);
    chamber_entrance_window_log_ = new G4LogicalVolume(entrance_window, Kapton, "entrance_window.log");
    chamber_conflat_log_ = new G4LogicalVolume(entrance_conflat, StainlessSteel, "chamber_conflat.log");
    chamber_conflat_vacuum_log_ = new G4LogicalVolume(conflat_vacuum, Vacuum, "conflat_vacuum.log");
    
    chamber_side_exit_window_mylar_log_ = NULL;
    chamber_side_exit_window_kevlar_log_ = NULL;
    chamber_side_exit_window_poly_log_ = NULL;
    chamber_exit_window_kapton_log_ = NULL;
    MakeWindow("side_exit_window_mylar", chamber_side_exit_window_mylar_log_, Mylar, side_exit_window_w, side_exit_window_h, mylar_window_thickness, side_exit_window_r);
    MakeWindow("side_exit_window_kevlar", chamber_side_exit_window_kevlar_log_, Kevlar, side_exit_window_w, side_exit_window_h, kevlar_window_thickness, side_exit_window_r);
    //MakeWindow("side_exit_window_poly", chamber_side_exit_window_poly_log_, Polyethylene, side_exit_window_w, side_exit_window_h, poly_window_thickness, side_exit_window_r);
    MakeWindow("downstream_exit_kapton", chamber_exit_window_kapton_log_, Kapton, downstream_exit_window_w, downstream_exit_window_h, downstream_exit_window_thickness, downstream_exit_window_r);
    
    // ------------ side exit window frames --------------
    G4VSolid *side_window_frame_ = NULL;
    G4VSolid *side_window_frame_vacuum_ = NULL;

    MakeWindow("side_frame", side_window_frame_, side_window_frame_w, side_window_frame_h, side_window_frame_bar_thickness, side_window_frame_r);
    MakeWindow("side_frame_minus", side_window_frame_vacuum_, side_exit_window_w, side_exit_window_h, side_window_frame_bar_thickness, side_exit_window_r);  
    chamber_left_side_exit_frame_log_ = new G4LogicalVolume(side_window_frame_, StainlessSteel, "left_side_exit_frame_log");
    chamber_right_side_exit_frame_log_ = new G4LogicalVolume(side_window_frame_, StainlessSteel, "right_side_exit_frame_log");
    chamber_side_exit_frame_vacuum_log_ = new G4LogicalVolume(side_window_frame_vacuum_, Vacuum, "side_exit_frame_vacuum_log");

    // ------ veto detectors  
    /*
    G4double sc_width = .504 * inches;
    G4double sc_bar_length = (downstream_exit_window_h/2 + frame_w) * 2;
    G4double sc_thickness = 5.0*mm;

    G4double sc_width_ds_frame = frame_w + .5 * inches - window_cover_t;
    G4double sc_width_us_frame = .5 * inches;
    
    G4Material* sc_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4Box* veto_bar = new G4Box("veto", sc_width/2, sc_thickness/2, sc_bar_length/2);
    chamber_veto_log_ = new G4LogicalVolume(veto_bar, sc_material, "veto_log", 0, 0, 0);

    G4Box* veto_bar_ds_frame = new G4Box("veto_ds_frame", sc_width_ds_frame/2, sc_thickness/2, sc_bar_length/2 );
    G4Box* veto_bar_us_frame = new G4Box("veto_us_frame", sc_width_us_frame/2, sc_thickness/2, sc_bar_length/2 );
    chamber_veto_log_ds_ = new G4LogicalVolume(veto_bar_ds_frame, sc_material, "veto_log_ds", 0, 0, 0);
    chamber_veto_log_us_ = new G4LogicalVolume(veto_bar_us_frame, sc_material, "veto_log_us", 0, 0, 0);
    */
    G4double sc_width = 12.1*mm;
    G4double sc_bar_length = 300*mm;
    //G4double sc_bar_length = side_window_frame_h;
    G4double sc_thickness = 2.5*mm;

    G4double sc_sidebar_width = 0; 

    G4Material* sc_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4Box* veto_bar = new G4Box("veto", sc_width/2, sc_thickness/2, sc_bar_length/2);
    chamber_veto_log_ = new G4LogicalVolume(veto_bar, sc_material, "veto_log", 0, 0, 0);


    // ------ placements ------------------------------------------------------------------------------------------------
    
    G4RotationMatrix *rot_chamber = new G4RotationMatrix;
    rot_chamber->rotateX(90 * deg);
    
    G4double sc_y_shift = -(entrance_window_sc_position - sc_half_height);

    //Comment this out if you want to do a simulation without the chamber
    new G4PVPlacement(rot_chamber, G4ThreeVector(0, sc_y_shift, -sc_z_shift), chamber_full_log_, "chamber_full", mother_volume_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0,0), chamber_full_vac_log_, "chamber_full_vac", chamber_full_log_, false, 1);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, -sc_half_y - input_stub_width - sc_z_shift - (conflat_one_width + conflat_two_width)/2), chamber_conflat_log_, "chamber_conflat", mother_volume_, false, 2);
    //-----------------------------------------------
    //Uncomment this if you want to do a simulation without the chamber
    //new G4PVPlacement(rot_chamber, G4ThreeVector(0,sc_y_shift, -sc_z_shift), chamber_full_vac_log_, "chamber_full_vac", mother_volume_, false, 0);
    //-----------------------------------------------
    new G4PVPlacement(NULL, G4ThreeVector(0,0, -sc_half_y - input_stub_width - sc_z_shift - (conflat_one_width + conflat_two_width)/2), chamber_conflat_vacuum_log_, "chamber_conflat_vacuum", mother_volume_, false, 57);
    
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, -(conflat_one_width)/2 + entrance_window_thickness/2), chamber_entrance_window_log_, "entrance_window", chamber_conflat_vacuum_log_, false, 3);
    
    new G4PVPlacement(rot_entrance_window, G4ThreeVector(0,-sc_half_y + downstream_exit_window_thickness/2.0, entrance_window_sc_position - sc_half_height), chamber_exit_window_kapton_log_, "downstream_exit_window_kapton", chamber_full_vac_log_, false, 4);
  
    //Place side exit windows

    //G4double wshift_mylar = (sc_wall_thickness - side_exit_window_thickness) / 2.0;
    G4double wshift_mylar = (sc_wall_thickness - mylar_window_thickness) / 2.0;
    G4double wshift_kevlar = (sc_wall_thickness - mylar_window_thickness - kevlar_window_thickness) / 2.0;
    G4double wshift_poly = (sc_wall_thickness - mylar_window_thickness + poly_window_thickness) / 2.0;

    new G4PVPlacement(rot_side_window_1, G4ThreeVector(xwin + wshift_mylar * cos(thetaR), ywin - wshift_mylar * sin(thetaR), zwin), chamber_side_exit_window_mylar_log_, "side_exit_window_mylar1", chamber_full_vac_log_, false, 5);
    new G4PVPlacement(rot_side_window_1, G4ThreeVector(xwin + wshift_kevlar * cos(thetaR), ywin - wshift_kevlar * sin(thetaR), zwin), chamber_side_exit_window_kevlar_log_, "side_exit_window_kevlar1", chamber_full_vac_log_, false, 6);
    //new G4PVPlacement(rot_side_window_1, G4ThreeVector(xwin + wshift_poly * cos(thetaR), ywin - wshift_poly * sin(thetaR), zwin), chamber_side_exit_window_poly_log_, "side_exit_window_poly1", chamber_full_vac_log_, false, 57);
    new G4PVPlacement(rot_side_window_2, G4ThreeVector(-xwin - wshift_mylar * cos(thetaR), ywin - wshift_mylar * sin(thetaR), zwin), chamber_side_exit_window_mylar_log_, "side_exit_window_mylar2", chamber_full_vac_log_, false, 7);
    new G4PVPlacement(rot_side_window_2, G4ThreeVector(-xwin - wshift_kevlar * cos(thetaR), ywin - wshift_kevlar * sin(thetaR), zwin), chamber_side_exit_window_kevlar_log_, "side_exit_window_kevlar2", chamber_full_vac_log_, false, 8);
    //new G4PVPlacement(rot_side_window_2, G4ThreeVector(-xwin - wshift_poly * cos(thetaR), ywin - wshift_poly * sin(thetaR), zwin), chamber_side_exit_window_poly_log_, "side_exit_window_poly2", chamber_full_vac_log_, false, 58);

    //Place side exit window frames
    G4double wshift_window_frame = (sc_wall_thickness + side_window_frame_bar_thickness)/2;
    
    //Comment this out if you want to do a simulation without the chamber
    new G4PVPlacement(rot_side_window_1 , G4ThreeVector(xwin + wshift_window_frame * cos(thetaR), ywin - wshift_window_frame * sin(thetaR), zwin), chamber_left_side_exit_frame_log_, "side_exit_window_frame1", chamber_full_vac_log_, false, 50);
    new G4PVPlacement(NULL , G4ThreeVector(0,0,0), chamber_side_exit_frame_vacuum_log_, "side_exit_window_frame_vacuum1", chamber_left_side_exit_frame_log_, false, 51);
    new G4PVPlacement(rot_side_window_2 , G4ThreeVector(-xwin - wshift_window_frame * cos(thetaR), ywin - wshift_window_frame * sin(thetaR), zwin), chamber_right_side_exit_frame_log_, "side_exit_window_frame2", chamber_full_vac_log_, false, 52);
    new G4PVPlacement(NULL , G4ThreeVector(0,0,0), chamber_side_exit_frame_vacuum_log_, "side_exit_window_frame_vacuum2", chamber_right_side_exit_frame_log_, false, 53);
    
    
    //Place veto scintillators on the sides of the chamber
    G4RotationMatrix *rot_veto_a = new G4RotationMatrix;
    rot_veto_a->rotateX(90 * deg);
    rot_veto_a->rotateZ(thetaR + 90 * deg);
    G4RotationMatrix *rot_veto_b = new G4RotationMatrix;
    rot_veto_b->rotateX(90 * deg);
    rot_veto_b->rotateZ(-thetaR + 90 * deg);
    G4double delta = 0.504 * inches;
    G4double x_ = +sc_half_x_small - cos(theta) * (sc_width / 2 - delta) + sin(theta) * (sc_thickness / 2);
    G4double y_ = -sc_y_shift - entrance_window_sc_position + sc_half_height;
    G4double z_ = -sc_z_shift + sc_half_y + sin(theta) * (sc_width / 2 - delta) + cos(theta) * (sc_thickness / 2);
    
    new G4PVPlacement(rot_veto_a, G4ThreeVector(x_ + side_window_frame_bar_thickness*sin(theta), y_, z_ + side_window_frame_bar_thickness*cos(theta)), chamber_veto_log_, "veto 0", mother_volume_, false, 0);
    new G4PVPlacement(rot_veto_b, G4ThreeVector(-x_ - side_window_frame_bar_thickness*sin(theta), y_, z_ + side_window_frame_bar_thickness*cos(theta)), chamber_veto_log_, "veto 1", mother_volume_, false, 1);

    //Veto scintillators on the downstream exit window frame
    /*
    G4RotationMatrix *rot_veto1 = new G4RotationMatrix;
    rot_veto1->rotateX(90*deg);

    G4double veto_angle1 = 15.039 * deg;
    G4RotationMatrix *rot_veto2 = new G4RotationMatrix;
    rot_veto2->rotateX(90*deg);
    rot_veto2->rotateZ(-veto_angle1);

    G4double veto_angle2 = 26.556 * deg;
    G4RotationMatrix *rot_veto3 = new G4RotationMatrix;
    rot_veto3->rotateX(90*deg);
    rot_veto3->rotateZ(-veto_angle1 - veto_angle2);

    G4RotationMatrix *rot_veto4 = new G4RotationMatrix;
    rot_veto4->rotateX(90*deg);
    rot_veto4->rotateZ(veto_angle1);

    G4RotationMatrix* rot_veto5 = new G4RotationMatrix;
    rot_veto5->rotateX(90*deg);
    rot_veto5->rotateZ(veto_angle1 + veto_angle2);


    G4double x = downstream_exit_window_w/2 + sc_width/2;
    G4double y = 0;
    G4double z = -sc_z_shift + sc_half_y + 14.985*mm + sc_thickness/2;

    //The shifts in the positions of the veto scintillators come from Tiko's autocad design, they are the distance between the centers of consecutive scintillator bars
    new G4PVPlacement(rot_veto2, G4ThreeVector(x + 12.282*mm,y,z -1.62123*mm), chamber_veto_log_, "veto 1", mother_volume_, false, 1);
    new G4PVPlacement(rot_veto2, G4ThreeVector(x + 12.282*mm + 11.68555*mm,y,z - 1.62123*mm - 3.13971*mm), chamber_veto_log_, "veto 2", mother_volume_, false, 2);
    new G4PVPlacement(rot_veto3, G4ThreeVector(x + 12.282*mm + 11.68555*mm + 10.97384 * mm,y,z -1.62123*mm - 3.13971*mm - 5.91312*mm), chamber_veto_log_, "veto 3", mother_volume_, false, 3);

    //new G4PVPlacement(rot_veto1, G4ThreeVector(-x,y,z), chamber_veto_log_, "veto 5", mother_volume_, false, 7);
    new G4PVPlacement(rot_veto4, G4ThreeVector(-x - 12.282*mm,y,z -1.62123*mm), chamber_veto_log_, "veto 5", mother_volume_, false, 5);
    new G4PVPlacement(rot_veto4, G4ThreeVector(-x - 12.282*mm - 11.68555*mm,y,z - 1.62123*mm - 3.13971*mm), chamber_veto_log_, "veto 6", mother_volume_, false, 6);
    new G4PVPlacement(rot_veto5, G4ThreeVector(-x - 12.282*mm - 11.68555*mm - 10.97384 * mm,y,z - 1.62123*mm - 3.13971*mm - 5.91312*mm), chamber_veto_log_, "veto 7", mother_volume_, false, 7);
    */

    //Comment this out if you want to do a simulation without the chamber
    //Place downstream exit window frame and reinforcement
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, - sc_z_shift + sc_half_y + (downstream_frame_thickness)/2), chamber_downstream_window_frame_log_, "downstream_window_frame", mother_volume_, false, 54);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), chamber_downstream_window_vacuum_log_, "downstream_window_vacuum", chamber_downstream_window_frame_log_, false, 55);

    G4RotationMatrix *frame_rot = new G4RotationMatrix;
    frame_rot->rotateZ(-90*deg);
    frame_rot->rotateX(180*deg);

    //Comment this out if you want to do a simulation without the chamber
    new G4PVPlacement(frame_rot, G4ThreeVector(0, downstream_frame_h/2, - sc_z_shift + sc_half_y + downstream_frame_thickness + downstream_block_thickness/2), chamber_downstream_window_reinforcement_log_, "reinforced frame", mother_volume_, false, 56);

    // ------ vis
    
    G4Colour colour_chamber_window;
    G4Colour colour_chamber;
    G4Colour colour_chamber_vacuum;
    G4Colour colour_veto;
    G4Colour colour_window_cover;
    G4Colour::GetColour("scattering_chamber_window", colour_chamber_window);
    G4Colour::GetColour("scattering_chamber", colour_chamber);
    G4Colour::GetColour("scattering_chamber_vacuum", colour_chamber_vacuum);
    G4Colour::GetColour("scwall", colour_veto);
    G4Colour::GetColour("gray", colour_window_cover);
    G4VisAttributes* scatWindowVisAtt = new G4VisAttributes(colour_chamber_window);
    G4VisAttributes* scatSteelAtt = new G4VisAttributes(colour_chamber);
    G4VisAttributes* scatVacuumVisAtt = new G4VisAttributes(colour_chamber_vacuum);
    G4VisAttributes* scatVetoVisAtt = new G4VisAttributes(colour_veto);
    //G4VisAttributes* scatWindowFrameAtt = new G4VisAttributes(colour_chamber_vacuum);
    
    scatVacuumVisAtt->SetVisibility(true);
    scatSteelAtt->SetVisibility(true);
    scatWindowVisAtt->SetVisibility(true);
    scatVetoVisAtt->SetVisibility(true);
    //scatWindowFrameAtt->SetVisibility(false);
    
    chamber_full_log_->SetVisAttributes(scatSteelAtt);
    chamber_full_vac_log_->SetVisAttributes(scatVacuumVisAtt);
    chamber_conflat_log_->SetVisAttributes(scatSteelAtt);
    chamber_conflat_vacuum_log_->SetVisAttributes(scatVacuumVisAtt);
    chamber_exit_window_kapton_log_->SetVisAttributes(scatWindowVisAtt);
    chamber_entrance_window_log_->SetVisAttributes(scatWindowVisAtt);
    chamber_side_exit_window_mylar_log_->SetVisAttributes(scatWindowVisAtt);
    chamber_side_exit_window_kevlar_log_->SetVisAttributes(scatWindowVisAtt);
    //chamber_side_exit_window_poly_log_->SetVisAttributes(scatWindowVisAtt);
    chamber_veto_log_->SetVisAttributes(scatVetoVisAtt);
    //chamber_veto_log_us_->SetVisAttributes(scatVetoVisAtt);
    //chamber_veto_log_ds_->SetVisAttributes(scatVetoVisAtt);
    chamber_downstream_window_reinforcement_log_->SetVisAttributes(scatSteelAtt);
    chamber_downstream_window_frame_log_->SetVisAttributes(scatSteelAtt);
    chamber_left_side_exit_frame_log_->SetVisAttributes(scatSteelAtt);
    chamber_right_side_exit_frame_log_->SetVisAttributes(scatSteelAtt);
    //chamber_side_exit_frame_vacuum_log_->SetVisAttributes(scatWindowFrameAtt);

    detector_parts->AttachBiasingOperator(chamber_full_log_);
    detector_parts->AttachBiasingOperator(chamber_full_vac_log_);
    detector_parts->AttachBiasingOperator(chamber_conflat_log_);
    detector_parts->AttachBiasingOperator(chamber_conflat_vacuum_log_);
    detector_parts->AttachBiasingOperator(chamber_entrance_window_log_);
    detector_parts->AttachBiasingOperator(chamber_exit_window_kapton_log_);
    detector_parts->AttachBiasingOperator(chamber_side_exit_window_mylar_log_);
    detector_parts->AttachBiasingOperator(chamber_side_exit_window_kevlar_log_);
    //detector_parts->AttachBiasingOperator(chamber_side_exit_window_poly_log_);
    detector_parts->AttachBiasingOperator(chamber_downstream_window_reinforcement_log_);
    detector_parts->AttachBiasingOperator(chamber_downstream_window_frame_log_);
    detector_parts->AttachBiasingOperator(chamber_left_side_exit_frame_log_);
    detector_parts->AttachBiasingOperator(chamber_right_side_exit_frame_log_);
    detector_parts->AttachBiasingOperator(chamber_side_exit_frame_vacuum_log_);
 
    chamber_full_log_->SetUserLimits(GetChamberStepLimit());
    chamber_full_vac_log_->SetUserLimits(GetChamberStepLimit());
    chamber_conflat_log_->SetUserLimits(GetChamberStepLimit());
    chamber_conflat_vacuum_log_->SetUserLimits(GetChamberStepLimit());
    chamber_veto_log_->SetUserLimits(GetChamberStepLimit());
    chamber_downstream_window_reinforcement_log_->SetUserLimits(GetChamberStepLimit());
    chamber_downstream_window_frame_log_->SetUserLimits(GetChamberStepLimit());
    chamber_left_side_exit_frame_log_->SetUserLimits(GetChamberStepLimit());
    chamber_right_side_exit_frame_log_->SetUserLimits(GetChamberStepLimit());
    chamber_side_exit_frame_vacuum_log_->SetUserLimits(GetChamberStepLimit());

    
    // determine where and how the target is to be placed:
    // --------------------------------------------------
    
    SetMotherLVforTarget(chamber_full_vac_log_);
    SetChamberBeamPosition(G4ThreeVector(0,-sc_z_shift,entrance_window_sc_position - sc_half_height));
    SetChamberBeamRotation(new G4RotationMatrix(G4ThreeVector(1,0,0), G4ThreeVector(0,0,-1), G4ThreeVector(0,1,0)));
}


void g4PSIChamberTrapezoid::SetSD(G4SDManager *SDman) {
    
    G4String ChamberSDname = "g4PSI/" + label_;
    G4VSensitiveDetector* sd = SDman->FindSensitiveDetector(ChamberSDname);
    if (sd == NULL) {
        sd = chamber_sd_ = new g4PSITargetSD(ChamberSDname, label_ + "_Collection", GetChamberSDpMin(), GetChamberSDthetaMin());
        SDman->AddNewDetector( sd );
    };
    
    if (chamber_full_log_) chamber_full_log_->SetSensitiveDetector( sd );
    if (chamber_full_vac_log_) chamber_full_vac_log_->SetSensitiveDetector( sd );
    if (chamber_exit_window_kapton_log_) chamber_exit_window_kapton_log_->SetSensitiveDetector( sd );
    if (chamber_side_exit_window_mylar_log_) chamber_side_exit_window_mylar_log_->SetSensitiveDetector( sd );
    if (chamber_side_exit_window_kevlar_log_) chamber_side_exit_window_kevlar_log_->SetSensitiveDetector( sd );
    //if (chamber_side_exit_window_poly_log_) chamber_side_exit_window_poly_log_->SetSensitiveDetector( sd );
    if (chamber_entrance_window_log_) chamber_entrance_window_log_->SetSensitiveDetector( sd );
    if (chamber_conflat_log_) chamber_conflat_log_->SetSensitiveDetector( sd );
    if (chamber_conflat_vacuum_log_) chamber_conflat_vacuum_log_->SetSensitiveDetector( sd );
    if (chamber_left_side_exit_frame_log_) chamber_left_side_exit_frame_log_->SetSensitiveDetector ( sd );
    if (chamber_right_side_exit_frame_log_) chamber_left_side_exit_frame_log_->SetSensitiveDetector ( sd );
    if (chamber_side_exit_frame_vacuum_log_) chamber_side_exit_frame_vacuum_log_->SetSensitiveDetector ( sd );
    if (chamber_downstream_window_reinforcement_log_) chamber_downstream_window_reinforcement_log_->SetSensitiveDetector( sd );
    if (chamber_downstream_window_frame_log_) chamber_downstream_window_frame_log_->SetSensitiveDetector ( sd );
    if (chamber_downstream_window_vacuum_log_) chamber_downstream_window_vacuum_log_->SetSensitiveDetector ( sd );
    
    if (chamber_veto_log_) {
        G4String ChamberSDnameV = "g4PSI/" + label_ + "_VETO";
        G4VSensitiveDetector* veto_sd = SDman->FindSensitiveDetector(ChamberSDnameV);
        if (veto_sd == NULL) {
            veto_sd = veto_sd_ = new g4PSIScintillatorSD(ChamberSDnameV, 2, label_ + "_VETO_Collection", true );
            SDman->AddNewDetector( veto_sd );
        };
        chamber_veto_log_->SetSensitiveDetector( veto_sd );
        //chamber_veto_log_us_->SetSensitiveDetector( veto_sd );
        //chamber_veto_log_ds_->SetSensitiveDetector( veto_sd );

    } 
    
}   


void g4PSIChamberTrapezoid::Write() {
    TVectorD v1(1);
    v1[0] = 0;
    v1.Write(label_);
}


void g4PSIChamberTrapezoid::InitTree(TTree *T) {
    if (chamber_sd_) chamber_sd_->InitTree(T);
    if (veto_sd_) veto_sd_->InitTree(T);
}


void g4PSIChamberTrapezoid::DeleteEventData() {
    if (chamber_sd_) chamber_sd_->DeleteEventData();
    if (veto_sd_) veto_sd_->DeleteEventData();
}



void g4PSIChamberTrapezoid::MakeWindow(G4String name, G4VSolid *&window, G4double dx, G4double dy, G4double t, G4double r) {
    
    using namespace CLHEP;
    
    window = NULL;
    if (t > 0) {
        G4Box *w1 = new G4Box(name + "_w1", dx/2, dy/2 - r, t/2);
        G4Box *w2 = new G4Box(name + "_w2", dx/2 - r, dy/2, t/2);
        window = new G4UnionSolid(name + "_c0", w1, w2);
        
        if (r > 0) {
            
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
            
        } else {
            window = w1;
        }
    }
}


void g4PSIChamberTrapezoid::MakeWindow(G4String name, G4LogicalVolume* &window_log, G4Material* mat, G4double dx, G4double dy, G4double t, G4double r) {
    
    G4VSolid *window = NULL;
    MakeWindow(name, window, dx, dy, t, r);
    window_log = new G4LogicalVolume(window, mat, name + "_log");
}

