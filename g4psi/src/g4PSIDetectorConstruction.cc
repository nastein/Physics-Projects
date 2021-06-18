#include "g4PSIDetectorConstruction.hh"
#include "g4PSIDetectorParts.hh"
#include "g4PSIDetectorMessenger.hh"
#include "g4PSIDetectorBase.hh"

//#include "g4PSIWC.hh"
#include "g4PSISTT.hh"
#include "g4PSITestPlane.hh"
#include "g4PSITestSphere.hh"
#include "g4PSISCWall.hh"
#include "g4PSISPS.hh"
#include "g4PSIGEM.hh"
#include "g4PSIBeamSC.hh"
#include "g4PSIVeto.hh"
#include "g4PSIBeamMonitor.hh"
#include "g4PSIBeamCherenkov.hh"
#include "g4PSITarget.hh"
#include "g4PSINaI.hh"
#include "g4PSICal.hh"
#include "g4PSIChamberCylinder.hh"
#include "g4PSIChamberTrapezoid.hh"
#include "g4PSITargetCylinder.hh"
#include "g4PSITargetJeru.hh"
#include "g4PSITargetGrid.hh"
#include "g4PSITargetUM.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "globals.hh"

/*
#include "G4GDMLParser.hh"
#include <sys/stat.h>
#include <stdlib.h>
*/

using namespace CLHEP;

g4PSIDetectorConstruction::g4PSIDetectorConstruction()
:  experimentalHall_log_(0), experimentalHall_phys_(0)
{
    // default detector parts to be constructed
    // ----------------------------------------
    
    /// \todo look into
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    // we should not hard-wire which detector components are being included;
    // this should be determined in the *.mac files.
    DetectorParts->AllOff();
    
    //-------------------//
    //      Colors       //
    //-------------------//
    
#ifdef OPAQUE
    const double alpha = 1;
#else
    const double alpha = 0.5;
#endif

    G4Colour::AddToMap("gem_detector", G4Colour(0.23,0.17,0.76));
    G4Colour::AddToMap("gem_detector_frame", G4Colour(166./255,90./255,0./255));
    G4Colour::AddToMap("gem_detector_gas", G4Colour(0.50,0.34,1.00));
    G4Colour::AddToMap("stc_rfshield", G4Colour(0.3,0.3,0.3,0.3));
    G4Colour::AddToMap("stc_testplane", G4Colour(1,0,0));
    G4Colour::AddToMap("stc_frame", G4Colour(0.00,0.29,0.67));
    G4Colour::AddToMap("stc_detector", G4Colour(64/255.,141/255.,173/255.));
    G4Colour::AddToMap("stc_detector_gas", G4Colour(43/255.,199/255.,214/255.));

    G4Colour::AddToMap("scattering_chamber", G4Colour(0.4, 0.7, 0.7, alpha));
    G4Colour::AddToMap("scattering_chamber_window", G4Colour(1.0, 0.0, 0.0, alpha));
    G4Colour::AddToMap("scattering_chamber_airgap", G4Colour(0.4, 0.8, 1.0, alpha));
    G4Colour::AddToMap("target", G4Colour(0.09, 0.52, 0.18, alpha));
    G4Colour::AddToMap("target_base", G4Colour(0.20,0.65,0.70, alpha));
    G4Colour::AddToMap("target_film", G4Colour(0.60,0.78,0.21, alpha));
    G4Colour::AddToMap("target_pipe", G4Colour(0.09,0.40,0.45, alpha));
    G4Colour::AddToMap("target_si", G4Colour(0.8, 1, 0.8, alpha));

    G4Colour::AddToMap("scattering_chamber_vacuum", G4Colour(0.53,0.80,0.98));
    
    G4Colour::AddToMap("scifi", G4Colour(0.99,0.72,0.12));
    G4Colour::AddToMap("scifi_cladding", G4Colour(0.75,0.17,0.10));
    G4Colour::AddToMap("scwall", G4Colour(0.99,0.72,0.12));
    //    G4Colour::AddToMap("bcc", G4Colour(0.52,0.84,1.00));
    G4Colour::AddToMap("bcc", G4Colour(0.99,0.72,0.12));
    //    G4Colour::AddToMap("bm", G4Colour(0.52,0.84,1.00));
    G4Colour::AddToMap("bm", G4Colour(0.99,0.72,0.12));
    G4Colour::AddToMap("floor", G4Colour(0.95,0.95,0.95));
    G4Colour::AddToMap("testplane", G4Colour(1,0,0));
    
    G4Colour::AddToMap("lead", G4Colour(0.40,0.40,0.44));
    G4Colour::AddToMap("aluminum", G4Colour(0.80,0.80,0.80));
    G4Colour::AddToMap("table", G4Colour(0.80,0.80,0.80,0.10));
    G4Colour::AddToMap("concrete", G4Colour(0.70,0.70,0.70));
    G4Colour::AddToMap("plastic", G4Colour(0.96,0.88,0.69));
    G4Colour::AddToMap("pmt", G4Colour(0.50,0.50,0.30));
    G4Colour::AddToMap("tedlar", G4Colour(0.1,0.1,0.1));
    G4Colour::AddToMap("carbon", G4Colour(0.1,0.1,0.1));
    G4Colour::AddToMap("kapton", G4Colour(0.75,0.20,0.02));
    G4Colour::AddToMap("copper", G4Colour(0.65,0.40,0.31));
    
    
    expHall_x_ = 0*m;
    expHall_y_ = 0*m;
    expHall_z_ = 0*m;
    experimentalHall_log_ = NULL;
    experimentalHall_phys_ = NULL;
    
    /// \todo look into
    StepLimit_ = new G4UserLimits(0.5*mm);
    
    /// \todo look into
    detectorMessenger_ = new g4PSIDetectorMessenger();
}

g4PSIDetectorConstruction::~g4PSIDetectorConstruction() {
    
    // need to delete all the other objects
    
    delete StepLimit_;
}

G4Colour g4PSIDetectorConstruction::G4HexColour(int red, int green, int blue) {
    return G4Colour( (double) red/256.,
                    (double) green/256.,
                    (double) blue/256. );
}


void g4PSIDetectorConstruction::Place_Target(G4double angle, G4bool in) {
    
    G4double tgt_z = 20*mm;
    G4double tgt_x = 16*cm;
    G4double tgt_y = 20*cm;
    
    G4double support_z = 21*cm;
    G4double support_x = 18.5*cm;
    G4double support_y = 20*mm;
    
    G4Material *tgt_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS");
    
    G4Box* support_geom = new G4Box("tgt_support_geom", support_x/2., support_y/2., support_z/2.);
    G4LogicalVolume* support_log = new G4LogicalVolume(support_geom, tgt_mat, "tgt_support_log", 0, 0, 0);
    
    G4RotationMatrix* yRot = new G4RotationMatrix;
    /// \todo Why rotating by the cos/sin of an angle?
    yRot->rotateY(support_z/sqrt(support_z*support_z + support_x*support_x));
    if (angle > 0) {
        new G4PVPlacement(yRot, G4ThreeVector(0,-tgt_y/2-support_y/2,0),
                          support_log, "tgt_support", experimentalHall_log_, false, 0);
    }
    
    if (in) {
        G4Box* tgt_geom = new G4Box("tgt_geom", tgt_x/2., tgt_y/2., tgt_z/2.);
        G4LogicalVolume* tgt_log = new G4LogicalVolume(tgt_geom, tgt_mat, "tgt_log", 0, 0, 0);
        
        G4RotationMatrix* yRotTgt = new G4RotationMatrix;
        yRotTgt->rotateY(-angle);
        new G4PVPlacement(yRotTgt, G4ThreeVector(0,0,0),
                          tgt_log, "target", experimentalHall_log_, false, 0);
    }
}


void g4PSIDetectorConstruction::Place_SciFi_Detectors(double z) {
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    if (DetectorParts->build[g4PSIDetectorParts::kBeamSciFi] ||
        DetectorParts->build[g4PSIDetectorParts::kBeamSciFi_Type0]) {
        // original detector with circular fibers
        double f = 2.0 * mm;
        double cl = 0.06 * mm;
        DetectorParts->Add(new g4PSIBeamSC("TSFX", 90*deg, f, cl, 0, 8*cm, 8*cm, 0, 0, z - f));
        DetectorParts->Add(new g4PSIBeamSC("TSFY",  0*deg, f, cl, 0, 8*cm, 8*cm, 0, 0, z    ));
        DetectorParts->Add(new g4PSIBeamSC("TSFU", 45*deg, f, cl, 0, 8*cm, 8*cm, 0, 0, z + f));
    } else if (DetectorParts->build[g4PSIDetectorParts::kBeamSciFi_Type1]) {
        // original detector with square fibers
        double f = 2.0 * mm;
        double cl = 0.06 * mm;
        DetectorParts->Add(new g4PSIBeamSC("TSFX", 90*deg, -f, cl, 0, 8*cm, 8*cm, 0, 0, z - f));
        DetectorParts->Add(new g4PSIBeamSC("TSFY",  0*deg, -f, cl, 0, 8*cm, 8*cm, 0, 0, z    ));
        DetectorParts->Add(new g4PSIBeamSC("TSFU", 45*deg, -f, cl, 0, 8*cm, 8*cm, 0, 0, z + f));
    } else if (DetectorParts->build[g4PSIDetectorParts::kBeamSciFi_Type2]) {
        // two double planes with gaps between fibers
        double f = 2.0 * mm;
        double cl = 0.06 * mm;
        double d = 1.4 * mm;
        double phi = 45 * deg;
        double ox = (f + d) / 2. * cos(phi);
        double oy = -(f + d) / 2. * sin(phi);
        DetectorParts->Add(new g4PSIBeamSC("TSFU1", phi, f, cl, d, 8*cm, 8*cm,  0, 0, z - 2*f));
        DetectorParts->Add(new g4PSIBeamSC("TSFU2", phi, f, cl, d, 8*cm, 8*cm, ox, oy, z - f));
        phi -= 90 * deg;
        ox = (f + d) / 2. * cos(phi);
        oy = -(f + d) / 2. * sin(phi);
        DetectorParts->Add(new g4PSIBeamSC("TSFV1", phi, f, cl, d, 8*cm, 8*cm,  0, 0, z + f));
        DetectorParts->Add(new g4PSIBeamSC("TSFV2", phi, f, cl, d, 8*cm, 8*cm, ox, oy, z + 2*f));
    } else if (DetectorParts->build[g4PSIDetectorParts::kBeamSciFi_Type3]) {
        // two double planes without gaps
        double f = 2.0 * mm;
        double cl = 0.06 * mm;
        double d = 0.0 * mm;
        double phi = 45 * deg;
        double ox = (f + d) / 2. * cos(phi);
        double oy = -(f + d) / 2. * sin(phi);
        DetectorParts->Add(new g4PSIBeamSC("TSFU1", phi, f, cl, d, 8*cm, 8*cm,  0,  0, z - 2*f));
        DetectorParts->Add(new g4PSIBeamSC("TSFU2", phi, f, cl, d, 8*cm, 8*cm, ox, oy, z - f));
        phi = -45 * deg;
        ox = (f + d) / 2. * cos(phi);
        oy = -(f + d) / 2. * sin(phi);
        DetectorParts->Add(new g4PSIBeamSC("TSFV1", phi, f, cl, d, 8*cm, 8*cm,  0,  0, z + f));
        DetectorParts->Add(new g4PSIBeamSC("TSFV2", phi, f, cl, d, 8*cm, 8*cm, ox, oy, z + 2*f));
    }
}


void g4PSIDetectorConstruction::Place_SC(G4String name,
                                         G4double r, // to the center of the SC
                                         G4double angle) {
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    
    G4double sc_wid = 5*cm;
    G4double sc_len = 50*cm;
    
    G4double rail_a = 3*cm;
    
    G4double scA_x = 50*cm;
    G4double scA_y = sc_len;
    G4double scA_z = sc_wid + rail_a;
    
    G4double del = scA_z/2. - sc_wid/2.;
    G4double r0 = r - del;
    
    G4Box* scA_box = new G4Box("scA_solid", scA_x/2., scA_y/2., scA_z/2.);
    
    G4LogicalVolume* scA_log = new G4LogicalVolume(scA_box, Air, name+"_log");
    G4RotationMatrix* RotSC = new G4RotationMatrix;
    RotSC->rotateY(-angle);
    new G4PVPlacement(RotSC, G4ThreeVector(r0*sin(angle), 0., r0*cos(angle)),
                      scA_log, name, experimentalHall_log_, false, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    VisAtt->SetForceWireframe(true);
    scA_log->SetVisAttributes(VisAtt);
    scA_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    
    G4Box* rail_22_box = new G4Box("al_22_solid", 22.*cm/2., rail_a/2., rail_a/2.);
    G4LogicalVolume* rail_22_log = new G4LogicalVolume(rail_22_box, Al, "al_22_log");
    new G4PVPlacement(0, G4ThreeVector(+5.5*cm, 15*cm, del-sc_wid/2. - rail_a/2.),
                      rail_22_log,"al_22up_phys", scA_log, false, 0);
    G4Box* rail_40_box = new G4Box("al_40_solid", 40.*cm/2., rail_a/2., rail_a/2.);
    G4LogicalVolume* rail_40_log = new G4LogicalVolume(rail_40_box, Al, "al_40_log");
    new G4PVPlacement(0, G4ThreeVector(+5.5*cm, -15*cm, del-sc_wid/2. - rail_a/2.),
                      rail_40_log,"al_40down_phys", scA_log, false, 0);
    G4Box* rail_30_box = new G4Box("al_30_solid", rail_a/2., 30*cm/2., rail_a/2.);
    G4LogicalVolume* rail_30_log = new G4LogicalVolume(rail_30_box, Al, "al_30_log");
    new G4PVPlacement(0, G4ThreeVector(-7*cm, 1.5*cm, del-sc_wid/2. - rail_a/2.),
                      rail_30_log,"al_30_phys", scA_log, false, 0);
    
    G4Colour col;
    G4Colour::GetColour("aluminum", col);
    G4VisAttributes* concreteVisAtt = new G4VisAttributes(col);
    rail_22_log->SetVisAttributes(concreteVisAtt);
    rail_40_log->SetVisAttributes(concreteVisAtt);
    rail_30_log->SetVisAttributes(concreteVisAtt);
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    DetectorParts->Add(scA_log, new g4PSISCWall("SC", sc_wid, sc_len, sc_wid, 0*cm, del));
}


void g4PSIDetectorConstruction::Place_GEMtelescope(G4String name,
                                                   G4double r, //
                                                   G4double z0, //  distance to front-face upstream GEM
                                                   G4double angle,
                                                   G4double r1, G4double r2, G4double r3, G4double rsc) {
    
    /// \todo double check dimensions and overlapping/spillover material
    
    G4double rail_a = 3*cm;
    G4double rail_l = 90*cm;
    G4double plate_z = 1.5*cm;
    G4double y0 = -(12.5*cm/2 + 6.4*cm + rail_a/2.);
    G4double x0 = 5*cm;
    G4double dGEM = -0.5*cm;  // front face of plate to front face of 1st GEM
    
    double telescope_x = 50.*cm;
    double telescope_y = 50.*cm;
    double telescope_z = dGEM < 0 ? rail_l + plate_z - dGEM : rail_l + plate_z;
    
    G4double dz = dGEM > 0 ? r - dGEM + telescope_z / 2. : r + telescope_z/2.;
    
    G4Box* telescope_box
    = new G4Box("al_solid", telescope_x/2., telescope_y/2., telescope_z/2.);
    G4LogicalVolume* telescope_log =
    new G4LogicalVolume(telescope_box,
                        G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), name+"_log");
    G4RotationMatrix* RotTele = new G4RotationMatrix;
    RotTele->rotateY(-angle);
    new G4PVPlacement(RotTele, G4ThreeVector(dz*sin(angle), 0., z0+dz*cos(angle)),
                      telescope_log, name, experimentalHall_log_, false, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    VisAtt->SetForceWireframe(true);
    telescope_log->SetVisAttributes(VisAtt);
    
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    
    // telescope support structure
    
    G4Box* rail_box
    = new G4Box("al_solid", rail_a/2., rail_a/2., rail_l/2.);
    G4LogicalVolume* rail_log = new G4LogicalVolume(rail_box, Al, "al_log");
    new G4PVPlacement(0, G4ThreeVector(+x0, y0, plate_z/2),
                      rail_log,"al_rail1_r", telescope_log, false, 0);
    new G4PVPlacement(0, G4ThreeVector(-x0, y0, plate_z/2),
                      rail_log,"al_rail2_r", telescope_log, false, 0);
    
    G4Box* plate_box = new G4Box("plate_solid", 13.*cm/2., 6.*cm/2., plate_z/2.);
    G4LogicalVolume* plate_log = new G4LogicalVolume(plate_box, Al, "plate_log");
    new G4PVPlacement(0, G4ThreeVector(0., y0 + 1.5*cm, -telescope_z/2+plate_z/2),
                      plate_log,"plate_r_phys", telescope_log, false, 0);
    
    G4Colour col;
    G4Colour::GetColour("aluminum", col);
    G4VisAttributes* floorVisAtt = new G4VisAttributes(col);
    rail_log->SetVisAttributes(floorVisAtt);
    plate_log->SetVisAttributes(floorVisAtt);
    
    //  telescope detectors
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    DetectorParts->Add(telescope_log, new g4PSIGEM(name+"GEM1", 10.0*cm, 12.5*cm, r1+dGEM-telescope_z/2.));
    DetectorParts->Add(telescope_log, new g4PSIGEM(name+"GEM2", 10.0*cm, 12.5*cm, r2+dGEM-telescope_z/2.));
    DetectorParts->Add(telescope_log, new g4PSIGEM(name+"GEM3", 10.0*cm, 12.5*cm, r3+dGEM-telescope_z/2.));
    if (rsc > 0) {
        DetectorParts->Add(telescope_log, new g4PSISCWall(name+"SC", 12*cm, 12*cm, 2*mm, 0*cm, rsc+dGEM-telescope_z/2.));
    }
}


void g4PSIDetectorConstruction::Place_TableAndFrames() {
    // Material
    
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    
    // ------------------Sliding table
    
    //  rail top
    const G4double railtop_a = 50*cm;
    const G4double railtop_b = 1*cm;
    const G4double railtop_l = 80*cm;
    const G4double x0_top = 0*cm;
    const G4double y0_top = -14*cm;
    
    const G4double z0_top = -60*cm;
    
    G4Box* railtop_box = new G4Box("al_railtop_solid", railtop_a/2., railtop_b/2., railtop_l/2.);
    G4LogicalVolume* railtop_log = new G4LogicalVolume(railtop_box, Al, "al_railtop_log");
    new G4PVPlacement(0, G4ThreeVector(x0_top, y0_top, z0_top), railtop_log,"al_railtop", experimentalHall_log_, false, 0);
    
    // rail
    
    const G4double rail_a = 3*cm;
    const G4double y0 = y0_top - railtop_b/2. - rail_a/2.;
    const G4double x0 = 15*cm;
//    const G4double z0_upstream = -27*cm;
//    const G4double z0_downstream = -55*cm;
    // const G4double rail_l = z0_upstream - z0_downstream;
    const G4double rail_l = 80*cm;
    
    G4Box* rail_box = new G4Box("al_rail_solid", rail_a/2., rail_a/2., rail_l/2.);
    G4LogicalVolume* rail_log = new G4LogicalVolume(rail_box, Al, "al_rail_log");
    //  new G4PVPlacement(0, G4ThreeVector(+x0, y0, (z0_upstream+z0_downstream)/2.), rail_log,"al_rail1_r", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(+x0, y0, z0_top - 46*cm), rail_log,"al_rail1_r", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(-x0, y0, z0_top - 46*cm), rail_log,"al_rail2_r", experimentalHall_log_, false, 0);
    
    // Rail bottom
    
    const G4double railbottom_a = 50*cm;
    const G4double railbottom_b = 1*cm;
    const G4double railbottom_l = rail_l;
    const G4double x0_bottom = 0*cm;
    const G4double y0_bottom = y0 - rail_a/2. - railbottom_b/2.;
    
    G4Box* railbottom_box = new G4Box("al_railb_solid", railbottom_a/2., railbottom_b/2., railbottom_l/2.);
    G4LogicalVolume* railbottom_log = new G4LogicalVolume(railbottom_box, Al, "al_railb_log");
    new G4PVPlacement(0, G4ThreeVector(x0_bottom, y0_bottom, z0_top - 46*cm), railbottom_log,"al_railbottom", experimentalHall_log_, false, 0);
    
    // Rail support
    
    const G4double railsupport_a = 3*cm;
    const G4double railsupport_b = 53*cm;
    const G4double railsupport_l = 3*cm;
    const G4double x0_support = 15*cm;
    const G4double y0_support = y0_bottom - railbottom_b/2. - railsupport_b/2.;
    
    G4Box* railsupport_box = new G4Box("al_rails_solid", railsupport_a/2., railsupport_b/2., railsupport_l/2.);
    G4LogicalVolume* railsupport_log = new G4LogicalVolume(railsupport_box, Al, "al_rails_log");
    new G4PVPlacement(0, G4ThreeVector(x0_support, y0_support, -82*cm), railsupport_log,"al_rails1_r", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(-x0_support, y0_support, -82*cm), railsupport_log,"al_rails2_r", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, y0_support, -130*cm), railsupport_log,"al_rails2_r", experimentalHall_log_, false, 0);
    
    
    
    const G4double railsupport2_a = 30*cm;
    const G4double railsupport2_b = 3*cm;
    const G4double railsupport2_l = 3*cm;
    
    const G4double y0_support2 = y0_support - railsupport_b/2. - railsupport2_b/2.;
    
    G4Box* railsupport2_box = new G4Box("al_rails2_solid", railsupport2_a/2., railsupport2_b/2., railsupport2_l/2.);
    G4LogicalVolume* railsupport2_log = new G4LogicalVolume(railsupport2_box, Al, "al_rails2_log");
    new G4PVPlacement(0, G4ThreeVector(0., y0_support2, -82*cm), railsupport2_log,"al_rails", experimentalHall_log_, false, 0);
    
    // ------------------- Rotating system
    
    // circle table
    
    const G4double circle_table_radius = 44*cm;
    const G4double circle_table_inner_radius = 18*cm;
    const G4double circle_table_height = 1*cm;
    const G4double x0_circletable = 0*cm;
    const G4double y0_circletable = -63.5*cm;
    const G4double z0_circletable = 0*cm;
    G4RotationMatrix* Rot= new G4RotationMatrix;
    Rot->rotateX(90*deg);
    
    G4Tubs* circle_table = new G4Tubs("circle_table", circle_table_inner_radius, circle_table_radius, circle_table_height/2.0, 0, twopi);
    
    G4LogicalVolume *circle_table_log = new G4LogicalVolume(circle_table, Al, "al_circletable_log");
    new G4PVPlacement(Rot, G4ThreeVector(x0_circletable, y0_circletable, z0_circletable), circle_table_log,"al_circle_table", experimentalHall_log_, false, 0);
    
    
    
    const G4double circle_table2_radius = 40*cm;
    const G4double circle_table2_inner_radius = 26*cm;
    const G4double circle_table2_height = 4*cm;
    const G4double x0_circletable2 = 0*cm;
    const G4double y0_circletable2 = y0_circletable - circle_table_height/2. - 1*cm - circle_table2_height/2.;
    const G4double z0_circletable2 = 0*cm;
    
    G4Tubs* circle_table2 = new G4Tubs("circle_table2", circle_table2_inner_radius, circle_table2_radius, circle_table2_height/2.0, 0, twopi);
    
    G4LogicalVolume *circle_table2_log = new G4LogicalVolume(circle_table2, Al, "al_circletable2_log");
    new G4PVPlacement(Rot, G4ThreeVector(x0_circletable2, y0_circletable2, z0_circletable2), circle_table2_log,"al_circle_table2", experimentalHall_log_, false, 0);
    
    
    
    const G4double circle_table3_radius = 50*cm;
    const G4double circle_table3_inner_radius = 25*cm;
    const G4double circle_table3_height = 2.5*cm;
    const G4double x0_circletable3 = 0*cm;
    const G4double y0_circletable3 = y0_circletable2 - circle_table2_height/2. - circle_table3_height/2.;
    const G4double z0_circletable3 = 0*cm;
    
    G4Tubs* circle_table3 = new G4Tubs("circle_table3", circle_table3_inner_radius, circle_table3_radius, circle_table3_height/2.0, 0, twopi);
    
    G4LogicalVolume *circle_table3_log = new G4LogicalVolume(circle_table3, Al, "al_circletable3_log");
    new G4PVPlacement(Rot, G4ThreeVector(x0_circletable3, y0_circletable3, z0_circletable3), circle_table3_log,"al_circle_table3", experimentalHall_log_, false, 0);
    
    
    // triangle
    
    G4double theta = 73.5*deg;
    G4RotationMatrix* zRot1 = new G4RotationMatrix;
    G4RotationMatrix* zRot2 = new G4RotationMatrix;
    zRot1->rotateY(theta);
    zRot2->rotateY(-theta);
    
    const G4double tri_bar_a = 92*cm;
    const G4double tri_bar2_a = 160*cm;
    const G4double tri_bar_b = 4*cm;
    const G4double tri_bar_l = 4*cm;
    //    const G4double distance_to_front = 30*cm;
    const G4double x0_tri1 = (23*cm) * sin(theta);
    //    const G4double x0_tri2 = (distance_to_front + tri_bar_l/2.) * sin(-theta);
    //    const G4double y0_tri = y0_circletable - circle_table_height/2. - tri_bar_b/2.;
    const G4double z0_tri1 = (23*cm) * cos(theta) - 65*cm;
    
    
    const G4double x0_tri = 0*cm;
    const G4double y0_tri = y0_circletable3 - circle_table3_height/2. - tri_bar_b/2.;
    const G4double z0_tri = 20*cm;
    
    G4Box* tri_box = new G4Box("al_tri_solid", tri_bar_a/2., tri_bar_b/2., tri_bar_l/2.);
    
    
    G4LogicalVolume* tri_log = new G4LogicalVolume(tri_box, Al, "al_tri_log");
    new G4PVPlacement(0, G4ThreeVector(x0_tri, y0_tri, z0_tri), tri_log,"al_tri1", experimentalHall_log_, false, 0);
    //    new G4PVPlacement(zRot2, G4ThreeVector(x0_tri2, y0_tri, z0_tri), tri_log,"al_tri2", experimentalHall_log_, false, 0);
    //    new G4PVPlacement(0, G4ThreeVector(0, y0_tri, -(distance_to_front + tri_bar_l/2.)), tri_log,"al_tri3", experimentalHall_log_, false, 0);
    
    G4Box* tri2_box = new G4Box("al_tri2_solid", tri_bar2_a/2., tri_bar_b/2., tri_bar_l/2.);
    G4LogicalVolume* tri2_log = new G4LogicalVolume(tri2_box, Al, "al_tri2_log");
    new G4PVPlacement(zRot1, G4ThreeVector(x0_tri1, y0_tri, z0_tri1), tri2_log,"al_tri2", experimentalHall_log_, false, 0);
    new G4PVPlacement(zRot2, G4ThreeVector(-x0_tri1, y0_tri, z0_tri1), tri2_log,"al_tri3", experimentalHall_log_, false, 0);
    
    const G4double tri_bar3_a = 50*cm;
    G4Box* tri3_box = new G4Box("al_tri3_solid", tri_bar3_a/2., tri_bar_b/2., tri_bar_l/2.);
    G4LogicalVolume* tri3_log = new G4LogicalVolume(tri3_box, Al, "al_tri3_log");
    new G4PVPlacement(0, G4ThreeVector(0., y0_tri, -40*cm), tri3_log,"al_tri", experimentalHall_log_, false, 0);
    
    // triangle feet
    
    const G4double tri_feet_a = 4.5*cm;
    const G4double tri_feet_b = 55*cm;
    const G4double tri_feet_l = 4.5*cm;
    
    const G4double x0_tf1 = 43*cm;
    const G4double x0_tf2 = -x0_tf1;
    const G4double x0_tf3 = 0*cm;
    const G4double y0_tf = y0_tri - tri_bar_b/2. - (tri_feet_b/2.);
    const G4double z0_tf1 = z0_tri - 2*cm;
    const G4double z0_tf2 = -126*cm;
    
    G4Box* tf_box = new G4Box("al_tf_solid", tri_feet_a/2., tri_feet_b/2., tri_feet_l/2.);
    G4LogicalVolume* tf_log = new G4LogicalVolume(tf_box, Al, "al_tf_log");
    new G4PVPlacement(0, G4ThreeVector(x0_tf1, y0_tf, z0_tf1), tf_log,"al_tf1", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(x0_tf2, y0_tf, z0_tf1), tf_log,"al_tf2", experimentalHall_log_, false, 0);
    new G4PVPlacement(0, G4ThreeVector(x0_tf3, y0_tf, z0_tf2), tf_log,"al_tf3", experimentalHall_log_, false, 0);
    
    
    
    // table feet support
    
    const G4double table_fs_a = 61*cm;
    const G4double table_fs_b = 3*cm;
    const G4double table_fs_l = 3*cm;
    
    G4double tt = 52.4*deg;
    G4RotationMatrix* rot1 = new G4RotationMatrix;
    G4RotationMatrix* rot2 = new G4RotationMatrix;
    rot1->rotateZ(tt);
    rot2->rotateZ(-tt);
    
    const G4double x0_tfs = table_fs_a * cos(tt)/2. + 4*cm;
    
    const G4double y0_tfs = y0_tri - tri_bar_b/2. - table_fs_a * sin(tt)/2.;
    const G4double z0_tfs = z0_tri;
    
    
    G4Box* tfs_box = new G4Box("al_tfs_solid", table_fs_a/2., table_fs_b/2., table_fs_l/2.);
    G4LogicalVolume* tfs_log = new G4LogicalVolume(tfs_box, Al, "al_tfs_log");
    new G4PVPlacement(rot1, G4ThreeVector(x0_tfs, y0_tfs, z0_tfs), tfs_log,"al_tfs1", experimentalHall_log_, false, 0);
    new G4PVPlacement(rot2, G4ThreeVector(-x0_tfs, y0_tfs, z0_tfs), tfs_log,"al_tfs2", experimentalHall_log_, false, 0);
    //    new G4PVPlacement(0, G4ThreeVector(x0_tf2, y0_tf, z0_tf1), tf_log,"al_tf2", experimentalHall_log_, false, 0);
    //    new G4PVPlacement(0, G4ThreeVector(x0_tf1, y0_tf, z0_tf2), tf_log,"al_tf3", experimentalHall_log_, false, 0);
    //    new G4PVPlacement(0, G4ThreeVector(x0_tf2, y0_tf, z0_tf2), tf_log,"al_tf4", experimentalHall_log_, false, 0);
    // -------------------------------------------------------------------------------------------------------------------
    const G4double table_fs1_a = 85*cm;
    G4double tt1 = 37*deg;
    G4RotationMatrix* rot3 = new G4RotationMatrix;
    G4RotationMatrix* rot4 = new G4RotationMatrix;
    rot3->rotateY(theta);
    rot3->rotateZ(tt1);
    rot4->rotateY(-theta);
    rot4->rotateZ(-tt1);
    
    const G4double x0_tfs1 = 34*cm;
    
    const G4double y0_tfs1 = y0_tfs;
    const G4double z0_tfs1 = -16*cm;
    G4Box* tfs1_box = new G4Box("al_tfs1_solid", table_fs1_a/2., table_fs_b/2., table_fs_l/2.);
    G4LogicalVolume* tfs1_log = new G4LogicalVolume(tfs1_box, Al, "al_tfs1_log");
    new G4PVPlacement(rot3, G4ThreeVector(x0_tfs1, y0_tfs1, z0_tfs1), tfs1_log,"al_tfs3", experimentalHall_log_, false, 0);
    new G4PVPlacement(rot4, G4ThreeVector(-x0_tfs1, y0_tfs1, z0_tfs1), tfs1_log,"al_tfs4", experimentalHall_log_, false, 0);
    // -------------------------------------------------------------------------------------------------------------------
    
    const G4double table_fs2_a = 94*cm;
    G4double tt2 = 35*deg;
    G4RotationMatrix* rot5 = new G4RotationMatrix;
    G4RotationMatrix* rot6 = new G4RotationMatrix;
    rot5->rotateY(theta);
    rot5->rotateZ(-tt2);
    rot6->rotateY(-theta);
    rot6->rotateZ(tt2);
    
    const G4double x0_tfs2 = 13*cm;
    
    const G4double y0_tfs2 = y0_tfs;
    const G4double z0_tfs2 = -87*cm;
    G4Box* tfs2_box = new G4Box("al_tfs2_solid", table_fs2_a/2., table_fs_b/2., table_fs_l/2.);
    G4LogicalVolume* tfs2_log = new G4LogicalVolume(tfs2_box, Al, "al_tfs2_log");
    new G4PVPlacement(rot5, G4ThreeVector(x0_tfs2, y0_tfs2, z0_tfs2), tfs2_log,"al_tfs5", experimentalHall_log_, false, 0);
    new G4PVPlacement(rot6, G4ThreeVector(-x0_tfs2, y0_tfs2, z0_tfs2), tfs2_log,"al_tfs6", experimentalHall_log_, false, 0);
    
    G4Colour col;
    G4Colour::GetColour("table", col);
    G4VisAttributes* alVisAttr = new G4VisAttributes(col);
    railtop_log->SetVisAttributes(alVisAttr);
    rail_log->SetVisAttributes(alVisAttr);
    railbottom_log->SetVisAttributes(alVisAttr);
    railsupport_log->SetVisAttributes(alVisAttr);
    
    circle_table_log->SetVisAttributes(alVisAttr);
    circle_table2_log->SetVisAttributes(alVisAttr);
    circle_table3_log->SetVisAttributes(alVisAttr);
    tri_log->SetVisAttributes(alVisAttr);
    tri2_log->SetVisAttributes(alVisAttr);
    tf_log->SetVisAttributes(alVisAttr);
    tfs_log->SetVisAttributes(alVisAttr);
    tfs1_log->SetVisAttributes(alVisAttr);
    tfs2_log->SetVisAttributes(alVisAttr);

}


void g4PSIDetectorConstruction::Place_Shielding() {
    
    G4Material* Lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
    G4Material* concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
    
    G4double floor_y = -150*cm;
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    
    // Concrete Floor
    if (DetectorParts->build[g4PSIDetectorParts::kShieldFloor]) {
        G4double floor_height = 150*cm;
        G4Box* floor_box
        = new G4Box("floor_solid",expHall_x_, floor_height/2., expHall_z_);
        G4LogicalVolume* floor_log = new G4LogicalVolume(floor_box,
                                                         concrete, "floor_log");
        new G4PVPlacement(0, G4ThreeVector(0., floor_y - floor_height/2., 0.),
                          floor_log,"floor", experimentalHall_log_, false, 0);
        G4Colour colour_floor;
        G4Colour::GetColour("floor", colour_floor);
        G4VisAttributes* floorVisAtt = new G4VisAttributes(colour_floor);
        floor_log->SetVisAttributes(floorVisAtt);
    }
    
    // Concrete Shielding tube in the area of the GEMs
    if (DetectorParts->build[g4PSIDetectorParts::kShieldBLConcreteTube]) {
        G4double shieldTubeR1 = 30.0*cm;
        G4double shieldTubeR2 = 40.0*cm;
        G4double shieldTubeL = 60.0*cm;
        G4double shieldTubeEndZ = -20.0*cm;
        G4Tubs* shield_tube
        = new G4Tubs("shieldTube_solid", shieldTubeR1, shieldTubeR2, shieldTubeL/2., 0, twopi);
        G4LogicalVolume* shieldTube_log = new G4LogicalVolume(shield_tube,
                                                              concrete, "shieldTube_log");
        //  G4double = -floor_height/2 + shieldTube_z/2
        G4double shieldTube_zPos = -shieldTubeL/2 + shieldTubeEndZ;
        new G4PVPlacement(0, G4ThreeVector( 0, 0, shieldTube_zPos),
                          shieldTube_log,"shieldTube", experimentalHall_log_, false, 0);
        G4Colour colour_concrete;
        G4Colour::GetColour("concrete", colour_concrete);
        G4VisAttributes* shieldVisAtt = new G4VisAttributes(colour_concrete);
        shieldTube_log->SetVisAttributes(shieldVisAtt);
    }
    
    // Concrete Shielding Blocks around the beamline
    if (DetectorParts->build[g4PSIDetectorParts::kShieldBLConcrete]) {
        G4double shieldwall_halfx = 150*cm;
        G4double shieldwall_halfy = 150*cm;
        G4double shieldwall_halfz = 50.0*cm;
        G4double shieldwall_offset = (shieldwall_halfx/cm + 20.32 + 2.0)*cm;
        G4ThreeVector shieldwall_position_1 = G4ThreeVector(shieldwall_offset,0,-(140.0 + shieldwall_halfz/cm)*cm);
        G4ThreeVector shieldwall_position_2 = G4ThreeVector(-shieldwall_offset,0,-(140.0 + shieldwall_halfz/cm)*cm);
        
        G4Box* shieldwall_solid = new G4Box("ShieldWall",shieldwall_halfx,shieldwall_halfy,shieldwall_halfz);
        G4LogicalVolume* shieldwall_logical = new G4LogicalVolume(shieldwall_solid,concrete,"ShieldWall",0,0,0);
        new G4PVPlacement(0, shieldwall_position_1, shieldwall_logical, "ShieldWallL", experimentalHall_log_, false, 0);
        new G4PVPlacement(0, shieldwall_position_2, shieldwall_logical, "ShieldWallR", experimentalHall_log_, false, 1);
        
        /*    G4double shieldBlock_x = 100*cm;
         G4double shieldBlock_y = 100*cm;
         G4double shieldBlock_z = 50*cm;
         G4Box* shieldBlock_box
         = new G4Box("shieldBlock_solid",shieldBlock_x/2,shieldBlock_y/2,shieldBlock_z/2);
         G4LogicalVolume* shieldBlock_log = new G4LogicalVolume(shieldBlock_box,
         concrete, "shieldBlock_log");
         //  G4double = -floor_height/2 + shieldBlock_z/2
         G4double shieldBlock_zPos = -shieldBlock_z/2 - 95*cm;
         new G4PVPlacement(0, G4ThreeVector( shieldBlock_x/2 + 5*cm, floor_y + (1./2.)*shieldBlock_y, shieldBlock_zPos),
         shieldBlock_log,"shieldBlock0", experimentalHall_log_, false, 0);
         new G4PVPlacement(0, G4ThreeVector( shieldBlock_x/2 + 5*cm, floor_y + (3./2.)*shieldBlock_y, shieldBlock_zPos),
         shieldBlock_log,"shieldBlock0", experimentalHall_log_, false, 0);
         new G4PVPlacement(0, G4ThreeVector(-shieldBlock_x/2 - 5*cm, floor_y + (1./2.)*shieldBlock_y, shieldBlock_zPos),
         shieldBlock_log,"shieldBlock1", experimentalHall_log_, false, 0);
         new G4PVPlacement(0, G4ThreeVector(-shieldBlock_x/2 - 5*cm, floor_y + (3./2.)*shieldBlock_y, shieldBlock_zPos),
         shieldBlock_log,"shieldBlock0", experimentalHall_log_, false, 0);
         G4Colour colour_concrete;
         G4Colour::GetColour("concrete", colour_concrete);
         G4VisAttributes* shieldVisAtt = new G4VisAttributes(colour_concrete);
         shieldBlock_log->SetVisAttributes(shieldVisAtt);
         */
    }
    
    // Downstream Concrete Shieldwall
    if (DetectorParts->build[g4PSIDetectorParts::kShieldWall]) {
        G4double backWall_x = 400*cm;
        G4double backWall_y = 300*cm;
        G4double backWall_z = 100*cm;
        G4double backWall_zPos = 400*cm;
        
        G4Box* backWall_box
        = new G4Box("backWall_solid",backWall_x/2,backWall_y/2,backWall_z/2);
        G4LogicalVolume* backWall_log = new G4LogicalVolume(backWall_box,
                                                            concrete, "backWall_log");
        new G4PVPlacement(0, G4ThreeVector(0., floor_y + backWall_y / 2, backWall_zPos + backWall_z /2),
                          backWall_log,"backWall0", experimentalHall_log_, false, 0);
        G4Colour colour_concrete;
        G4Colour::GetColour("concrete", colour_concrete);
        G4VisAttributes* backWallVisAtt = new G4VisAttributes(colour_concrete);
        backWall_log->SetVisAttributes(backWallVisAtt);
    }
    
    // Lead Bricks around the GEMs
    if (DetectorParts->build[g4PSIDetectorParts::kShieldBLLead]) {
        G4double inch = 2.54*cm;
        G4double leadBrick_x = 2*inch;
        G4double leadBrick_y = 8*inch;
        G4double leadBrick_z = 4*inch;
        G4RotationMatrix* RotPb = new G4RotationMatrix;
        RotPb->rotateZ(90*deg);
        G4RotationMatrix* RotPbR = NULL; //new G4RotationMatrix;
        G4RotationMatrix* RotPbL = NULL; // new G4RotationMatrix;
        
        G4Box* leadBrick_box
        = new G4Box("leadBrick_solid",leadBrick_x/2,leadBrick_y/2,leadBrick_z/2);
        G4LogicalVolume* leadBrick_log = new G4LogicalVolume(leadBrick_box,
                                                             Lead, "leadBrick_log");
        G4double zPb = -20*cm;
        for (G4int i = 0; i < 3; i++) {
            new G4PVPlacement(RotPb, G4ThreeVector(0*cm, leadBrick_y/2 + leadBrick_x/2, zPb),
                              leadBrick_log,"leadBrickA", experimentalHall_log_, false, 0);
            new G4PVPlacement(RotPb, G4ThreeVector(0*cm, -(leadBrick_y/2 + leadBrick_x/2), zPb),
                              leadBrick_log,"leadBrickB", experimentalHall_log_, false, 0);
            new G4PVPlacement(RotPbR, G4ThreeVector(5.*cm, 0., zPb),
                              leadBrick_log,"leadBrickC", experimentalHall_log_, false, 0);
            new G4PVPlacement(RotPbL, G4ThreeVector(-5.*cm, 0., zPb),
                              leadBrick_log,"leadBrickD", experimentalHall_log_, false, 0);
            zPb -= 20*cm;
        };
        
        new G4PVPlacement(RotPb, G4ThreeVector(0*cm, leadBrick_y/2 + leadBrick_x/2, zPb),
                          leadBrick_log,"leadBrickA", experimentalHall_log_, false, 0);
        new G4PVPlacement(RotPb, G4ThreeVector(0*cm, -(leadBrick_y/2 + leadBrick_x/2), zPb),
                          leadBrick_log,"leadBrickB", experimentalHall_log_, false, 0);
        new G4PVPlacement(RotPbR, G4ThreeVector(8.*cm, 0., zPb),
                          leadBrick_log,"leadBrickC", experimentalHall_log_, false, 0);
        new G4PVPlacement(RotPbL, G4ThreeVector(-8.*cm, 0., zPb),
                          leadBrick_log,"leadBrickD", experimentalHall_log_, false, 0);
        G4Colour colour_lead;
        G4Colour::GetColour("lead", colour_lead);
        G4VisAttributes* leadBrickVisAtt = new G4VisAttributes(colour_lead);
        leadBrick_log->SetVisAttributes(leadBrickVisAtt);
    }
    
}

void g4PSIDetectorConstruction::Place_Beamline(G4double bl_zPos = -140*cm) {
    
    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* Aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    //G4Material* Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    G4Material* Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
    
    // --- Beamline --- //
    // Estimated dimensions from tape measurements at PSI 26/July/2012
    // Vacuum: 12'' outer diameter, 2'' thick aluminum walls
    // Beamline extends from start of downstream world to Z = ~1.4 m (wrt target)
    // Beam exit window is 200 microns of Kapton
    G4double bl_rMin = 15.24*cm;
    G4double bl_rMax = 20.32*cm;
    G4double bl_rMinVac = bl_rMin; // - 0.1*mm;
    //    G4double bl_L = 360*cm; // why so long?
    G4double bl_L = 100*cm;
    //  G4double bl_zPos = -90*cm;  // end of the beam line
    //  G4double bl_dExitWindow = 0.2*mm; // vacuum exit window, added downstream on the tube
    
    // Konrad Deiters: e-mail 02/09/15: The beam exit window is 190 micro Mylar
    
    G4double bl_dExitWindow = 0.190*mm; // vacuum exit window, added downstream on the tube
    
    // don't let beamline extend out-of-world
    bl_L = std::min(bl_L, bl_zPos + expHall_z_ );
    
    G4Tubs* bl_tub = new G4Tubs("beam_tub", bl_rMin, bl_rMax, bl_L/2., 0, twopi);
    G4LogicalVolume* bl_log = new G4LogicalVolume(bl_tub, Aluminum, "beamline_log", 0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0.,0., bl_zPos - bl_L/2.), bl_log, "beamline_phys", experimentalHall_log_,false,0);
    
    G4Tubs* bl_tub_vac = new G4Tubs("beam_tub_vac", 0.0, bl_rMinVac, bl_L/2., 0, twopi);
    G4LogicalVolume* bl_vac_log = new G4LogicalVolume(bl_tub_vac, Vacuum, "beamline_vacuum_log", 0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0.,0., bl_zPos - bl_L/2.), bl_vac_log, "beamline_vacuum_phys", experimentalHall_log_,false,0);
    
    G4Tubs* bl_exit_win = new G4Tubs("beam_exit_win", 0.0, bl_rMax, bl_dExitWindow/2., 0, twopi);
    //    G4LogicalVolume* bl_exit_log = new G4LogicalVolume(bl_exit_win, Kapton, "beamline_exit_log", 0,0,0);
    G4LogicalVolume* bl_exit_log = new G4LogicalVolume(bl_exit_win, Mylar, "beamline_exit_log", 0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0.,0., bl_zPos + bl_dExitWindow/2.), bl_exit_log, "beamline_exit_phys", experimentalHall_log_, false,0);
    
    G4Colour color_beamline;
    G4Colour::GetColour("aluminum", color_beamline);
    bl_log->SetVisAttributes(new G4VisAttributes(color_beamline));
    bl_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    // G4Colour color_beamline_vac;
    // G4Colour::GetColour("target", color_beamline_vac);
    // bl_vac_log->SetVisAttributes(new G4VisAttributes(color_beamline_vac));
    
    G4Colour color_beamline_exit;
    G4Colour::GetColour("plastic", color_beamline_exit);
    bl_exit_log->SetVisAttributes(new G4VisAttributes(color_beamline_exit));
    //bl_vac_log->SetVisAttributes(G4VisAttributes::Invisible);
    
    //bl_exit_log->SetVisAttributes(G4VisAttributes::Invisible);
}

void g4PSIDetectorConstruction::DefineExpHall(G4double dx = 4*m,
                                              G4double dy = 3*m,
                                              G4double dz = 5*m,
                                              G4String material = "G4_AIR")
{
    // --- Experimental Hall (world volume) --- //
    // --- beam line along z axis --- //
    
    G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (mat == 0) {
        G4Exception("g4PSIDetectorConstruction::DefineExpHall",
                    ("Error: unknown material " + material).c_str(),
                    FatalException, "");
    };
    
    expHall_x_ = dx;
    expHall_y_ = dy;
    expHall_z_ = dz;
    G4Box* experimentalHall_box = new G4Box("expHall_solid",expHall_x_,expHall_y_,expHall_z_);
    experimentalHall_log_ = new G4LogicalVolume(experimentalHall_box, mat, "expHall_log");
    experimentalHall_phys_ = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log_,"expHall",0,false,0);
    G4VisAttributes *visAtt = new G4VisAttributes(G4Colour(0.,1.,1.));
    visAtt->SetForceWireframe(true);
    experimentalHall_log_->SetVisAttributes(visAtt);
    experimentalHall_log_->SetVisAttributes (G4VisAttributes::Invisible);
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    DetectorParts->SetDefaultMotherVolume(experimentalHall_log_);
}

void g4PSIDetectorConstruction::ConstructMUSE2015()
{
    DefineExpHall();
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    // -------------------------------- //
    //  Define beamline and shielding   //
    // -------------------------------- //
    
    if (DetectorParts->build[g4PSIDetectorParts::kBeamline]) {
        Place_Beamline();
    }
    if (DetectorParts->build[g4PSIDetectorParts::kStructure]) {
        Place_TableAndFrames();
    }
    
    //    Place_Shielding();
    
    // ------------------------------ //
    //  Define the target/detectors   //
    // ------------------------------ //
    
    
    //------------------------------ Scattering chamber and Target
    
    // Target Parameters go : Target Radius, Kapton thickness, pipe distance, pipe angle
    
//    g4PSITarget *tgt = new g4PSITarget("TGT", DetectorParts->GetUserT2Par(), DetectorParts->GetUserT1Par(), DetectorParts->GetUserT3Par(), DetectorParts->GetUserT4Par(), DetectorParts->GetUserT5Par());
//    double zoffset = tgt->GetUpstreamPos();
//    double target_radius = tgt->GetTargetRadius();
//    DetectorParts->Add(tgt);
    
    double zoffset = -20.0038 * cm;
    double target_radius = 3.0 * cm;
    
    // Scattering chamber
    // ------------------
    
    g4PSIChamberBase* chamber = NULL;
    if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeCylinder]) {
        chamber = new g4PSIChamberCylinder("CHAMBER");
    } else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeWindowlessCylinder]) {
        G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        G4double outer_radius = 16 * cm;
        G4double wall_thickness = 4 * mm;
        chamber = new g4PSIChamberCylinder("CHAMBER", mat, outer_radius, wall_thickness);
    } else if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_TypeTrapezoid]) {
        chamber = new g4PSIChamberTrapezoid("CHAMBER");
    }
    if (chamber) DetectorParts->Add(chamber);
    
    
    // Target
    // ------
    
    g4PSITargetBase* target = NULL;
    if (DetectorParts->build[g4PSIDetectorParts::kTarget_Type0]) {
        target = new g4PSITargetCylinder("TGT", chamber);
    } else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeUMich]) {
        
        target = new g4PSITargetUM("TGT", chamber);
        
    } else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeru]) {
        
        G4double tgt_R = 30*mm;
        G4double tgt_H = 80*mm;
        
        G4double base_r1 = 61*mm / 2;
        G4double base_r2 = 43.8*mm / 2;
        G4double base_r3 = 60*mm / 2;
        
        G4double base_h1 = 3*mm;
        G4double base_h2 = 9*mm;
        G4double base_h3 = 5*mm;
        G4double tgt_wall_t = 75*um;  // Dany Horovitz e-mail, 03/13/16
        
        G4Material* base_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        G4Material* tgt_wall_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        
        target = new g4PSITargetJeru("TGT", tgt_R, tgt_H, base_r1, base_r2, base_r3, base_h1, base_h2, base_h3, base_material, tgt_wall_t, tgt_wall_material, chamber);
    
    } else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeJeruCu]) {

        G4double tgt_R = 30*mm;
        G4double tgt_H = 80*mm;
        
        G4double base_r1 = 61*mm / 2;
        G4double base_r2 = 43.8*mm / 2;
        G4double base_r3 = 60*mm / 2;
        
        G4double base_h1 = 3*mm;
        G4double base_h2 = 9*mm;
        G4double base_h3 = 5*mm;
        G4double tgt_wall_t = 75*um;  // Dany Horovitz e-mail, 03/13/16
        
        G4Material* base_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
        G4Material* tgt_wall_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        
        target = new g4PSITargetJeru("TGT", tgt_R, tgt_H, base_r1, base_r2, base_r3, base_h1, base_h2, base_h3, base_material, tgt_wall_t, tgt_wall_material, chamber);
        
    } else if (DetectorParts->build[g4PSIDetectorParts::kTarget_TypeGrid]) {
        target = new g4PSITargetGrid("TGT", chamber);
    }
    if (target) DetectorParts->Add(target);
    
    //------------------------------ Beam Cherenkov Detectors
    if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_quartz]) {
        G4Material *quartz = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        DetectorParts->Add(new g4PSIBeamCherenkov("BCC", quartz,  -65*cm + zoffset, +45*deg, 16*cm, 10*cm, 3*mm));
    } else if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc]) {
        G4Material *scintillator = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        DetectorParts->Add(new g4PSIBeamCherenkov("BCC", scintillator,  -65*cm + zoffset, 90*deg, 16*cm, 10*cm, 2*mm));
    } else if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc2] ||
               DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc3] ||
               DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc4]) {
        G4Material *sc = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

        G4double gap = 6 * um;   // updated from 12 um, 01/09/17
        G4double w1 = 4 * mm;    // updated from 8 mm, 01/09/17
        G4double w2 = 8 * mm;    // updated from 8 mm, 01/09/17
        
        const G4double l = 10. * cm;
        const G4double t = 2 * mm;
        const G4double zsep = 2 * cm;
        G4double zdown = -0*cm  -33 * cm + zoffset;
        if (DetectorParts->GetUserI1Par() == 1000) {
            w1 = DetectorParts->GetUserD1Par();
            zdown = DetectorParts->GetUserD2Par() - 33*cm + zoffset;
        }
        G4int N1 = 6;  // inner paddles
        G4int N2 = 5;  // N2 outer paddles on each side

        if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc2]) {
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC1", sc,   -zsep + zdown, 90*deg, l,   N1, w1, N2, w2, t,  0*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC2", sc,           zdown, 90*deg, l,   N1, w1, N2, w2, t, 90*deg, gap));
        } else if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc3]) {
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC1", sc, -2*zsep + zdown, 90*deg, l,   N1, w1, N2, w2, t,  0*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC2", sc,   -zsep + zdown, 90*deg, l,   N1, w1, N2, w2, t, 90*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC3", sc,           zdown, 90*deg, l, N1+1, w1, N2, w2, t,  0*deg, gap));
        } else if (DetectorParts->build[g4PSIDetectorParts::kBeamCherenkov_sc4]) {
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC1", sc, -3*zsep + zdown, 90*deg, l,   N1, w1, N2, w2, t,  0*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC2", sc, -2*zsep + zdown, 90*deg, l,   N1, w1, N2, w2, t, 90*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC3", sc,   -zsep + zdown, 90*deg, l, N1+1, w1, N2, w2, t,  0*deg, gap));
            DetectorParts->Add(new g4PSIBeamCherenkov("BCC4", sc,           zdown, 90*deg, l, N1+1, w1, N2, w2, t, 90*deg, gap));
        }
    }
    
    //------------------------------ Target Scintillating Fiber Array
    
    double zSciFi = -45 * cm + zoffset;
    Place_SciFi_Detectors(zSciFi);
    
    //------------------------------ Beam GEM Detectors
    if (DetectorParts->build[g4PSIDetectorParts::kBeamGEM1])
        DetectorParts->Add(new g4PSIGEM("GEM1", 10.0*cm, 12.5*cm, -26.8*cm + zoffset));
    if (DetectorParts->build[g4PSIDetectorParts::kBeamGEM2])
        DetectorParts->Add(new g4PSIGEM("GEM2", 10.0*cm, 12.5*cm, -18.4*cm + zoffset));
    if (DetectorParts->build[g4PSIDetectorParts::kBeamGEM3])
        DetectorParts->Add(new g4PSIGEM("GEM3", 10.0*cm, 12.5*cm, -10.0*cm + zoffset));
    
    //------------------------------ Upstream beam scintillator ring, veto detector
    if (DetectorParts->build[g4PSIDetectorParts::kBeamVetoSC]) {
        G4double dz = 1*cm; // 4*mm;
        G4double z = -5*cm + zoffset; // previously -10*cm;
        G4int nSC = 8;
        G4int nPlane = 1;
        G4double rMin = target_radius > 0 ? target_radius : 2.000*cm; // rMin = 30 mm with Jeru Tgt.
        G4double rMax = 16*cm; // previously 7.243*cm;
        DetectorParts->Add(new g4PSIVeto("VSC", dz, z, nSC, rMin, rMax, nPlane));
        DetectorParts->Add(new g4PSITestPlane("VSC_TEST", 0, z-dz/2-1.0*mm,  100*cm, 100*cm));
    }
    
    //------------------------------ Downstream beam scintillators
    
    if (DetectorParts->build[g4PSIDetectorParts::kBeamMonitorSC]) {
        G4double z = 120*cm;
        DetectorParts->Add(new g4PSIBeamMonitor("BLSC", z));
        DetectorParts->Add(new g4PSITestPlane("BLSC_TEST", 0, z-6.0*cm,  150*cm, 150*cm));
    }
    
    //------------------------------ Downstream Calorimeter
    
    if (DetectorParts->build[g4PSIDetectorParts::kCal]) {
        G4double z = 130*cm;
        G4double tPb = 0. * mm;
        G4double tCal = 20. * cm;
        G4int nLayers = 1;
        g4PSICal::CalType ct = g4PSICal::kPbGlass;
        
        if (DetectorParts->GetUserI1Par() == 1110) {
            ct = g4PSICal::kPbSC;
            nLayers = DetectorParts->GetUserI2Par();
            tPb = DetectorParts->GetUserD1Par();
            tCal = DetectorParts->GetUserD2Par();
        }
        DetectorParts->Add(new g4PSICal("CAL", ct, z, tPb, tCal, nLayers));
    }
    
    //------------------------------ Downstream NaI detector
    
    if (DetectorParts->build[g4PSIDetectorParts::kNaI]) {
        G4double z = 130*cm;
        G4double t = 0. * mm;
        if (DetectorParts->GetUserI1Par() == 1100) {
            t = DetectorParts->GetUserD1Par();
        }
        DetectorParts->Add(new g4PSINaI("NAI", z, t));
    }
    
    //------------------------------ The SC Scintillator Walls
    
    G4double thetaC = 60*deg;
    G4double sc_al_thickness = 0.0;
    
    if (DetectorParts->build[g4PSIDetectorParts::kSCWall1]) {
        
        G4double sc_dx = 6.0*cm;
//        G4double sc_r = 52.0*cm;  // original
        G4double sc_r = 57.0*cm;  // e-mail Tom O'Connor 08/29/2017
        G4double sc_small_gap = 0.68*mm;
        G4double sc_wide_gap = 0.92*mm;
        
        G4double sc_dz = 3.0*cm;
        G4double sc_height = 120*cm;
//        G4double sc_n_up = 8;
//        G4double sc_n_down = 9;
//        G4double shift = -3.0 * cm;
        G4double sc_units = 9;
        
        G4double shift = DetectorParts->GetUserI1Par() == 100 ? DetectorParts->GetUserD1Par() : 0.0 * cm;
        
        DetectorParts->Add(new g4PSISPS("SPSLF", g4PSISPS::front, sc_dx,sc_small_gap, sc_wide_gap, sc_dz, sc_height, sc_r, thetaC, shift, sc_units));
        DetectorParts->Add(new g4PSISPS("SPSRF", g4PSISPS::front, sc_dx, sc_small_gap, sc_wide_gap, sc_dz, sc_height, sc_r, -thetaC, shift, sc_units));
    }
    if (DetectorParts->build[g4PSIDetectorParts::kSCWall2]) {
        
        G4double sc_dx = 6.0*cm;
        G4double sc_dz = 6.0*cm;
//        G4double sc_r = 74.0*cm;   // original
//        G4double sc_r = 77.5*cm;     // e-mail Tom O'Connor 08/29/2017
        G4double sc_r = 79*cm;     // e-mail Tom O'Connor 01/05/2018
        
        G4double sc_small_gap = 0.53*mm;
        G4double sc_wide_gap = 0.77*mm;
        
        G4double sc_height = 220*cm;
//        G4double sc_n_up = 15;
//        G4double sc_n_down = 12;
//        G4double shift = 9.0 * cm;
        G4double sc_units = 14;

        G4double shift = DetectorParts->GetUserI2Par() == 200 ? DetectorParts->GetUserD2Par() : 9.0 * cm;

        DetectorParts->Add(new g4PSISPS("SPSLR", g4PSISPS::rear, sc_dx, sc_small_gap, sc_wide_gap, sc_dz, sc_height, sc_r, thetaC, shift, sc_units));
        DetectorParts->Add(new g4PSISPS("SPSRR", g4PSISPS::rear, sc_dx, sc_small_gap, sc_wide_gap, sc_dz, sc_height, sc_r, -thetaC, shift, sc_units));
    }
    
    //------------------------------ The Straw Tube Chambers
    
    if (DetectorParts->build[g4PSIDetectorParts::kWireChamber1]) {
        DetectorParts->Add(new g4PSISTT("STTL", thetaC));
        DetectorParts->Add(new g4PSISTT("STTR", -thetaC));
    }
    
//    if (DetectorParts->build[g4PSIDetectorParts::kWireChamber1]) {
//        DetectorParts->Add(new g4PSISTC("STCL1", thetaC, g4PSISTC::kFront));
//        DetectorParts->Add(new g4PSISTC("STCR1", -thetaC, g4PSISTC::kFront));
//    }
//    if (DetectorParts->build[g4PSIDetectorParts::kWireChamber2]) {
//        DetectorParts->Add(new g4PSISTC("STCR2", -thetaC, g4PSISTC::kRear));
//        DetectorParts->Add(new g4PSISTC("STCL2", thetaC, g4PSISTC::kRear));
//    }
    
    //------------------------------ The Testplanes
    
    if (DetectorParts->build[g4PSIDetectorParts::kTestPlanes]) {
        double tp_sz = 1000*cm;
        DetectorParts->Add(new g4PSITestPlane("TSTBL_A", 0, -100*cm + zoffset, tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_B", 0, -75*cm + zoffset,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_C", 0, -35*cm + zoffset,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_D", 0, -15*cm + zoffset,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_E", 0, -0.5*cm + zoffset,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_F", 0, +25*cm,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_G", 0, +55*cm,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_H", 0, +95*cm,  tp_sz, tp_sz));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_I", 0, +150*cm, tp_sz, tp_sz));
        //
        //        //	DetectorParts->Add(new g4PSITestPlane("TSTBL_m1490", 0, -149*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_m1200", 0, -120*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_m810", 0, -81*cm, tp_sz, tp_sz));
        //        //	DetectorParts->Add(new g4PSITestPlane("TSTBL_m790", 0, -79*cm, tp_sz, tp_sz));
        //        //	DetectorParts->Add(new g4PSITestPlane("TSTBL_m710", 0, -71*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_m620", 0, -62*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_m420", 0, -42*cm, tp_sz, tp_sz));
        //		DetectorParts->Add(new g4PSITestPlane("TSTBL_m220", 0, -22*cm, tp_sz, tp_sz));
        //        // DetectorParts->Add(new g4PSITestPlane("TSTBL_m21", 0, -2.1*cm, 12*cm, 12*cm, mother_for_target_log));
        //        // DetectorParts->Add(new g4PSITestPlane("TSTBL_p21", 0, +2.1*cm, 12*cm, 12*cm, mother_for_target_log));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_p200", 0, +20*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_p500", 0, +50*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_p970", 0, +97*cm, tp_sz, tp_sz));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_p1500", 0, +150*cm, tp_sz, tp_sz));
    } else if (DetectorParts->build[g4PSIDetectorParts::kTestPlanes_full_setup]) {
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_B", 0, -36*cm + zoffset, 300*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_C", 0, -31*cm + zoffset, 90*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_D", 0, -15*cm + zoffset, 60*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_E", 0, +20*cm + zoffset,  30*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_F", 0, +55*cm + zoffset,  30*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_G", 0, +95*cm + zoffset,  55*cm, 300*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_H", 0, +150*cm + zoffset, 500*cm, 500*cm));
    } else if (DetectorParts->build[g4PSIDetectorParts::kTestPlanes_beamline_detectors]) {
        DetectorParts->Add(new g4PSITestPlane("TSTBL_BCC", 0, DetectorParts->GetDetector("BCC1")->GetZPos() - 1*cm, 90*cm, 300*cm));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM1", 0, DetectorParts->GetDetector("GEM1")->GetZPos() - 1.5*cm, 80*cm, 300*cm));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM2", 0, DetectorParts->GetDetector("GEM2")->GetZPos() - 1.5*cm, 70*cm, 300*cm));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM3", 0, DetectorParts->GetDetector("GEM3")->GetZPos() - 1.5*cm, 60*cm, 300*cm));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_TGT", 0, zoffset - 2*cm, 50*cm, 300*cm));
        DetectorParts->Add(new g4PSITestPlane("TSTBL_BLSC", 0, DetectorParts->GetDetector("BLSC")->GetZPos() - 4*cm, 50*cm, 300*cm));
        
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_BCC", 0, -75*cm + zoffset, 200*cm, 200*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_BCC_SciFi", 0, -55*cm + zoffset, 200*cm, 200*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_SciFi_GEM1", 0, -35*cm + zoffset, 90*cm, 200*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM1_GEM2", 0, -25*cm +zoffset, 75*cm, 200*cm));
        //		DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM2_GEM3", 0, -15*cm +zoffset, 60*cm, 200*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_GEM3_Veto", 0, -7*cm +zoffset, 50*cm, 200*cm));
        //        DetectorParts->Add(new g4PSITestPlane("TSTBL_Veto_Target", 0, -2*cm +zoffset,  40*cm, 200*cm));
        
    }
    
    //------------------------------------------------------------------
    DetectorParts->Info();
    DetectorParts->Placement();  // placement into the default volume, if not specified otherwise
    DetectorParts->SetSD(SDman);
    
    // \todo
    // DetectorParts->AddKillPrimaryLog(DetectorParts->GetDetector("TSTBL_BLSC")->GetDetectorLVolume());
    //------------------------------------------------------------------
}

void g4PSIDetectorConstruction::ConstructBeamtest2013(double angle, bool tgtIN)
{
    DefineExpHall();
    
    G4double GEMangle = angle > 40*deg ? 40*deg : angle;
    G4double SCangle = angle;
    G4double TgtAngle = angle > 40*deg ? 45*deg : 0*deg;
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    //------------------------------ Beamline
    Place_Beamline(-136.8*cm);
    Place_Shielding();
    Place_Target(TgtAngle, tgtIN);
    
    //  Place_BeamtestSupport();
    //------------------------------ Cherenkov
    //DetectorParts->Add(new g4PSIBeamCherenkov("BCC", -12.5*cm));
    G4Material *sapphir = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    G4Material *quartz = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    DetectorParts->Add(new g4PSIBeamCherenkov("BCC1", sapphir, -12.5*cm, -45*deg, 10, 8*cm, 3*mm));
    DetectorParts->Add(new g4PSIBeamCherenkov("BCC2", quartz,   5.5*cm, +45*deg, 10, 8*cm, 3*mm));
    //------------------------------ SC Bar
    if (angle == 0) {Place_SC("SC", 0*cm - 2.5*cm, SCangle);}
    else {Place_SC("SC", 150*cm - 2.5*cm, SCangle);}
    //------------------------------ GEM telescopes
    Place_GEMtelescope("TeleUp", 0*cm, -126.3*cm, 0.*deg, 1.*cm, 33.0*cm, 76.3*cm, 81.8*cm);
    Place_GEMtelescope("TeleDown", 30.5*cm, 0.*cm, GEMangle, 1.*cm, 34.5*cm, 76.5*cm, 82.*cm);
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructTOFTest2014()
{
    DefineExpHall();
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    G4Material *quartz = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* Aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    
    //------------------------------ Beamline
    Place_Beamline();
    Place_Shielding();
    
    double zBCC = -30 * cm;
    //    double zFrontSC = zBCC + 335 * cm;
    //    double zFrontSC = zBCC + 235 * cm;
    double zFrontSC = zBCC + 236.5 * cm;
    
    DetectorParts->Add(new g4PSISCWall("SCsmall", 20*mm, 20*mm, 2*mm, 0.0, zBCC-30*cm));
    DetectorParts->Add(new g4PSIBeamCherenkov("BCC", quartz, zBCC, +45*deg, 24, 3*cm, 3*mm));
    Place_SC("SC", zFrontSC + 2.5*cm, 0.0);
    
    const bool do_structure = true;
    if (do_structure) {
        G4Box* structure_box
        = new G4Box("structure", 0.5*cm, 5*cm, 0.5*cm);  // half lengths
        //	= new G4Box("structure", 10*cm, 2*cm, 10*cm);
        G4LogicalVolume* structure_log = new G4LogicalVolume(structure_box,
                                                             Aluminum, "structure_log");
        new G4PVPlacement(0, G4ThreeVector(0*cm, -6.0*cm, zBCC - 30*cm),
                          //	new G4PVPlacement(0, G4ThreeVector(0., 0*cm, zFrontSC + 60.*cm),
                          structure_log,"floor", experimentalHall_log_, false, 0);
        G4Colour colour_structure;
        G4Colour::GetColour("aluminum", colour_structure);
        G4VisAttributes* structureVisAtt = new G4VisAttributes(colour_structure);
        structure_log->SetVisAttributes(structureVisAtt);
    }
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructTOFTest2015()
{
    DefineExpHall();
    double beam_line_window = -254.8*cm;
    
    Place_Beamline(beam_line_window);
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    G4bool checkOverlaps = true;
    G4NistManager* nist = G4NistManager::Instance();
    
    //
    // scintillator
    //
    G4Material* plastic = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    //    G4Box* scBox = new G4Box("Scint", 2.5*cm, 12.5*cm, 2.5*cm);
    //    G4LogicalVolume* scLog = new G4LogicalVolume(scBox, plastic, "Scint");
    //    G4ThreeVector scPos  = G4ThreeVector(0*cm, 0*cm, (-57.2-2.5)*cm);
    //    new G4PVPlacement(0, scPos, scLog, "Scint", experimentalHall_log_, false, 0, checkOverlaps);
    DetectorParts->Add(new g4PSISCWall("Scint", 50*mm, 250*mm, 50*mm, 0.0, 2.5*cm));
    
    // Window
    G4Box* windowSolid = new G4Box("WindowSolid", 15.*cm, 15.*cm, 0.15*mm);
    G4LogicalVolume* windowLog = new G4LogicalVolume(windowSolid, plastic, "WindowLog");
    
    G4ThreeVector wn1Pos = G4ThreeVector(0*cm, 0*cm, -(193.1+4.5)*cm);
    G4ThreeVector wn2Pos = G4ThreeVector(0*cm, 0*cm, -(193.1+4.5+30)*cm);
    new G4PVPlacement(0, wn1Pos, windowLog, "windowPhys1", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(0, wn2Pos, windowLog, "windowPhys2", experimentalHall_log_, false, 0, checkOverlaps);
    
    // SIPM
    //    G4Box* sipmBox = new G4Box("Scint", 5*cm, 1.5*cm, 1*mm);
    //    G4LogicalVolume* sipmPhys = new G4LogicalVolume(sipmBox, plastic, "Scint");
    
    //    G4ThreeVector sipm1  = G4ThreeVector(0*cm, 0*cm, 7.5*mm+move);
    //    G4ThreeVector sipm2  = G4ThreeVector(0*cm, 0*cm, 22.5*mm+move);
    //    G4ThreeVector sipm3  = G4ThreeVector(0*cm, 0*cm, 37.5*mm+move);
    //    new G4PVPlacement(rot, sipm1, sipmPhys, "SIPM1", experimentalHall_log_, false, 0, checkOverlaps);
    //    new G4PVPlacement(rot, sipm2, sipmPhys, "SIPM2", experimentalHall_log_, false, 0, checkOverlaps);
    //    new G4PVPlacement(rot, sipm3, sipmPhys, "SIPM3", experimentalHall_log_, false, 0, checkOverlaps);
    
    
    // SiPM Detectors
    // ==============
    
    // rotation matrix
    // tan(theta) = 7.5mm/50cm, theta ~ 0.85937224364468 degree
    G4double theta = 0.8593722*deg, phi = 0.*deg;
    G4double x1 = cos(theta)*cos(phi), x2 = -sin(phi), x3 = sin(theta)*cos(phi);
    G4double y1 = cos(theta)*sin(phi), y2 =  cos(phi), y3 = sin(theta)*sin(phi);
    G4double z1 = -sin(theta),         z2 =  0,        z3 = cos(theta);
    G4ThreeVector xAxis(x1, x2, x3);
    G4ThreeVector yAxis(y1, y2, y3);
    G4ThreeVector zAxis(z1, z2, z3);
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);
    rot->invert();
    
    G4double move = -111.70*cm + 50.0*cm - DetectorParts->GetUserD1Par();
    
    DetectorParts->Add(new g4PSISCWall("SIPMU", 100*mm, 5*mm, 2*mm, 0.0, ( 7.5*mm+move), rot));
    DetectorParts->Add(new g4PSISCWall("SIPMM", 100*mm, 5*mm, 2*mm, 0.0, (22.5*mm+move), rot));
    DetectorParts->Add(new g4PSISCWall("SIPMD", 100*mm, 5*mm, 2*mm, 0.0, (37.5*mm+move), rot));
    
    // Tedlar material
    G4Element* C = new G4Element("Carbon"  , "C", 6., 12.011*g/mole);
    G4Element* H = new G4Element("Hydrogen", "H", 1., 1.0079*g/mole);
    G4Element* F = new G4Element("Fluoride", "F", 9., 18.998*g/mole);
    G4Material* Tedlar = new G4Material("PolyvinylFluoride", 1.5*g/cm3, 3);
    Tedlar->AddElement(C, 2);
    Tedlar->AddElement(H, 3);
    Tedlar->AddElement(F, 1);
    
    // Tedlar
    G4Box* tedBox = new G4Box("TedBox", 5*cm, 1.5*cm, 0.05*mm);
    G4LogicalVolume* tedlog = new G4LogicalVolume(tedBox, Tedlar, "tedLog");
    G4Colour colour;
    G4Colour::GetColour("tedlar", colour);
    G4VisAttributes* tedVisAtt = new G4VisAttributes(colour);
    tedlog->SetVisAttributes(tedVisAtt);
    
    G4ThreeVector ted1Pos = G4ThreeVector(0*cm, 0*cm, ( 2.+move)*mm);
    G4ThreeVector ted2Pos = G4ThreeVector(0*cm, 0*cm, (13.+move)*mm);
    G4ThreeVector ted3Pos = G4ThreeVector(0*cm, 0*cm, (17.+move)*mm);
    G4ThreeVector ted4Pos = G4ThreeVector(0*cm, 0*cm, (28.+move)*mm);
    G4ThreeVector ted5Pos = G4ThreeVector(0*cm, 0*cm, (32.+move)*mm);
    G4ThreeVector ted6Pos = G4ThreeVector(0*cm, 0*cm, (43.+move)*mm);
    new G4PVPlacement(rot, ted1Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(rot, ted2Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(rot, ted3Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(rot, ted4Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(rot, ted5Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    new G4PVPlacement(rot, ted6Pos, tedlog, "tedLog", experimentalHall_log_, false, 0, checkOverlaps);
    
    double dz = DetectorParts->GetUserD2Par();
    if (dz > 0) {
        //        DetectorParts->Add(new g4PSISCWall("Degrader", 10*cm, 10*cm, dz, 0.0, beam_line_window + 1*cm + 0.5*dz));
        DetectorParts->Add(new g4PSISCWall("Degrader", 10*cm, 10*cm, dz, 0.0, move - 1*cm + 0.5*dz));
    }
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructBeamProfile()
{
    DefineExpHall();
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    //------------------------------ Beamline
    Place_Beamline();
    //    Place_Shielding();
    
    //    DetectorParts->Add(new g4PSISCWall("SCsmall", 40*cm, 40*cm, 1*cm, 0.0, -99.*cm));
    
    //------------------------------ Beam GEM Detectors
    // DetectorParts->Add(new g4PSIGEM("GEM1", 10.0*cm, 12.5*cm, -40.0*cm));
    // DetectorParts->Add(new g4PSIGEM("GEM2", 10.0*cm, 12.5*cm,   0.0*cm));
    // DetectorParts->Add(new g4PSIGEM("GEM3", 10.0*cm, 12.5*cm,  40.0*cm));
    // Place_SC("SC", 66.*cm, 0.0);
    
    //------------------------------ The Testplanes
    
    double dx = 100*cm;
    DetectorParts->Add(new g4PSITestPlane("TSTBL_m240", 0, -240*mm, dx, dx));
    DetectorParts->Add(new g4PSITestPlane("TSTBL_m120", 0, -120*mm, dx, dx));
    DetectorParts->Add(new g4PSITestPlane("TSTBL_m0", 0, 0*cm, dx, dx));
    DetectorParts->Add(new g4PSITestPlane("TSTBL_p120", 0,  120*mm, dx, dx));
    DetectorParts->Add(new g4PSITestPlane("TSTBL_p240", 0,  240*mm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_m041", 0,  -41*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_m020", 0,  -20*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_m001", 0,    1*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_p020", 0,   20*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_p041", 0,   41*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_p060", 0,   60*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_p080", 0,   80*cm, dx, dx));
    // DetectorParts->Add(new g4PSITestPlane("TSTBL_p100", 0,  100*cm, dx, dx));
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructScatteringTest()
{
    DefineExpHall(5*m, 5*m, 10*m, "G4_Galactic");
    //    DefineExpHall();
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    DetectorParts->AllOff();
    DetectorParts->build[g4PSIDetectorParts::kTarget_Type0] = true;
    DetectorParts->Add(new g4PSITarget("TGT"));
    DetectorParts->Add(new g4PSITestSphere("DET", 100*CLHEP::cm, 0*CLHEP::MeV));
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructNaITest()
{
    //DefineExpHall(5*m, 5*m, 10*m, "G4_Galactic");
    //    DefineExpHall();
    DefineExpHall(5*m, 5*m, 10*m, "G4_AIR");
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    // -------------------------------- //
    //  Define beamline and shielding   //
    // -------------------------------- //
    
    // Place_Beamline();
    //Place_Shielding();
    
    // ------------------------------ //
    //  Define the target/detectors   //
    // ------------------------------ //
    
    //    G4LogicalVolume* mother_for_target_log = experimentalHall_log_;
    
    //------------------------------ Scattering chamber
    // if (DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type1] ||
    // 	DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type2] ||
    // 	DetectorParts->build[g4PSIDetectorParts::kScatteringChamber_Type3]
    // 	) Place_ScatteringChamber(mother_for_target_log);
    
    //    DetectorParts->Add(new g4PSIGEM("GEM", 10.0*cm, 12.5*cm, 10.0*cm));
    
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* NaI = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    
    const bool lin = false;
    if (lin) {
        G4bool checkOverlaps = true;
        
        //
        // Envelope
        //
        
        // rotation matrix
        G4double theta = 5.*deg;
        G4double ax1 = 1., ax2 = 0., ax3 = 0.;
        G4double ay1 = 0., ay2 = cos(theta), ay3 = -sin(theta);
        G4double az1 = 0., az2 = sin(theta), az3 = cos(theta);
        
        G4ThreeVector xAxis(ax1, ax2, ax3);
        G4ThreeVector yAxis(ay1, ay2, ay3);
        G4ThreeVector zAxis(az1, az2, az3);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateAxes(xAxis, yAxis, zAxis);
        rot->invert();
        
        G4double x1 = (10.16/2.)*cm, y1 = 5.08*cm, z1 = (40.64/2.)*cm;
        G4ThreeVector pos1 = G4ThreeVector(0*cm, 0*cm, -25*cm);
        G4Box* Solid1 = new G4Box("solid1", x1, y1, z1);
        G4LogicalVolume* Log1 = new G4LogicalVolume(Solid1, Al, "log1");
        new G4PVPlacement(rot, pos1, Log1, "Phys1", experimentalHall_log_, false, 0, checkOverlaps);
        
        G4double x2 = (10.16/2.)*cm, y2 = (5.08/2.)*cm, z2= (40.64/2.)*cm;
        G4Box* Solid2 = new G4Box("solid2", x2, y2, z2);
        G4LogicalVolume* Log2 = new G4LogicalVolume(Solid2, Al, "log2");
        G4LogicalVolume* Log3 = new G4LogicalVolume(Solid2, Al, "log2");
        G4ThreeVector pos2 = G4ThreeVector(0*cm, (5.08/2.)*cm, 0*cm);
        G4ThreeVector pos3 = G4ThreeVector(0*cm, -(5.08/2.)*cm, 0*cm);
        new G4PVPlacement(0, pos2, Log2, "Phys2", Log1, false, 0, checkOverlaps);
        new G4PVPlacement(0, pos3, Log3, "Phys3", Log1, false, 0, checkOverlaps);
        
        G4ThreeVector pos4 = G4ThreeVector(0*cm, 0*cm, 0*cm);
        G4double x3 = (10.16/2.-0.1)*cm, y3 = (5.08/2.-0.1)*cm, z3= (40.64/2.-0.1)*cm;
        G4Box* Solid3 = new G4Box("solid3", x3, y3, z3);
        G4LogicalVolume* Log4 = new G4LogicalVolume(Solid3, NaI, "log4");
        new G4PVPlacement(0, pos4, Log4, "Phys4", Log2, false, 0, checkOverlaps);
        new G4PVPlacement(0, pos4, Log4, "Phys4", Log3, false, 0, checkOverlaps);
        
        //
        // 2x2 Plastic plate
        //
        G4Material* shape2_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        G4Box* solidb = new G4Box("Block", 1*cm, 1*cm, 0.5*mm);
        G4ThreeVector posPP = G4ThreeVector(0*cm, 0*cm, 0*cm);
        G4LogicalVolume* sBlock = new G4LogicalVolume(solidb, shape2_mat, "Block");
        new G4PVPlacement(0, posPP, sBlock, "Block", experimentalHall_log_, false, 0, checkOverlaps);
        
    } else {
        
        // rotation matrix
        G4double theta = 5.*deg;
        G4double ax1 = 1., ax2 = 0., ax3 = 0.;
        G4double ay1 = 0., ay2 = cos(theta), ay3 = -sin(theta);
        G4double az1 = 0., az2 = sin(theta), az3 = cos(theta);
        
        G4ThreeVector xAxis(ax1, ax2, ax3);
        G4ThreeVector yAxis(ay1, ay2, ay3);
        G4ThreeVector zAxis(az1, az2, az3);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateAxes(xAxis, yAxis, zAxis);
        rot->invert();
        
        const double inch = 2.54 * cm;
        const double xNaI =  4 * inch;
        const double yNaI =  2 * inch;
        const double zNaI = 16 * inch;
        const double dx = 1.016 * mm;
        const double dy = 1.016 * mm;
        const double dz = 3.175 * mm;
        G4Box* NaIHousing1 = new G4Box("NaI1_geom", xNaI/2.+dx, yNaI/2.+dy, zNaI/2.+dz);
        G4LogicalVolume* NaILog1 = new G4LogicalVolume(NaIHousing1, Al, "NaIHousing1_log");
        G4Box* NaIHousing2 = new G4Box("NaI2_geom", xNaI/2.+dx, yNaI/2.+dy, zNaI/2.+dz);
        G4LogicalVolume* NaILog2 = new G4LogicalVolume(NaIHousing2, Al, "NaIHousing2_log");
        
        const double x = 0 * cm;
        const double y = (yNaI + dy) / 2.0;
        const double z = 25 * cm;
        const double y0 = 1.0 * yNaI/2.;
        
        new G4PVPlacement(rot, G4ThreeVector(x, y0+y, z), NaILog1, "NaIHousing1_phys", experimentalHall_log_, false, 0);
        g4PSISCWall* NaIDetector1 = new g4PSISCWall("NaI1", xNaI, yNaI, zNaI, 0.0, 0.0, NULL, NaI);
        NaIDetector1->SetMotherVolume(NaILog1);
        DetectorParts->Add(NaIDetector1);
        new G4PVPlacement(rot, G4ThreeVector(x, y0-y, z), NaILog2, "NaIHousing2_phys", experimentalHall_log_, false, 0);
        g4PSISCWall* NaIDetector2 = new g4PSISCWall("NaI2", xNaI, yNaI, zNaI, 0.0, 0.0, NULL, NaI);
        NaIDetector2->SetMotherVolume(NaILog2);
        DetectorParts->Add(NaIDetector2);
        
        DetectorParts->Add(new g4PSISCWall("sci", 2*cm, 2*cm, 1*mm, 0*cm, 0*cm));
        
        DetectorParts->Info();
        DetectorParts->Placement();
        DetectorParts->SetSD(SDman);
    }
    
    // G4int idet = 3;
    
    // DefineExpHall(5*m, 5*m, 10*m, "G4_Galactic");
    
    // g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    // G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    // G4Material *sapphir = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE"); // sapphir
    // G4Material *quartz = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); // quartz
    
    // switch (idet) {
    // case 1:
    // 	DetectorParts->Add(new g4PSIBeamCherenkov("BCC", sapphir, 0*cm, +45*deg, 10, 8*cm, 3*mm));
    // 	break;
    // case 2:
    // 	DetectorParts->Add(new g4PSIBeamCherenkov("BCC", quartz, 0*cm, +45*deg, 10, 8*cm, 3*mm));
    // 	break;
    // case 3:
    // 	Place_GEMtelescope("TeleUp", 0*cm, -8*cm, 0.*deg, 5.*cm, 8*cm, 11*cm, -13*cm);
    // 	break;
    // default:;
    // }
    // DetectorParts->Add(new g4PSITestPlane("TSTBL", 0, +500*cm, 100*cm, 100*cm));
    
    // DetectorParts->Info();
    // DetectorParts->Placement();
    // DetectorParts->SetSD(SDman);
}

void g4PSIDetectorConstruction::ConstructGeantTest()
{
    //DefineExpHall(5*m, 5*m, 10*m, "G4_Galactic");
    //    DefineExpHall();
    DefineExpHall(5*m, 5*m, 10*m, "G4_AIR");
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    // -------------------------------- //
    //  Define beamline and shielding   //
    // -------------------------------- //
    
    Place_Beamline();
    //Place_Shielding();
    
    // G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    // G4Material* NaI = G4NistManager::Instance()->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    
    
    DetectorParts->Add(new g4PSISCWall("sci", 5*cm, 50*cm, 5*cm, 0*cm, 0*cm));
    
    G4double thetaC = 0*deg;
    G4double sc_al_thickness = 0.0;
    
    G4double sc_dx = 6.0*cm;
    G4double sc_r = 52.0*cm;
    G4double sc_gap = 0.0*cm;
    
    G4double sc_dz = 3.0*cm;
    G4double sc_height = 120*cm;
    G4double sc_n_up = 1;
    G4double sc_n_down = 1;
    
    DetectorParts->Add(new g4PSISCWall("SC1", sc_dx, sc_gap, sc_dz, sc_height,
                                       sc_r, thetaC, sc_n_up, sc_n_down, sc_al_thickness));

    sc_dx = 6.0*cm;
    sc_dz = 6.0*cm;
    sc_r = 74.0*cm;
    sc_gap = 0.0*cm;
    sc_height = 220*cm;
    sc_n_up = 1;
    sc_n_down = 1;
    
    DetectorParts->Add(new g4PSISCWall("SC2", sc_dx, sc_gap, sc_dz, sc_height,
                                       sc_r, thetaC, sc_n_up, sc_n_down, sc_al_thickness));
    
    DetectorParts->Info();
    DetectorParts->Placement();
    DetectorParts->SetSD(SDman);
    
}

/*
void g4PSIDetectorConstruction::WriteGDML()
{

    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();

    //if (DetectorParts->GetUserB1Par() == true) {
        
        //Instatiate the GDML parser
        G4GDMLParser parser;        

        //Segfault is thrown if GDML file alredy exists, lets check if the file exists
        struct stat buffer;   
        if (!(stat ("MUSE.gdml", &buffer) == 0)){
            G4cout << "Preparing to write GDML file!" << G4endl;
            parser.Write("MUSE.gdml", experimentalHall_log_);
        }

        else G4cout << "MUSE.gdml file already exists. Not writing geometry." << G4endl;

    //}

    //#endif

}
*/

G4VPhysicalVolume* g4PSIDetectorConstruction::Construct()
{
    //-------------------//
    //      Volumes      //
    //-------------------//
    
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    
    switch (DetectorParts->GetSetup()) {
        case g4PSIDetectorParts::kStandard2015:
            ConstructMUSE2015();
            break;
        case g4PSIDetectorParts::kTest2013_0:
            ConstructBeamtest2013(0.*deg, false);
            break;
        case g4PSIDetectorParts::kTest2013_20:
            ConstructBeamtest2013(20.*deg, false);
            break;
        case g4PSIDetectorParts::kTest2013_20t:
            ConstructBeamtest2013(20.*deg, true);
            break;
        case g4PSIDetectorParts::kTest2013_40:
            ConstructBeamtest2013(40.*deg, false);
            break;
        case g4PSIDetectorParts::kTest2013_40t:
            ConstructBeamtest2013(40.*deg, true);
            break;
        case g4PSIDetectorParts::kTOFTest2014:
            ConstructTOFTest2014();
            break;
        case g4PSIDetectorParts::kTOFTest2015:
            ConstructTOFTest2015();
            break;
        case g4PSIDetectorParts::kBeamProfile:
            ConstructBeamProfile();
            break;
        case g4PSIDetectorParts::kScatteringTest:
            ConstructScatteringTest();
            break;
        case g4PSIDetectorParts::kNaITest:
            ConstructNaITest();
            break;
        case g4PSIDetectorParts::kGeantTest:
            ConstructGeantTest();
            break;
        default:;
    }
    
    //WriteGDML();
    //do_the_gdml();

    //DetectorParts->Print();
    return experimentalHall_phys_;
}
