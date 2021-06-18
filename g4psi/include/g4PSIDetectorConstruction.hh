#ifndef TrekDetectorConstruction_h
#define TrekDetectorConstruction_h 1

#include <map>
#include "G4VUserDetectorConstruction.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class g4PSIDetectorMessenger;

class g4PSIDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    g4PSIDetectorConstruction();
    ~g4PSIDetectorConstruction();
  
    G4VPhysicalVolume* Construct();

private:
  
    G4double expHall_x_;
    G4double expHall_y_;
    G4double expHall_z_;
    G4LogicalVolume* experimentalHall_log_;
    G4VPhysicalVolume* experimentalHall_phys_;
    g4PSIDetectorMessenger* detectorMessenger_;
  
    G4UserLimits* StepLimit_;            // pointer to user step limits

    G4Colour G4HexColour(int r, int g, int b);
  
    void DefineExpHall(G4double dx, G4double dy, G4double dz, G4String material);
    //void WriteGDML();

    void Place_Beamline(G4double end_of_beamline_z);
    void Place_LH2_Target(G4LogicalVolume* mother_volume);
    void Place_Target(G4double &target_r, G4double &target_z, G4LogicalVolume* mother_volume);
    void Place_Target(G4double angle, G4bool in); // beam-test target
    void Place_ScatteringChamber(G4LogicalVolume* &lvolume);
    void Place_TableAndFrames();
    void Place_Shielding();
    void Place_BeamtestSupport();
    void Place_SciFi_Detectors(G4double z);
    void Place_GEMtelescope(G4String name, G4double r, G4double z0, G4double angle,
			    G4double r1, G4double r2, G4double r3, G4double rsc);
    void Place_SC(G4String name, G4double r, G4double angle);

    void ConstructMUSE();
    void ConstructMUSE2014();
    void ConstructMUSE2015();
    void ConstructBeamtest2013(double angle, bool tgt_in);
    void ConstructTOFTest2014();
    void ConstructTOFTest2015();
    void ConstructBeamProfile();
    void ConstructScatteringTest();
    void ConstructNaITest();
    void ConstructGeantTest();
};

#endif
