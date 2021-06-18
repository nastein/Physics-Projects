#ifndef g4PSIDetectorParts_h
#define g4PSIDetectorParts_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIDetectorParts
//!
//! This soliton class stores the global access to setup parts
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIBOptrMultiParticleChangeCrossSection.hh"
#include <iostream>
#include <fstream>

class g4PSIDetectorParts {
    
public:
    
    /// Selection of predefined setups
    
    enum Setup {
        kStandard2015,
        kTest2013_0,
        kTest2013_20,
        kTest2013_20t,
        kTest2013_40,
        kTest2013_40t,
        kTOFTest2014,
        kTOFTest2015,
        kBeamProfile,
        kScatteringTest,
        kNaITest,
        kGeantTest
    };
    
    Setup setup_;
    
    /// BuildModes lists detector geometries which can be individually turned on or
    /// off at the init phase of the simulation
    /// \todo include more parts of the setup
    
    enum BuildModes {
        // keep kShieldFloor always first in the list and kTestPlanes always last, or update Write(), Print(), AllOff()
        kShieldFloor,
        kShieldWall,
        kShieldBLConcrete,
        kShieldBLLead,
        kShieldBLConcreteTube,
        kStructure,
        kBeamline,
        
        kScatteringChamber_Type1,
        kScatteringChamber_Type2,
        kScatteringChamber_Type3,
        kScatteringChamber_Type4,  //dummy Indiana model
        kScatteringChamber_TypeUM,  // UMich model
        kScatteringChamber_TypeJeru, //Jerusalem configuration
        kScatteringChamber_TypeJeru2, //Jerusalem configuration cylindrical
        
        kScatteringChamber_TypeCylinder,
        kScatteringChamber_TypeWindowlessCylinder,
        kScatteringChamber_TypeTrapezoid,
        
        kWireChamber1,
        kWireChamber2,
        kWireChamber_noframe,
        kWireChamber_rfshield,
        kOutgoingGEM1,
        kOutgoingGEM2,
        kOutgoingGEM3,
        kSCWall1,
        kSCWall2,
        kPMTs,
        kNaI,
        kCal,
        kTarget_Type0,
        kTarget_Type1,
        kTarget_Type2,
        kTarget_Type3,  // upright coffee-can
        kTarget_Type4,  // tuna-can
        kTarget_TypeUM,  // UMich target
        kTarget_TypeUM2,
        
        kTarget_TypeCylinder,
        kTarget_TypeJeru, //Jerusalem configuration
        kTarget_TypeJeruCu, //Jerusalem configuration
        kTarget_TypeUMich,  // UMich target
        kTarget_TypeUMich10cm,  // UMich target
        kTarget_TypeUMich11cm,  // UMich target
        kTarget_TypeGrid,

        kEmptyTarget_Type1,
        kEmptyTarget_Type2,
        kEmptyTarget_Type3,
        kEmptyTarget_Type4,
        kEmptyTarget_TypeUM, // UMich target
        kEmptyTarget_TypeUM2,
        kEmptyTarget_TypeJeru, //Jerusalem configuration
        
        kBeamGEM1,
        kBeamGEM2,
        kBeamGEM3,
        kBeamSciFi,
        kBeamSciFi_Type0,
        kBeamSciFi_Type1,
        kBeamSciFi_Type2,
        kBeamSciFi_Type3,
        kBeamVetoSC,
        kBeamMonitorSC,
        kBeamCherenkov_vmode,
        kBeamCherenkov_quartz,
        kBeamCherenkov_sapphire,
        kBeamCherenkov_sc,
        kBeamCherenkov_sc2,   // two planes
        kBeamCherenkov_sc3,   // three planes
        kBeamCherenkov_sc4,   // four planes
        kTestPlanes_beamline_detectors,
        kTestPlanes_full_setup,
        kTestPlanes
    };
    
    std::map<BuildModes,bool> build;
    
    enum TargetState {
        kFull,
        kEmpty,
        kFoil
    };
    
    
    
    
public:
    static g4PSIDetectorParts* getInstance();
    static void dispose();
    
    g4PSIDetectorBase* GetDetector(G4String name);
    
    void Print();
    void AllOff();
    void Add(g4PSIDetectorBase* d);
    void Add(G4LogicalVolume *m, g4PSIDetectorBase* d);
    
    void SetTriggerDetectors(G4String name);
    void AddRequiredDetector(G4String name);
    void BeginOfRun();
    
    void Info();
    void Placement(G4LogicalVolume *m);
    void Placement();
    void SetSD(G4SDManager *sdman);
    void SetDefaultMotherVolume(G4LogicalVolume *mother) {default_mother_ = mother;};
    
    void RootInitTree(TTree *t);   // called at the begin of run
    void RootWrite();   // called at the end of run
    void DeleteEventData();   // called at the end of run
    
    TargetState GetTargetState();
    void SetTargetState(TargetState state);
    void SetSetup(Setup setup) {setup_ = setup;};
    Setup GetSetup() {return setup_;};
    G4bool HasTrigger();  /// if true, fill root tree with event
    
    void AttachBiasingOperator(G4LogicalVolume *lvolume);
    
    void SetWikiMarkupFile(G4String s);
    
    void AddKillPrimaryLog(G4LogicalVolume *lvolume);
    
    void SetUserD1Par(G4double p);
    void SetUserD2Par(G4double p);
    G4double GetUserD1Par() {return userDouble1Par_;};
    G4double GetUserD2Par() {return userDouble2Par_;};
    void SetUserI1Par(G4int i);
    G4int GetUserI1Par() {return userInt1Par_;};
    void SetUserI2Par(G4int i);
    G4int GetUserI2Par() {return userInt2Par_;};

    void SetUserT1Par(G4double p);
    void SetUserT2Par(G4double p);
    void SetUserT3Par(G4double p);
    void SetUserT4Par(G4double p);
    void SetUserT5Par(G4double p);

    G4double GetUserT1Par() {return userDoubleT1Par_;};
    G4double GetUserT2Par() {return userDoubleT2Par_;};
    G4double GetUserT3Par() {return userDoubleT3Par_;};
    G4double GetUserT4Par() {return userDoubleT4Par_;};
    G4double GetUserT5Par() {return userDoubleT5Par_;};
    
    
private:
    
    TargetState target_state_;
    

    typedef std::vector <g4PSIDetectorBase*> detector_list;
    
    detector_list detector_;
    std::vector <G4String> trigger_detectors_string_;
    std::vector <detector_list> trigger_detectors_;
    
    G4LogicalVolume *default_mother_;
    
    static g4PSIDetectorParts* fManager;
    g4PSIDetectorParts();    // private constructor
    ~g4PSIDetectorParts();
    g4PSIDetectorParts(const g4PSIDetectorParts&);   // Prevent copy-construction
    g4PSIDetectorParts& operator=(const g4PSIDetectorParts&);  // Prevent assignment}
    
    GB01BOptrMultiParticleChangeCrossSection* testMany_;

    G4LogicalVolume *target_;        // target material (hydrogen)
    G4LogicalVolume *target_wall_;   // target flask
    G4LogicalVolume *target_si_;     // super insulation
    
    G4String wiki_markup_filename_;
    G4bool do_wiki_markup_;
    std::ofstream wiki_markup_;
    
    std::vector <G4LogicalVolume*> kill_primary_log_;
    
    G4double userDouble1Par_;
    G4double userDouble2Par_;
    G4int userInt1Par_;
    G4int userInt2Par_;

    G4double userDoubleT1Par_;
    G4double userDoubleT2Par_;
    G4double userDoubleT3Par_;
    G4double userDoubleT4Par_;
    G4double userDoubleT5Par_;
};

#endif
