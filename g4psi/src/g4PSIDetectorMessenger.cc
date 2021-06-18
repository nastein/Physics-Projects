#include "g4PSIDetectorMessenger.hh"

#include "g4PSIDetectorConstruction.hh"
#include "g4PSIDetectorParts.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4ExceptionSeverity.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIDetectorMessenger::g4PSIDetectorMessenger()
{
    g4PSIDir = new G4UIdirectory("/g4PSI/");
    g4PSIDir->SetGuidance("UI commands of MUSE Simulation");
    
    detDir = new G4UIdirectory("/g4PSI/det/");
    detDir->SetGuidance("Detector control");
    
    SetupCmd = new G4UIcmdWithAString("/g4PSI/det/setup",this);
    SetupCmd->SetGuidance("Choose default experimental setup [standard, standard2014, test2013_0, test2013_20, test2013_20t, test2013_40, test2013_40t].");
    SetupCmd->SetParameterName("choice",false);
    //  SetupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    SetupCmd->AvailableForStates(G4State_PreInit);
    
    DetCompOnCmd = new G4UIcmdWithAString("/g4PSI/det/component_on",this);
    DetCompOnCmd->SetGuidance("Turn on selected detector component.");
    DetCompOnCmd->SetParameterName("choice",false);
    //  DetCompOnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    DetCompOnCmd->AvailableForStates(G4State_PreInit);
    
    DetCompOffCmd = new G4UIcmdWithAString("/g4PSI/det/component_off",this);
    DetCompOffCmd->SetGuidance("Turn off selected detector component.");
    DetCompOffCmd->SetParameterName("choice",false);
    //  DetCompOffCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    DetCompOffCmd->AvailableForStates(G4State_PreInit);

    TargetStateCmd = new G4UIcmdWithAString("/g4PSI/det/target_state",this);
    TargetStateCmd->SetGuidance("Give target state (full / empty / foil)");
    TargetStateCmd->SetParameterName("choice",false);
    //  TargetStateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    TargetStateCmd->AvailableForStates(G4State_PreInit);
    
    DetTrigCmd = new G4UIcmdWithAString("/g4PSI/det/trigger",this);
    DetTrigCmd->SetGuidance("Name detectors in the trigger.");
    DetTrigCmd->SetParameterName("choice",false);
    //  DetTrigCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    DetTrigCmd->AvailableForStates(G4State_PreInit);
    
    DetInfo = new G4UIcmdWithAString("/g4PSI/det/info",this);
    DetInfo->SetGuidance("[wiki-markup-file].");
    DetInfo->SetParameterName("choice",false);
    //  DetInfo->AvailableForStates(G4State_PreInit,G4State_Idle);
    DetInfo->AvailableForStates(G4State_PreInit);
    
    DetUserD1Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/user_parameter1",this);
    DetUserD1Cmd->SetGuidance("User parameter with unit");
    DetUserD1Cmd->SetParameterName("user_parameterD1", true);
    DetUserD1Cmd->SetDefaultValue(0.0*cm);
    DetUserD1Cmd->SetUnitCategory("Length");
    DetUserD1Cmd->AvailableForStates(G4State_PreInit);

    DetUserD2Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/user_parameter2",this);
    DetUserD2Cmd->SetGuidance("User parameter with unit");
    DetUserD2Cmd->SetParameterName("user_parameterD2", true);
    DetUserD2Cmd->SetDefaultValue(0.0*cm);
    DetUserD2Cmd->SetUnitCategory("Length");
    DetUserD2Cmd->AvailableForStates(G4State_PreInit);
    
    DetUserI1Cmd = new G4UIcmdWithAnInteger("/g4PSI/det/user_integer1",this);
    DetUserI1Cmd->SetGuidance("User integer parameter");
    DetUserI1Cmd->SetParameterName("user_parameterI1", true);
    DetUserI1Cmd->SetDefaultValue(0);
    DetUserI1Cmd->AvailableForStates(G4State_PreInit);

    DetUserI2Cmd = new G4UIcmdWithAnInteger("/g4PSI/det/user_integer2",this);
    DetUserI2Cmd->SetGuidance("User integer parameter");
    DetUserI2Cmd->SetParameterName("user_parameterI2", true);
    DetUserI2Cmd->SetDefaultValue(0);
    DetUserI2Cmd->AvailableForStates(G4State_PreInit);

    DetUserT1Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/kapton_thickness",this);
    DetUserT1Cmd->SetGuidance("Target paramater with unit");
    DetUserT1Cmd->SetParameterName("kapton_thickness", true);
    DetUserT1Cmd->SetDefaultValue(120*um);
    DetUserT1Cmd->SetUnitCategory("Length");
    DetUserT1Cmd->AvailableForStates(G4State_PreInit);

    DetUserT2Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/target_radius",this);
    DetUserT2Cmd->SetGuidance("Target paramater with unit");
    DetUserT2Cmd->SetParameterName("target_radius", true);
    DetUserT2Cmd->SetDefaultValue(30*mm);
    DetUserT2Cmd->SetUnitCategory("Length");
    DetUserT2Cmd->AvailableForStates(G4State_PreInit);

    DetUserT3Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/pipe_distance",this);
    DetUserT3Cmd->SetGuidance("Target paramater with unit");
    DetUserT3Cmd->SetParameterName("pipe_distance", true);
    DetUserT3Cmd->SetDefaultValue(60*mm);
    DetUserT3Cmd->SetUnitCategory("Length");
    DetUserT3Cmd->AvailableForStates(G4State_PreInit);

    DetUserT4Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/pipe_angle",this);
    DetUserT4Cmd->SetGuidance("Target paramater with unit");
    DetUserT4Cmd->SetParameterName("pipe_angle", true);
    DetUserT4Cmd->SetDefaultValue(25.0*deg);
    DetUserT4Cmd->SetUnitCategory("Angle");
    DetUserT4Cmd->AvailableForStates(G4State_PreInit);

    DetUserT5Cmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/det/entrance_window_thickness",this);
    DetUserT5Cmd->SetGuidance("Target parameter with unit");
    DetUserT5Cmd->SetParameterName("entrance_window_thickness", true);
    DetUserT5Cmd->SetDefaultValue(200*um);
    DetUserT5Cmd->SetUnitCategory("Length");
    DetUserT5Cmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIDetectorMessenger::~g4PSIDetectorMessenger()
{
    delete SetupCmd;
    delete DetCompOnCmd;
    delete DetCompOffCmd;
    delete DetTrigCmd;
    delete detDir;
    delete g4PSIDir;
    delete DetInfo;
    delete DetUserD1Cmd;
    delete DetUserD2Cmd;
    delete DetUserI1Cmd;
    delete DetUserI2Cmd;
    delete DetUserT1Cmd;
    delete DetUserT2Cmd;
    delete DetUserT3Cmd;
    delete DetUserT4Cmd;
    delete DetUserT5Cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4PSIDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    g4PSIDetectorParts *DetectorParts = g4PSIDetectorParts::getInstance();
    if ( command == DetCompOnCmd || command == DetCompOffCmd ) {
        
        /// \todo make this safer agains inconsistent setups
        
        g4PSIDetectorParts::BuildModes mode;
        if (newValue == "ShieldFloor") mode = g4PSIDetectorParts::kShieldFloor;
        else if (newValue == "ShieldWall") mode = g4PSIDetectorParts::kShieldWall;
        else if (newValue == "ShieldBLConcrete") mode = g4PSIDetectorParts::kShieldBLConcrete;
        else if (newValue == "ShieldBLLead") mode = g4PSIDetectorParts::kShieldBLLead;
        else if (newValue == "ShieldBLConcreteTube") mode = g4PSIDetectorParts::kShieldBLConcreteTube;
        else if (newValue == "Beamline") mode = g4PSIDetectorParts::kBeamline;
        else if (newValue == "Structure") mode = g4PSIDetectorParts::kStructure;
        else if (newValue == "ScatteringChamber_Type1") mode = g4PSIDetectorParts::kScatteringChamber_Type1;
        else if (newValue == "ScatteringChamber_Type2") mode = g4PSIDetectorParts::kScatteringChamber_Type2;
        else if (newValue == "ScatteringChamber_Type3") mode = g4PSIDetectorParts::kScatteringChamber_Type3;
        else if (newValue == "ScatteringChamber_Type4") mode = g4PSIDetectorParts::kScatteringChamber_Type4;
        else if (newValue == "ScatteringChamber_TypeUM") mode = g4PSIDetectorParts::kScatteringChamber_TypeUM;
        else if (newValue == "ScatteringChamber_TypeJeru") mode = g4PSIDetectorParts::kScatteringChamber_TypeJeru; //Jersualem configuration
        else if (newValue == "ScatteringChamber_TypeJeru2") mode = g4PSIDetectorParts::kScatteringChamber_TypeJeru2; //Jersualem configuration (cylinder)
        else if (newValue == "ScatteringChamber_TypeTrapezoid") mode = g4PSIDetectorParts::kScatteringChamber_TypeTrapezoid;
        else if (newValue == "ScatteringChamber_TypeCylinder") mode = g4PSIDetectorParts::kScatteringChamber_TypeCylinder;
        else if (newValue == "ScatteringChamber_TypeWindowlessCylinder") mode = g4PSIDetectorParts::kScatteringChamber_TypeWindowlessCylinder;

        
        else if (newValue == "Target_Type0") mode = g4PSIDetectorParts::kTarget_Type0;
        else if (newValue == "Target_Type1") mode = g4PSIDetectorParts::kTarget_Type1;
        else if (newValue == "Target_Type2") mode = g4PSIDetectorParts::kTarget_Type2;
        else if (newValue == "Target_Type3") mode = g4PSIDetectorParts::kTarget_Type3;
        else if (newValue == "Target_Type4") mode = g4PSIDetectorParts::kTarget_Type4;
        else if (newValue == "Target_TypeUM") mode = g4PSIDetectorParts::kTarget_TypeUM;
        else if (newValue == "Target_TypeUM2") mode = g4PSIDetectorParts::kTarget_TypeUM2;//Secondary UM config
        
        else if (newValue == "Target_TypeUMich") mode = g4PSIDetectorParts::kTarget_TypeUMich; // UMich configuration
        else if (newValue == "Target_TypeUMich10cm") mode = g4PSIDetectorParts::kTarget_TypeUMich10cm; // UMich configuration
        else if (newValue == "Target_TypeUMich11cm") mode = g4PSIDetectorParts::kTarget_TypeUMich11cm; // UMich configuration
        else if (newValue == "Target_TypeJeru") mode = g4PSIDetectorParts::kTarget_TypeJeru; //Jerusalem configuration
        else if (newValue == "Target_TypeJeruCu") mode = g4PSIDetectorParts::kTarget_TypeJeruCu; // -"- Cu base
        else if (newValue == "Target_TypeCylinder") mode = g4PSIDetectorParts::kTarget_TypeCylinder;
        else if (newValue == "Target_TypeGrid") mode = g4PSIDetectorParts::kTarget_TypeGrid;
        
        else if (newValue == "EmptyTarget_Type1") mode = g4PSIDetectorParts::kEmptyTarget_Type1;
        else if (newValue == "EmptyTarget_Type2") mode = g4PSIDetectorParts::kEmptyTarget_Type2;
        else if (newValue == "EmptyTarget_Type3") mode = g4PSIDetectorParts::kEmptyTarget_Type3;
        else if (newValue == "EmptyTarget_Type4") mode = g4PSIDetectorParts::kEmptyTarget_Type4;
        else if (newValue == "EmptyTarget_TypeUM") mode = g4PSIDetectorParts::kEmptyTarget_TypeUM;
        else if (newValue == "EmptyTarget_TypeUM2") mode = g4PSIDetectorParts::kEmptyTarget_TypeUM2;//Secondary UM config
        else if (newValue == "EmptyTarget_TypeJeru") mode = g4PSIDetectorParts::kEmptyTarget_TypeJeru; //Jerusalem configuration

        else if (newValue == "WireChamber_noframe") mode = g4PSIDetectorParts::kWireChamber_noframe;
        else if (newValue == "WireChamber_rfshield") mode = g4PSIDetectorParts::kWireChamber_rfshield;
        else if (newValue == "WireChamber1") mode = g4PSIDetectorParts::kWireChamber1;
        else if (newValue == "WireChamber2") mode = g4PSIDetectorParts::kWireChamber2;
        else if (newValue == "OutgoingGEM1") mode = g4PSIDetectorParts::kOutgoingGEM1;
        else if (newValue == "OutgoingGEM2") mode = g4PSIDetectorParts::kOutgoingGEM2;
        else if (newValue == "OutgoingGEM3") mode = g4PSIDetectorParts::kOutgoingGEM3;
        else if (newValue == "SCWall1") mode = g4PSIDetectorParts::kSCWall1;
        else if (newValue == "SCWall2") mode = g4PSIDetectorParts::kSCWall2;
        else if (newValue == "PMTs") mode = g4PSIDetectorParts::kPMTs;
        else if (newValue == "NaI") mode = g4PSIDetectorParts::kNaI;
        else if (newValue == "Cal") mode = g4PSIDetectorParts::kCal;
        else if (newValue == "BeamGEM1") mode = g4PSIDetectorParts::kBeamGEM1;
        else if (newValue == "BeamGEM2") mode = g4PSIDetectorParts::kBeamGEM2;
        else if (newValue == "BeamGEM3") mode = g4PSIDetectorParts::kBeamGEM3;
        else if (newValue == "BeamSciFi") mode = g4PSIDetectorParts::kBeamSciFi;
        else if (newValue == "BeamSciFi_Type0") mode = g4PSIDetectorParts::kBeamSciFi_Type0;
        else if (newValue == "BeamSciFi_Type1") mode = g4PSIDetectorParts::kBeamSciFi_Type1;
        else if (newValue == "BeamSciFi_Type2") mode = g4PSIDetectorParts::kBeamSciFi_Type2;
        else if (newValue == "BeamSciFi_Type3") mode = g4PSIDetectorParts::kBeamSciFi_Type3;
        else if (newValue == "BeamVetoSC") mode = g4PSIDetectorParts::kBeamVetoSC;
        else if (newValue == "BeamMonitorSC") mode = g4PSIDetectorParts::kBeamMonitorSC;
        
        else if (newValue == "BeamCherenkov_vmode") mode = g4PSIDetectorParts::kBeamCherenkov_vmode;
        else if (newValue == "BeamCherenkov_quartz") mode = g4PSIDetectorParts::kBeamCherenkov_quartz;
        else if (newValue == "BeamCherenkov_sapphire") mode = g4PSIDetectorParts::kBeamCherenkov_sapphire;
        else if (newValue == "BeamCherenkov_sc") mode = g4PSIDetectorParts::kBeamCherenkov_sc;
        else if (newValue == "BeamCherenkov_sc2") mode = g4PSIDetectorParts::kBeamCherenkov_sc2;
        else if (newValue == "BeamCherenkov_sc3") mode = g4PSIDetectorParts::kBeamCherenkov_sc3;
        else if (newValue == "BeamCherenkov_sc4") mode = g4PSIDetectorParts::kBeamCherenkov_sc4;
        
        else if (newValue == "TestPlanes_full_setup") mode = g4PSIDetectorParts::kTestPlanes_full_setup;
        else if (newValue == "TestPlanes_beamline_detectors") mode = g4PSIDetectorParts::kTestPlanes_beamline_detectors;
        else if (newValue == "TestPlanes") mode = g4PSIDetectorParts::kTestPlanes;
        else
            G4Exception("g4PSIDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)",
                        ("Error: unknown detector component in component_on/off: " + newValue).c_str(),
                        FatalException, "");
        
        DetectorParts->build[mode] = (command == DetCompOnCmd);

        
    } else if ( command == TargetStateCmd ) {
        if (newValue == "full") DetectorParts->SetTargetState(g4PSIDetectorParts::kFull);
        if (newValue == "empty") DetectorParts->SetTargetState(g4PSIDetectorParts::kEmpty);
        if (newValue == "foil") DetectorParts->SetTargetState(g4PSIDetectorParts::kFoil);
        
        
    } else if ( command == DetTrigCmd ) {
        DetectorParts->SetTriggerDetectors(newValue);
        
        
    } else if ( command == SetupCmd ) {
        if (newValue == "standard2015") DetectorParts->SetSetup(g4PSIDetectorParts::kStandard2015);
        else if (newValue == "test2013_0") DetectorParts->SetSetup(g4PSIDetectorParts::kTest2013_0);
        else if (newValue == "test2013_20") DetectorParts->SetSetup(g4PSIDetectorParts::kTest2013_20);
        else if (newValue == "test2013_20t") DetectorParts->SetSetup(g4PSIDetectorParts::kTest2013_20t);
        else if (newValue == "test2013_40") DetectorParts->SetSetup(g4PSIDetectorParts::kTest2013_40);
        else if (newValue == "test2013_40t") DetectorParts->SetSetup(g4PSIDetectorParts::kTest2013_40t);
        else if (newValue == "toftest2014") DetectorParts->SetSetup(g4PSIDetectorParts::kTOFTest2014);
        else if (newValue == "toftest2015") DetectorParts->SetSetup(g4PSIDetectorParts::kTOFTest2015);
        else if (newValue == "beamprofile") DetectorParts->SetSetup(g4PSIDetectorParts::kBeamProfile);
        else if (newValue == "scatteringtest") DetectorParts->SetSetup(g4PSIDetectorParts::kScatteringTest);
        else if (newValue == "nai_test") DetectorParts->SetSetup(g4PSIDetectorParts::kNaITest);
        else if (newValue == "geant_test") DetectorParts->SetSetup(g4PSIDetectorParts::kGeantTest);
        else
            G4Exception("g4PSIDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)",
                        ("Error: unknown setup: " + newValue).c_str(),
                        FatalException, "");
    }   else if ( command == DetInfo ) {
        DetectorParts->SetWikiMarkupFile(newValue);
    }   else if (command == DetUserD1Cmd) {
        DetectorParts->SetUserD1Par(DetUserD1Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserD2Cmd) {
        DetectorParts->SetUserD2Par(DetUserD2Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserI1Cmd) {
        DetectorParts->SetUserI1Par(DetUserI1Cmd->GetNewIntValue(newValue));
    }   else if (command == DetUserI2Cmd) {
        DetectorParts->SetUserI2Par(DetUserI2Cmd->GetNewIntValue(newValue));
    }   else if (command == DetUserT1Cmd) {
        DetectorParts->SetUserT1Par(DetUserT1Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserT2Cmd) {
        DetectorParts->SetUserT2Par(DetUserT2Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserT3Cmd) {
        DetectorParts->SetUserT3Par(DetUserT3Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserT4Cmd) {
        DetectorParts->SetUserT4Par(DetUserT4Cmd->GetNewDoubleValue(newValue));
    }   else if (command == DetUserT5Cmd) {
        DetectorParts->SetUserT5Par(DetUserT5Cmd->GetNewDoubleValue(newValue));
    }
    
}

