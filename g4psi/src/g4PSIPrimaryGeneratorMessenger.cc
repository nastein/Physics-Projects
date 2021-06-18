#include "g4PSIPrimaryGeneratorMessenger.hh"

#include "g4PSIPrimaryGeneratorAction.hh"
#include "g4PSIAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIPrimaryGeneratorMessenger::g4PSIPrimaryGeneratorMessenger(g4PSIPrimaryGeneratorAction* g4PSIGun)
:pParticleGun(g4PSIGun)
{
    BeamMomentumCmd = new G4UIcmdWithADoubleAndUnit("/gun/beam_momentum",this);
    BeamMomentumCmd->SetGuidance("Set momentum of beam particle");
    BeamMomentumCmd->SetParameterName("p_in", false);
    BeamMomentumCmd->SetUnitCategory("Energy");
    BeamMomentumCmd->AvailableForStates(G4State_Idle);
    
    BeamMomentumSpreadCmd = new G4UIcmdWithADouble("/gun/beam_momentum_spread",this);
    BeamMomentumSpreadCmd->SetGuidance("Set relative momentum spread");
    BeamMomentumSpreadCmd->SetParameterName("p_in", false);
    BeamMomentumSpreadCmd->AvailableForStates(G4State_Idle);
    
    ExternalSeedCmd = new G4UIcmdWithAString("/gun/seeds",this);
    ExternalSeedCmd->SetGuidance("Set external seeds.  Enter two (long) integers.");
    ExternalSeedCmd->SetParameterName("seeds", false);
    ExternalSeedCmd->AvailableForStates(G4State_Idle);
    
    ExternalSeedFileCmd = new G4UIcmdWithAString("/gun/seedfile",this);
    ExternalSeedFileCmd->SetGuidance("Set external seed file.  Enter filename.");
    ExternalSeedFileCmd->SetParameterName("filename", false);
    ExternalSeedFileCmd->AvailableForStates(G4State_Idle);
    
    TurtleFileCmd = new G4UIcmdWithAString("/gun/turtle",this);
    TurtleFileCmd->SetGuidance("Set external turtle file.  Enter filename.");
    TurtleFileCmd->SetParameterName("filename", false);
    TurtleFileCmd->AvailableForStates(G4State_Idle);
    
    GunModeCmd = new G4UIcmdWithAString("/gun/mode",this);
    GunModeCmd->SetGuidance("Gun mode.");
    GunModeCmd->SetGuidance("  Choice : default_gun, walluniform, walluniform_momdist, walluniform_from_center, wallcenter, wallcenter_momdist, beamline, beamline_pencil, beamline_box, full4pi, cosmic");
    GunModeCmd->SetParameterName("choice", false);
    GunModeCmd->SetDefaultValue("default_gun");
    GunModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    NParticleCmd = new G4UIcmdWithAnInteger("/gun/nparticles",this);
    NParticleCmd->SetGuidance("Set number of primary particles");
    NParticleCmd->SetParameterName("max_nPrimaryParticle", false);
    NParticleCmd->AvailableForStates(G4State_Idle);
    
    BeamMixModeCmd = new G4UIcmdWithAString("/gun/beam_mix",this);
    BeamMixModeCmd->SetGuidance("Beam-mix mode.");
    BeamMixModeCmd->SetGuidance("  Choice : default_particle, beam_mix_plus, beam_mix_minus");
    BeamMixModeCmd->SetParameterName("choice", false);
    BeamMixModeCmd->SetDefaultValue("default_particle");
    BeamMixModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    PathLengthForRF = new G4UIcmdWithADoubleAndUnit("/gun/set_pathlength",this);
    PathLengthForRF->SetGuidance("Set flight path of channel to target center");
    PathLengthForRF->SetParameterName("path_length", false);
    PathLengthForRF->SetUnitCategory("Length");
    PathLengthForRF->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    SetZOriginCmd = new G4UIcmdWithADoubleAndUnit("/gun/set_vertexz",this);
    SetZOriginCmd->SetGuidance("Set Z origin for beamline modes wrt to target (0)");
    SetZOriginCmd->SetParameterName("vertz", false);
    SetZOriginCmd->SetUnitCategory("Length");
    SetZOriginCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    BeamRFJitterCmd = new G4UIcmdWithADoubleAndUnit("/gun/set_protonjitter",this);
    BeamRFJitterCmd->SetGuidance("Set proton pulse width for intrinsic RF jitter");
    BeamRFJitterCmd->SetParameterName("proton_jitter", true);
    BeamRFJitterCmd->SetDefaultValue(0.0*ns);
    BeamRFJitterCmd->SetUnitCategory("Time");
    BeamRFJitterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    TargetRadiusCmd = new G4UIcmdWithADoubleAndUnit("/gun/set_radius",this);
    TargetRadiusCmd->SetGuidance("Set radius of full4pi beam origin");
    TargetRadiusCmd->SetParameterName("target_radius",true);
    TargetRadiusCmd->SetDefaultValue(3.0*cm);
    TargetRadiusCmd->SetUnitCategory("Length");
    TargetRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    // Create "analysis" directory for the resolution controls
    AnalysisDir = new G4UIdirectory("/g4PSI/analysis/");
    AnalysisDir->SetGuidance("Analysis Manager Control");
    
    SetResSCWallCmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/analysis/setres_scwall",this);
    SetResSCWallCmd->SetGuidance("Set the tube resolution for the SC wall");
    SetResSCWallCmd->SetParameterName("scwall_res", false);
    SetResSCWallCmd->SetUnitCategory("Time");
    SetResSCWallCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    SetResBeamCerenCmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/analysis/setres_beamceren",this);
    SetResBeamCerenCmd->SetGuidance("Set the tube resolution for the beam Cerenkov");
    SetResBeamCerenCmd->SetParameterName("beamceren_res", false);
    SetResBeamCerenCmd->SetUnitCategory("Time");
    SetResBeamCerenCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    SetResSFPlaneCmd = new G4UIcmdWithADoubleAndUnit("/g4PSI/analysis/setres_sfplane",this);
    SetResSFPlaneCmd->SetGuidance("Set the plane resolution for the Sci-Fi");
    SetResSFPlaneCmd->SetParameterName("sfplane_res", false);
    SetResSFPlaneCmd->SetUnitCategory("Time");
    SetResSFPlaneCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIPrimaryGeneratorMessenger::~g4PSIPrimaryGeneratorMessenger() {
    
    delete BeamMixModeCmd;
    delete BeamMomentumCmd;
    delete BeamMomentumSpreadCmd;
    delete NParticleCmd;
    delete ExternalSeedCmd;
    delete ExternalSeedFileCmd;
    delete TurtleFileCmd;
    delete GunModeCmd;
    delete PathLengthForRF;
    delete SetZOriginCmd;
    delete BeamRFJitterCmd;
    delete TargetRadiusCmd;
    delete SetResSCWallCmd;
    delete SetResBeamCerenCmd;
    delete SetResSFPlaneCmd;
    delete AnalysisDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4PSIPrimaryGeneratorMessenger::SetNewValue(
                                                 G4UIcommand* command, G4String newValues) {
    
    if (command == GunModeCmd) {
        if (newValues.compare("walluniform") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::WALLUNIFORM);
        } else if (newValues.compare("walluniform_momdist") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::WALLUNIFORM_MOMDIST);
        } else if (newValues.compare("walluniform_from_center") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::WALLUNIFORM_FROM_CENTER);
        } else if (newValues.compare("wallcenter") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::WALLCENTER);
        } else if (newValues.compare("wallcenter_momdist") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::WALLCENTER_MOMDIST);
        } else if (newValues.compare("beamline_pencil") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_PENCIL);
        } else if (newValues.compare("beamline_pencil_flat") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_PENCIL_FLAT);
        } else if (newValues.compare("beamline_box") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_BOX);
        } else if (newValues.compare("beamline") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_2015_CENTER);
        } else if (newValues.compare("beamline_2014") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_2014);
        } else if (newValues.compare("beamline_2015") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_2015);
        } else if (newValues.compare("beamline_2015_center") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::BEAMLINE_2015_CENTER);
        } else if (newValues.compare("full4pi") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::FULL4PI);
        } else if (newValues.compare("full4pi_momdist") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::FULL4PI_MOMDIST);
        } else if (newValues.compare("full4pi_from_center") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::FULL4PI_FROM_CENTER);
        } else if (newValues.compare("full4pi_default") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::FULL4PI_DEFAULT);
        } else if (newValues.compare("target_muon_decay") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::TARGET_MUON_DECAY);
        } else if (newValues.compare("cosmic") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::COSMIC);
        } else if (newValues.compare("default") == 0) {
            pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::DEFAULT_GUN);
        } else {
            G4Exception("g4PSIPrimaryGeneratorMessenger::SetNewValue",
                        ("Unknown gun mode: " + newValues).c_str(),
                        FatalException, "");
        }
        
        
    } else if (command == TurtleFileCmd) {
        pParticleGun->SetTurtleFile(newValues);
        pParticleGun->SetGunMode(g4PSIPrimaryGeneratorAction::TURTLE);
        
    } else if (command == BeamMixModeCmd) {
        std::cout << ">>> " << newValues << "\n";
        if (newValues.compare("beam_mix_plus") == 0) {
            pParticleGun->SetBeamMixMode(g4PSIPrimaryGeneratorAction::BEAM_MIX_PLUS);
        } else if (newValues.compare("beam_mix_minus") == 0) {
            pParticleGun->SetBeamMixMode(g4PSIPrimaryGeneratorAction::BEAM_MIX_MINUS);
        } else if (newValues.compare("default_particle") == 0) {
            pParticleGun->SetBeamMixMode(g4PSIPrimaryGeneratorAction::DEFAULT_PARTICLE);
        } else {
            std::cout << "ERROR: unknown particle mode; set to default\n";
            pParticleGun->SetBeamMixMode(g4PSIPrimaryGeneratorAction::DEFAULT_PARTICLE);
        }
        
    } else if (command == BeamMomentumCmd) {
        G4double p_in = BeamMomentumCmd->GetNewDoubleValue(newValues);
        pParticleGun->SetBeamMomentumCenter(p_in);
        
    } else if (command == BeamMomentumSpreadCmd) {
        G4double dpp = BeamMomentumSpreadCmd->GetNewDoubleValue(newValues);
        pParticleGun->SetBeamMomentumSpread(dpp);
        
    } else if (command == NParticleCmd) {
        G4int n = NParticleCmd->GetNewIntValue(newValues);
        pParticleGun->SetNPrimaryParticles(n);
        
    } else if (command == ExternalSeedCmd) {
        char *pEnd = NULL;
        G4long s1 = strtol(newValues.c_str(), &pEnd, 10);
        G4long s2 = strtol(pEnd, &pEnd, 10);
        pParticleGun->SetEventSeeds(s1, s2, false, "");
        
    } else if (command == ExternalSeedFileCmd) {
        pParticleGun->SetEventSeeds(0L, 0L, true, newValues);
        
    } else if (command == PathLengthForRF) {
        G4double path_length = PathLengthForRF->GetNewDoubleValue(newValues);
        pParticleGun->SetPathLengthForRF(path_length);
        
    } else if (command == SetZOriginCmd) {
        G4double vertz = SetZOriginCmd->GetNewDoubleValue(newValues);
        pParticleGun->SetBeamlineVertexZ(vertz);
        
    } else if (command == BeamRFJitterCmd) {
        G4double width = BeamRFJitterCmd->GetNewDoubleValue(newValues);
        pParticleGun->SetRFJitter(width);
    } else if (command == TargetRadiusCmd) {
        G4double r = TargetRadiusCmd->GetNewDoubleValue(newValues);
        pParticleGun->SetTargetRadius(r);
    }
    
    if (command == SetResSCWallCmd) {
        g4PSIAnalysisManager::getInstance()->SetResSCWall( SetResSCWallCmd->GetNewDoubleValue(newValues) );
    }
    if (command == SetResBeamCerenCmd) {
        g4PSIAnalysisManager::getInstance()->SetResBeamCeren( SetResBeamCerenCmd->GetNewDoubleValue(newValues) );
    }
    if (command == SetResSFPlaneCmd) {
        g4PSIAnalysisManager::getInstance()->SetResSFPlane( SetResSFPlaneCmd->GetNewDoubleValue(newValues) );
    }
}
