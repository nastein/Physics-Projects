#ifndef g4PSIPrimaryGeneratorMessenger_h
#define g4PSIPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class g4PSIPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class g4PSIPrimaryGeneratorMessenger: public G4UImessenger
{
public:
    g4PSIPrimaryGeneratorMessenger(g4PSIPrimaryGeneratorAction*);
    virtual ~g4PSIPrimaryGeneratorMessenger();
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    g4PSIPrimaryGeneratorAction* pParticleGun;
    G4UIdirectory* AnalysisDir;
    
    G4UIcmdWithADoubleAndUnit* BeamMomentumCmd;
    G4UIcmdWithADouble* BeamMomentumSpreadCmd;
    G4UIcmdWithAnInteger* NParticleCmd;
    G4UIcmdWithAString* GunModeCmd;
    G4UIcmdWithAString* BeamMixModeCmd;
    G4UIcmdWithAString* ExternalSeedCmd;
    G4UIcmdWithAString* ExternalSeedFileCmd;
    G4UIcmdWithAString* TurtleFileCmd;
    G4UIcmdWithADoubleAndUnit* PathLengthForRF;
    G4UIcmdWithADoubleAndUnit* SetZOriginCmd;
    G4UIcmdWithADoubleAndUnit* BeamRFJitterCmd;
    G4UIcmdWithADoubleAndUnit* TargetRadiusCmd;
    
    G4UIcmdWithADoubleAndUnit* SetResSCWallCmd;
    G4UIcmdWithADoubleAndUnit* SetResBeamCerenCmd;
    G4UIcmdWithADoubleAndUnit* SetResSFPlaneCmd;
    
};

///
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

