#ifndef g4PSIDetectorMessenger_h
#define g4PSIDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class g4PSIDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class g4PSIDetectorMessenger: public G4UImessenger
{
  public:
    g4PSIDetectorMessenger();
   ~g4PSIDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    G4UIdirectory*             g4PSIDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        DetCompOnCmd;
    G4UIcmdWithAString*        DetCompOffCmd;
    G4UIcmdWithAString*        DetInfo;
    G4UIcmdWithAString*        DetTrigCmd;
    G4UIcmdWithAString*        TargetStateCmd;
    G4UIcmdWithAString*        SetupCmd;
    G4UIcmdWithADoubleAndUnit* DetUserD1Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserD2Cmd;
    G4UIcmdWithAnInteger*      DetUserI1Cmd;
    G4UIcmdWithAnInteger*      DetUserI2Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserT1Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserT2Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserT3Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserT4Cmd;
    G4UIcmdWithADoubleAndUnit* DetUserT5Cmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

