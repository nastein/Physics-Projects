#include "g4PSIRunActionMessenger.hh"

#include "g4PSIRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIRunActionMessenger::g4PSIRunActionMessenger(g4PSIRunAction* r)
:runAction(r)
{
    runActionDir = new G4UIdirectory("/g4PSI/run/");
    runActionDir->SetGuidance("RunAction control");
    
    SetFileCmd = new G4UIcmdWithAString("/g4PSI/run/rootfile",this);
    SetFileCmd->SetGuidance("Root filename.");
    SetFileCmd->SetGuidance("  Choice : g4PSI.root(default), off");
    SetFileCmd->SetParameterName("choice",true);
    SetFileCmd->SetDefaultValue("g4PSI.root");
    SetFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4PSIRunActionMessenger::~g4PSIRunActionMessenger()
{
    delete SetFileCmd;
    delete runActionDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4PSIRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if ( command == SetFileCmd ) {
        runAction->SetRootFileName(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

