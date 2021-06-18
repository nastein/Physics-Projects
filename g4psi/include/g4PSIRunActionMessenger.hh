#ifndef g4PSIRunActionMessenger_h
#define g4PSIRunActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class g4PSIRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class g4PSIRunActionMessenger: public G4UImessenger
{
public:
  g4PSIRunActionMessenger(g4PSIRunAction*);
  virtual ~g4PSIRunActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  g4PSIRunAction*         runAction;
  G4UIdirectory*        runActionDir; 
  G4UIcmdWithAString*   SetFileCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

