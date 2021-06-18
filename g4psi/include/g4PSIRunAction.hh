#ifndef g4PSIRunAction_h
#define g4PSIRunAction_h 1

#include "g4PSIRunActionMessenger.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class g4PSIRunAction : public G4UserRunAction
{
public:
    g4PSIRunAction();
    ~g4PSIRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    void SetRootFileName(G4String name) {root_file_name_ = name;};
    
private:
    G4String root_file_name_;
    g4PSIRunActionMessenger* runActionMessenger;
};

#endif

