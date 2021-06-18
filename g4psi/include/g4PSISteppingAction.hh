#ifndef g4PSISteppingAction_h
#define g4PSISteppingAction_h 1

#define MAKEHBOOK 1

#include "G4UserSteppingAction.hh"

class g4PSIDetectorConstruction;
class g4PSIEventAction;
class g4PSIRunAction;

class g4PSISteppingAction : public G4UserSteppingAction
{
public:
  g4PSISteppingAction();
  // g4PSISteppingAction(g4PSIDetectorConstruction*,
  // 			  g4PSIEventAction*,
  // 			  g4PSIRunAction*);
  ~g4PSISteppingAction();
  
  void UserSteppingAction(const G4Step*);
  
private:
  
  // g4PSIDetectorConstruction* detector;
  // g4PSIEventAction* eventaction;
  // g4PSIRunAction* runaction;
};

#endif
