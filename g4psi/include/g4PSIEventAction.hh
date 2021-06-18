#ifndef g4PSIEventAction_h
#define g4PSIEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class g4PSIRunAction;

class g4PSIEventAction : public G4UserEventAction
{
public:
  g4PSIEventAction();
  virtual ~g4PSIEventAction();
  
  /// This method is invoked before converting the primary particles to G4Track objects. 
  /// A typical use of this method would be to initialize and/or book histograms for a 
  /// particular event.
  void  BeginOfEventAction(const G4Event*);

  /// This method is invoked at the very end of event processing. 
  /// It is typically used for a simple analysis of the processed event.
  void    EndOfEventAction(const G4Event*);
};

#endif

    
