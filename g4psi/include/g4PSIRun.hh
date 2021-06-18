#ifndef g4PSIRun_h
#define g4PSIRun_h 1

#include "G4Run.hh"
#include "globals.hh"

class G4Event;

/// Run class
///

class g4PSIRun : public G4Run
{
public:
    g4PSIRun();
    virtual ~g4PSIRun();    
};


#endif
