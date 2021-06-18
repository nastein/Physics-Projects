#include "g4PSIRunAction.hh"
#include "g4PSIRunActionMessenger.hh"
#include "g4PSIAnalysisManager.hh"
#include "g4PSIDetectorParts.hh"
#include "g4PSIRun.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "G4VSteppingVerbose.hh"

g4PSIRunAction::g4PSIRunAction() {
    root_file_name_ = "default_run_action.root";
    runActionMessenger = new g4PSIRunActionMessenger(this);
}

g4PSIRunAction::~g4PSIRunAction() {
    G4cout << "ending g4PSIRunAction" << G4endl;
}

void g4PSIRunAction::BeginOfRunAction(const G4Run*)
{
    // long seeds[2];
    // seeds[0] = (long) 956686776;
    // seeds[1] = (long) 1208586836;
    // CLHEP::HepRandom::setTheSeeds(seeds);
    // CLHEP::HepRandom::showEngineStatus();
    // G4cout << "^^^ RunAction ^^^\n";
    
    G4cout << "g4PSIRunAction: start of run" << G4endl;
    g4PSIAnalysisManager::getInstance()->BeginOfRun(root_file_name_);
    g4PSIDetectorParts::getInstance()->BeginOfRun();
}

void g4PSIRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << "g4PSIRunAction: end of run" << G4endl;
    g4PSIAnalysisManager::getInstance()->EndOfRun(aRun);
}

