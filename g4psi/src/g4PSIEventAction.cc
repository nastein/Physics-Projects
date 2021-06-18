#include "g4PSIEventAction.hh"
#include "g4PSIAnalysisManager.hh"

g4PSIEventAction::g4PSIEventAction() {
}


g4PSIEventAction::~g4PSIEventAction() {
    G4cout << "ending g4PSIEventAction\n";
}


void g4PSIEventAction::BeginOfEventAction(const G4Event* evt) {
    g4PSIAnalysisManager::getInstance()->BeginOfEvent(evt);
}


void g4PSIEventAction::EndOfEventAction(const G4Event* evt) {
    g4PSIAnalysisManager::getInstance()->EndOfEvent(evt);
}

