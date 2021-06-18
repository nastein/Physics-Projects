#include "g4PSISteppingAction.hh"
#include "g4PSIRunAction.hh"
#include "g4PSIAnalysisManager.hh"
#include "g4PSIDetectorParts.hh"

#include "G4SteppingManager.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4UnitsTable.hh"

// g4PSISteppingAction::g4PSISteppingAction(g4PSIDetectorConstruction* det,
// 						 g4PSIEventAction* evt,
// 						 g4PSIRunAction* run)
//   :detector(det),
//    eventaction(evt),
//    runaction(run)
// {
// }

g4PSISteppingAction::g4PSISteppingAction()
{
}

g4PSISteppingAction::~g4PSISteppingAction(){
    G4cout << "ending g4PSISteppingAction\n";
}

void g4PSISteppingAction::UserSteppingAction(const G4Step* theStep)
{
    G4Track * theTrack = theStep->GetTrack();
    
    G4StepPoint *pre_pt  = theStep->GetPreStepPoint(); //  the current step
    G4StepPoint *post_pt  = theStep->GetPostStepPoint();
    
    const G4LogicalVolume *lvolume = pre_pt->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    
    const G4VProcess *post_proc = post_pt->GetProcessDefinedStep();  // get the process which has limited the current step.
    G4ParticleDefinition * particleType = theTrack->GetDefinition();
    
    if (particleType == G4MuonPlus::MuonPlusDefinition() ||
        particleType == G4MuonMinus::MuonMinusDefinition()) {
        //
        // Record the muon decay vertex and momentum direction of secondary particle (e-/e+)
        //
        if (fDecay == post_proc->GetProcessType()) {
            
            G4double x = theTrack->GetPosition().x();
            G4double y = theTrack->GetPosition().y();
            G4double z = theTrack->GetPosition().z();
            
            /// \todo understand this better; make it right. How is n = -1 possible?
            const G4TrackVector *secondary = theStep->GetSecondary();
            G4int n = (*secondary).size() - 1;
            G4ThreeVector pDir(0,0,0);
            G4double p = 0;
            
            if (n >= 0) {
                pDir = (*secondary)[n]->GetMomentumDirection();
                p = (*secondary)[n]->GetMomentum().mag();
            }
            
            // std::cout << "MOM: " << theTrack->GetDynamicParticle()->GetTotalMomentum() << "\n";
            // std::cout << (*secondary)[n]->GetDefinition()->GetParticleName() << "\n";
            // std::cout << x << " " << y << " " << z << " --- " << p << " size = " << (*secondary).size()<< "\n";
            
            G4double w = theTrack->GetWeight();
            
            g4PSIAnalysisManager* analysis = g4PSIAnalysisManager::getInstance();
            analysis->SetMuonDecay(w, x, y, z, p, pDir.x(), pDir.y(), pDir.z());
        }
    }

    
    
    
    // Info:
    // -----
    // http://geant4.web.cern.ch/geant4/support/faq.shtml#TRACK-1
    
//    int vol_index = g4PSIDetectorParts::getInstance()->TestTargetVolume(lvolume);
//    if (
//        //fElectromagnetic == post_proc->GetProcessType()
//        theTrack->GetTrackID() == 1   // this is the beam particle
//        && vol_index > 0  // process within the lH2 target or walls or heat shield
//        ) {
//        
//        G4ThreeVector p_in = pre_pt->GetMomentum();
//        G4ThreeVector p_out = post_pt->GetMomentum();
//        G4double theta = p_in.polarAngle(p_out);
//        
//        if (theta > 10*CLHEP::deg && p_out.mag() > 20 * CLHEP::MeV) {
//            // only record larger momentum scattering events to
//            // suppress the chance of having two such processes in one
//            // event
//            
//            // this is clunky, but we need to convert G4ThreeVectors to root's TVector3.
//            
//            Double_t w = theTrack->GetWeight();
//            TVector3 vertex(theTrack->GetPosition().x(), theTrack->GetPosition().y(), theTrack->GetPosition().z());
//            TVector3 mom_in(p_in.x(), p_in.y(), p_in.z());
//            TVector3 mom_out(p_out.x(), p_out.y(), p_out.z());
//            g4PSIAnalysisManager::getInstance()->SetTargetProcess(w, vertex, mom_in, mom_out, vol_index);
//
//            // For electrons, the relevant processes involved are "Transportation" and "eIoni" and "CoulombScat"
//            // CoulombScat is the process at larger angles
////            if (p_out.mag() > 30 * CLHEP::MeV && theta > 30*CLHEP::deg) {
////                std::cout << "p_in = " <<  p_in.mag() << " p_out = " << p_out.mag() << " theta = " << theta * 180/3.1416;
////                std::cout << " w  = " << w;
////                std::cout << " vertex = " << vertex.x() << " " << vertex.y() << " " << vertex.z();
////                std::cout << " process: " << post_proc->GetProcessName() << "\n";
////                std::cout << g4PSIAnalysisManager::getInstance()->GetEventSeed1() << " ";
////                std::cout << g4PSIAnalysisManager::getInstance()->GetEventSeed2() << "\n";
////                std::cout << " ------------- \n\n";
////            }
//        }
//    }
}
