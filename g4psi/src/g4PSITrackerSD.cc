#include "g4PSITrackerSD.hh"
#include "g4PSITrackerHit.hh"
#include "g4PSIAnalysisManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

g4PSITrackerSD::g4PSITrackerSD( G4String name,
                               G4String colName,
                               G4bool isTP,
                               G4double pth)
: G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname=colName);
    HCID_ = -1;
    TrackHit_ = NULL;
    isTestPlane_ = isTP;
    p_threshold_ = pth;
    depth_ = 0;
}

g4PSITrackerSD::~g4PSITrackerSD() {;}

void g4PSITrackerSD::Initialize(G4HCofThisEvent* HCE)
{
    trackerCollection_ = new g4PSITrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if(HCID_<0) { HCID_ = GetCollectionID(0); }
    HCE->AddHitsCollection(HCID_,trackerCollection_);
    
    hit_ = false;
    if (TrackHit_) {
        *(TrackHit_) = 0;
        TrackParticleID_->clear();
        TrackTrackID_->clear();
        TrackCopyID_->clear();
        Trackw_->clear();
        Trackp_->clear();
        TrackEdep_->clear();
        TrackPosHitX_->clear();
        TrackPosHitY_->clear();
        TrackPosHitZ_->clear();
        TrackDirHitX_->clear();
        TrackDirHitY_->clear();
        TrackDirHitZ_->clear();
    } else
        G4Exception("g4PSITrackerSD::Initialize",
                    "InitTree has not been called...",
                    FatalException, "");
}

G4bool g4PSITrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    //        G4cout << "Debug: " << SensitiveDetectorName << " edep = " << edep << " ID=" << aStep->GetTrack()->GetDefinition()->GetPDGEncoding() << " track=" << aStep->GetTrack()->GetTrackID() << " step=" << aStep->GetPreStepPoint()->GetStepStatus() << "\n";
    G4StepPoint* preStep = aStep->GetPreStepPoint();
    
    G4double p = (preStep->GetMomentum()).mag();
    
    G4Track * theTrack = aStep->GetTrack();
    G4int pid = theTrack->GetDefinition()->GetPDGEncoding();
    
    // see http://geant4.slac.stanford.edu/tutorial/material06/HandsOn5/HandsOn5.html
    
    if ( pid == 0   // geantino events
        || edep != 0  // events with energy deposition
        || (isTestPlane_ && p > p_threshold_)  // test plane events above threshold
        ){
        
        G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
        G4int copyID = touchable->GetReplicaNumber(depth_);
        
        //std::cout << SensitiveDetectorName << " ID=" << copyID << " (" << touchable->GetVolume()->GetName()<< "), Pos = [" << preStep->GetPosition().x() << ", " << preStep->GetPosition().y() << ", " << preStep->GetPosition().z() << "], PID=" << aStep->GetTrack()->GetDefinition()->GetPDGEncoding() << ", edep = " << edep << "\n";
        
        hit_ = true;  // Detector got hit
        g4PSITrackerHit* newHit = new g4PSITrackerHit();
        newHit->SetEdep( edep );
        newHit->SetCopyID( copyID );
        newHit->SetTrackID( aStep->GetTrack()->GetTrackID() );
        newHit->SetParticleID( pid );
        newHit->SetWeight( preStep->GetWeight() );
        newHit->SetPos( preStep->GetPosition() );
        newHit->SetDir( preStep->GetMomentumDirection() );
        newHit->SetMomentum ( p );
        trackerCollection_->insert( newHit );
        
        //debug !!!
        //        G4AffineTransform aTrans = touchable->GetHistory()->GetTopTransform();
        //        aTrans.Invert();
        //        std::cout <<  SensitiveDetectorName << ": " << preStep->GetPosition()  << "   " << aTrans.NetTranslation() << "\n";
        
        
        //        double x0 = aStep->GetPreStepPoint()->GetPosition().x();
        //        double y0 = aStep->GetPreStepPoint()->GetPosition().y();
        //        double z0 = aStep->GetPreStepPoint()->GetPosition().z();
        //        double vx = aStep->GetPreStepPoint()->GetMomentumDirection().x();
        //        double vy = aStep->GetPreStepPoint()->GetMomentumDirection().y();
        //        double vz = aStep->GetPreStepPoint()->GetMomentumDirection().z();
        //        double z = z0 - vz * (vx*x0 + vy*y0) / (vx*vx + vy*vy);
        //        if (z < -4 * CLHEP::mm) {
        //           std::cout << "debug: reconstructed vertex z = " << z << "\n";
        //        };
        
        if(verboseLevel>0) {
            G4cout << " New hit " << SensitiveDetectorName << " " << newHit->GetPos() << G4endl;
            newHit->Print();
            newHit->Draw();
        }
    }
    
    return true;
}

void g4PSITrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    if(HCID_<0)
    { HCID_ = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
    HCE->AddHitsCollection( HCID_, trackerCollection_ );
    
    *(TrackHit_) = trackerCollection_->entries();
    for (int n = 0; n < trackerCollection_->entries(); n++) {
        TrackCopyID_->push_back( (*trackerCollection_)[n]->GetCopyID() );
        TrackParticleID_->push_back( (*trackerCollection_)[n]->GetParticleID() );
        TrackTrackID_->push_back( (*trackerCollection_)[n]->GetTrackID() );
        Trackw_->push_back( (*trackerCollection_)[n]->GetWeight() );
        Trackp_->push_back( (*trackerCollection_)[n]->GetMomentum() );
        TrackEdep_->push_back( (*trackerCollection_)[n]->GetEdep() );
        TrackPosHitX_->push_back( ((*trackerCollection_)[n]->GetPos()).x() );
        TrackPosHitY_->push_back( ((*trackerCollection_)[n]->GetPos()).y() );
        TrackPosHitZ_->push_back( ((*trackerCollection_)[n]->GetPos()).z() );
        TrackDirHitX_->push_back( ((*trackerCollection_)[n]->GetDir()).x() );
        TrackDirHitY_->push_back( ((*trackerCollection_)[n]->GetDir()).y() );
        TrackDirHitZ_->push_back( ((*trackerCollection_)[n]->GetDir()).z() );
    };
}

void g4PSITrackerSD::clear()
{
}

void g4PSITrackerSD::DrawAll()
{
}

void g4PSITrackerSD::PrintAll()
{
}

void g4PSITrackerSD::InitTree(TTree *T, std::string suffix) {
    
    TrackCopyID_ = new VectorInt();
    TrackParticleID_ = new VectorInt();
    TrackTrackID_ = new VectorInt();
    Trackp_ = new VectorDouble();
    Trackw_ = new VectorDouble();
    TrackEdep_ = new VectorDouble();
    TrackPosHitX_ = new VectorDouble();
    TrackPosHitY_ = new VectorDouble();
    TrackPosHitZ_ = new VectorDouble();
    TrackDirHitX_ = new VectorDouble();
    TrackDirHitY_ = new VectorDouble();
    TrackDirHitZ_ = new VectorDouble();
    TrackHit_ = new G4int();
    
    G4String label = GetName() + suffix;
    T->Branch((label+"_Hit").c_str(), TrackHit_);
    T->Branch((label+"_CopyID").c_str(), TrackCopyID_);
    T->Branch((label+"_ParticleID").c_str(), TrackParticleID_);
    T->Branch((label+"_TrackID").c_str(), TrackTrackID_);
    T->Branch((label+"_w").c_str(), Trackw_);
    T->Branch((label+"_p").c_str(), Trackp_);
    T->Branch((label+"_Edep").c_str(), TrackEdep_);
    T->Branch((label+"_PosHitX").c_str(), TrackPosHitX_);
    T->Branch((label+"_PosHitY").c_str(), TrackPosHitY_);
    T->Branch((label+"_PosHitZ").c_str(), TrackPosHitZ_);
    T->Branch((label+"_DirHitX").c_str(), TrackDirHitX_);
    T->Branch((label+"_DirHitY").c_str(), TrackDirHitY_);
    T->Branch((label+"_DirHitZ").c_str(), TrackDirHitZ_);
}

void g4PSITrackerSD::DeleteEventData() {
    
    delete TrackHit_;
    delete TrackCopyID_;
    delete TrackParticleID_;
    delete TrackTrackID_;
    delete Trackw_;
    delete Trackp_;
    delete TrackPosHitX_;
    delete TrackPosHitY_;
    delete TrackPosHitZ_;
    delete TrackDirHitX_;
    delete TrackDirHitY_;
    delete TrackDirHitZ_;
}
