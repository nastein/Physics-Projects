#include "g4PSITargetSD.hh"
#include "g4PSITargetHit.hh"
#include "g4PSIAnalysisManager.hh"

#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

g4PSITargetSD::g4PSITargetSD(G4String name,
                             G4String colName,
                             G4double min_p,
                             G4double min_theta)
: G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname=colName);
    HCID_ = -1;
    TargetHit_ = NULL;
    min_p_ = min_p;
    min_theta_ = min_theta;
}

g4PSITargetSD::~g4PSITargetSD() {;}

void g4PSITargetSD::Initialize(G4HCofThisEvent* HCE)
{
    trackerCollection_ = new g4PSITargetHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if(HCID_<0) { HCID_ = GetCollectionID(0); }
    HCE->AddHitsCollection(HCID_,trackerCollection_);
    
    hit_ = false;
    if (TargetHit_) {
        *(TargetHit_) = 0;
        TargetCopyID_->clear();
        TargetParticleID_->clear();
        TargetTrackID_->clear();
        TargetProcessID_->clear();
        Targetw_->clear();
        TargetVertexX_->clear();
        TargetVertexY_->clear();
        TargetVertexZ_->clear();
        TargetTheta_->clear();
        MultipleScatterTheta_->clear();
        TargetPIn_->clear();
        TargetPOut_->clear();
        TargetInScatter_->clear();
    } else
        G4Exception("g4PSITargetSD::Initialize",
                    "InitTree has not been called...",
                    FatalException, "");
}

G4bool g4PSITargetSD::ProcessHits(G4Step* theStep, G4TouchableHistory*)
{
    G4Track * theTrack = theStep->GetTrack();
    
    G4StepPoint *preStep  = theStep->GetPreStepPoint(); //  the current step
    G4StepPoint *postStep  = theStep->GetPostStepPoint();
    
    G4ThreeVector p_in = preStep->GetMomentum();
    G4ThreeVector p_out = postStep->GetMomentum();
    double deltaP = (p_in - p_out).mag();
    G4double theta = acos(p_in.cosTheta(p_out));
    G4ProcessType pType = postStep->GetProcessDefinedStep()->GetProcessType();
    G4int pSubType = postStep->GetProcessDefinedStep()->GetProcessSubType();
    G4int processID = pType * 1000 + pSubType;
    G4int trackID = theTrack->GetTrackID();
    G4int pid = theTrack->GetDefinition()->GetPDGEncoding();
    G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
    G4int copyID = touchable->GetReplicaNumber();
    
    
    if ( pid == 0   // geantino
        || (theta > min_theta_ && p_out.mag() > min_p_ &&
            pType != G4ProcessType::fGeneral &&
            processID != 1091 // Transportation
            )
        || processID == 2001 // Coulomb
        || (processID == 2002 && deltaP > min_p_ && (abs(pid) == 11 || abs(pid) == 13)) // process 2002 = eIoni + hIoni
        || (processID == 2003 && deltaP > min_p_) // process 2003 = eBrem
        || (trackID == 1 && (copyID == 0 || copyID == 2 || copyID == 50|| copyID == 52|| copyID == 54 || copyID == 56)) //Looking at inscattering from chamber 
        ){

        hit_ = true;  // Detector got hit
        g4PSITargetHit* newHit = new g4PSITargetHit();
        newHit->SetCopyID( copyID );
        newHit->SetParticleID( pid );
        newHit->SetTrackID( trackID );
        newHit->SetProcessID(processID);
        newHit->SetWeight( theTrack->GetWeight() );
        newHit->SetVertex( preStep->GetPosition() );
        newHit->SetTheta ( theta );
        newHit->SetMultipleScatterTheta( theta );
        newHit->SetPIn( p_in.mag() );
        newHit->SetPOut( p_out.mag() );

        if (trackID == 1 && copyID == 10 && processID == 2001 && theta*(180/3.1415926) > 20.0 && theta*(180/3.1415926) < 100.0 ) {
            g4PSIAnalysisManager::getInstance()->is_signal(true);
        } 

        //record hits in the chamber body
        if (trackID == 1 && (copyID == 0 || copyID == 2 || copyID == 50 || copyID == 52 || copyID == 54 || copyID == 56)) { newHit->SetInScatter(true);}
        else (newHit->SetInScatter(false));
        
        trackerCollection_->insert( newHit );
        
        if(verboseLevel>0 && trackID == 1 && (copyID == 0 || copyID == 50 || copyID == 52 || copyID == 54 || copyID == 56)) {
            G4cout << " New hit " << SensitiveDetectorName << " " << newHit->GetVertex() << G4endl;
            G4cout << " CopyID = " << copyID << G4endl;
            newHit->Print();
            newHit->Draw();
        }
    }
    
    
    // http://geant4.in2p3.fr/2005/Workshop/ShortCourse/session3/M.Asai1.pdf
    // "... Currently, returning boolean value is not used."
    
    return true;
}

void g4PSITargetSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    if(HCID_<0)
    { HCID_ = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
    HCE->AddHitsCollection( HCID_, trackerCollection_ );
    
    *(TargetHit_) = trackerCollection_->entries();
    for (int n = 0; n < trackerCollection_->entries(); n++) {
        TargetCopyID_->push_back( (*trackerCollection_)[n]->GetCopyID() );
        TargetParticleID_->push_back( (*trackerCollection_)[n]->GetParticleID() );
        TargetTrackID_->push_back( (*trackerCollection_)[n]->GetTrackID() );
        TargetProcessID_->push_back( (*trackerCollection_)[n]->GetProcessID() );
        Targetw_->push_back( (*trackerCollection_)[n]->GetWeight() );
        TargetVertexX_->push_back( ((*trackerCollection_)[n]->GetVertex()).x() );
        TargetVertexY_->push_back( ((*trackerCollection_)[n]->GetVertex()).y() );
        TargetVertexZ_->push_back( ((*trackerCollection_)[n]->GetVertex()).z() );
        TargetTheta_->push_back( (*trackerCollection_)[n]->GetTheta() );
        MultipleScatterTheta_->push_back( (*trackerCollection_)[n]->GetMultipleScatterTheta() );
        TargetPIn_->push_back( (*trackerCollection_)[n]->GetPIn() );
        TargetPOut_->push_back( (*trackerCollection_)[n]->GetPOut() );
        TargetInScatter_->push_back( (*trackerCollection_)[n]->GetInScatter() );
    };
}

void g4PSITargetSD::clear()
{
}

void g4PSITargetSD::DrawAll()
{
}

void g4PSITargetSD::PrintAll()
{
}

void g4PSITargetSD::InitTree(TTree *T, std::string suffix) {
    
    TargetHit_ = new G4int();
    TargetCopyID_ = new VectorInt();
    TargetParticleID_ = new VectorInt();
    TargetTrackID_ = new VectorInt();
    TargetProcessID_ = new VectorInt();
    Targetw_ = new VectorDouble();
    TargetVertexX_ = new VectorDouble();
    TargetVertexY_ = new VectorDouble();
    TargetVertexZ_ = new VectorDouble();
    TargetTheta_ = new VectorDouble();
    MultipleScatterTheta_ = new VectorDouble();
    TargetPIn_ = new VectorDouble();
    TargetPOut_ = new VectorDouble();
    TargetInScatter_ = new VectorBool();
    
    G4String label = GetName() + suffix;
    T->Branch((label+"_Hit").c_str(), TargetHit_);
    T->Branch((label+"_CopyID").c_str(), TargetCopyID_);
    T->Branch((label+"_ParticleID").c_str(), TargetParticleID_);
    T->Branch((label+"_TrackID").c_str(), TargetTrackID_);
    T->Branch((label+"_ProcessID").c_str(), TargetProcessID_);
    T->Branch((label+"_w").c_str(), Targetw_);
    T->Branch((label+"_VertexX").c_str(), TargetVertexX_);
    T->Branch((label+"_VertexY").c_str(), TargetVertexY_);
    T->Branch((label+"_VertexZ").c_str(), TargetVertexZ_);
    T->Branch((label+"_Theta").c_str(), TargetTheta_);
    T->Branch((label+"_MultipleScatterTheta").c_str(), MultipleScatterTheta_);
    T->Branch((label+"_PIn").c_str(), TargetPIn_);
    T->Branch((label+"_Pout").c_str(), TargetPOut_);
    T->Branch((label+"_InScatter").c_str(), TargetInScatter_);
}

void g4PSITargetSD::DeleteEventData() {
    
    delete TargetHit_;
    delete TargetCopyID_;
    delete TargetParticleID_;
    delete TargetTrackID_;
    delete TargetProcessID_;
    delete Targetw_;
    delete TargetVertexX_;
    delete TargetVertexY_;
    delete TargetVertexZ_;
    delete TargetTheta_;
    delete MultipleScatterTheta_;
    delete TargetPIn_;
    delete TargetPOut_;
    delete TargetInScatter_;
}
