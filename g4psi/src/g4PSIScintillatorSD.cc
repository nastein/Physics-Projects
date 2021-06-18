#include "g4PSIScintillatorSD.hh"
#include "g4PSIScintillatorHit.hh"

#include "g4PSIAnalysisManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"

using namespace CLHEP;

g4PSIScintillatorSD::g4PSIScintillatorSD( G4String name,
                                         G4int nCells,
                                         G4String colName )
: G4VSensitiveDetector(name),
numberOfCells(nCells),
HCID(-1)
{
    G4String HCname;
    collectionName.insert(HCname=colName);
    CellID = new G4int[numberOfCells];
    SCHit_ = NULL;
}
g4PSIScintillatorSD::g4PSIScintillatorSD( G4String name,
                                         G4int nCells,
                                         G4String colName,
                                         G4bool collect_PID )
: G4VSensitiveDetector(name),
numberOfCells(nCells),
HCID(-1)
{
    G4String HCname;
    collectionName.insert(HCname=colName);
    CellID = new G4int[numberOfCells];
    SCHit_ = NULL;
    collect_PID_ = collect_PID;
}

g4PSIScintillatorSD::~g4PSIScintillatorSD()
{
    delete [] CellID;
}

void g4PSIScintillatorSD::Initialize(G4HCofThisEvent*)
{
    CalCollection = new g4PSIScintillatorHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    for(G4int j=0;j<numberOfCells;j++)
    {
        CellID[j] = -1;
    }
    
    hit_ = false;
    if (SCHit_) {
        *(SCHit_) = 0;
        SCCopyID_->clear();
        SCParticle_->clear();
        SCW_->clear();
        SCEdep_->clear();
        SCTime_->clear();
        SCPosX_->clear();
        SCPosY_->clear();
        SCPosZ_->clear();
        SCPosHitX_->clear();
        SCPosHitY_->clear();
        SCPosHitZ_->clear();
        SCDirHitX_->clear();
        SCDirHitY_->clear();
        SCDirHitZ_->clear();
        SCPosOriginX_->clear();
        SCPosOriginY_->clear();
        SCPosOriginZ_->clear();
        if (collect_PID_) {
            SCMuonEdep_->clear();
            SCElectronEdep_->clear();
            SCPhotonEdep_->clear();
        }
        //        SCProcess_->clear();
    } else
        G4Exception("g4PSIScintillatorSD::Initialize",
                    ("InitTree has not been called..." + collectionName[0]).c_str(),
                    FatalException, "");
    
}

G4bool g4PSIScintillatorSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    verboseLevel = 0;
    
    //    G4bool get_process = true;
    
    G4Track * theTrack = aStep->GetTrack();
    G4int pid = theTrack->GetDefinition()->GetPDGEncoding();
    
    G4double edep = aStep->GetTotalEnergyDeposit();
    if(edep<=0. && pid!=0) return false;
    
    G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    const G4VPhysicalVolume* physVol = hist->GetVolume();
    G4int copyID = hist->GetReplicaNumber();
    G4double hitTime = aStep->GetTrack()->GetGlobalTime();
    // they are the same: std::cout << copyID << " " << hist->GetVolume()->GetCopyNo() << "\n";

    if(CellID[copyID]==-1)
    {
        //std::cout << SensitiveDetectorName << " " << copyID << " ---- W = " << aStep->GetPreStepPoint()->GetWeight() << "\n";
        hit_ = true;  // Detector got hit
        g4PSIScintillatorHit* calHit = new g4PSIScintillatorHit(physVol->GetLogicalVolume());
        calHit->SetCopyID( copyID );  // the scintillator bar IDs are: 0, 1, 2, ...
        calHit->SetWeight(aStep->GetPreStepPoint()->GetWeight());
        calHit->SetEdep(edep);
        calHit->SetGlobalTime(hitTime);
        calHit->SetPosHit(aStep->GetPreStepPoint()->GetPosition());
        calHit->SetDirHit(aStep->GetPreStepPoint()->GetMomentumDirection());
        calHit->SetPosOrigin(aStep->GetTrack()->GetVertexPosition());
        calHit->SetParticleID(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());

        if (collect_PID_) {
            if(abs(pid) == 11)calHit->SetElectronEdep(edep);
            if(abs(pid) == 13)calHit->SetMuonEdep(edep);
            if(abs(pid) == 22)calHit->SetPhotonEdep(edep);
        }
        
//        if (get_process) {
//            const G4VProcess* creatorprocess = aStep->GetTrack()->GetCreatorProcess();
//            G4String process = "None";
//            if (creatorprocess != 0) process = creatorprocess->GetProcessName();
//            calHit->SetCreatorProcess(GetProcessCode(process));
//        }
        
        G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
        aTrans.Invert();
        calHit->SetPos(aTrans.NetTranslation());   // this is the position of the hit scintillator
        calHit->SetRot((aTrans.NetRotation()).inverse());
        
        G4int icell = CalCollection->insert(calHit);
        //G4cout << "copyID = " << copyID << ", icell = " << icell << " , so the CellID[copyID] = " << icell - 1 << G4endl;
        CellID[copyID] = icell - 1;
        if(verboseLevel>0) {
            G4cout << " New scintillator hit "<< SensitiveDetectorName << " in bar " << copyID << G4endl;
            G4cout << "    PosHit " << aStep->GetPreStepPoint()->GetPosition() << G4endl;
            calHit->Print();
            calHit->Draw();
        }
    }
    else
    {
        // \todo
        // Check why weight can be >1.
        // In one example with /tracking/verbose 3, the output includes:
        // **PostStepDoIt (after all invocations):
        //        ++List of invoked processes
        //        1) Transportation
        //        2) msc
        //        3) biasWrapper(CoulombScat) (Forced)
        //        ++G4Step Information
        // The (Forced) makes this different from previous steps, the weight changed from 1 to 1.035..
        // the apparent process was "Transportation"
        //
//        if ((*CalCollection)[CellID[copyID]]->GetWeight() != aStep->GetPreStepPoint()->GetWeight()) {
//            G4cout << "This should not happen: " << SensitiveDetectorName << " " << copyID << " W = " << aStep->GetPreStepPoint()->GetWeight() << " " << (*CalCollection)[CellID[copyID]]->GetWeight() << " seed: " << g4PSIAnalysisManager::getInstance()->GetEventSeed1() << " " << g4PSIAnalysisManager::getInstance()->GetEventSeed2() << " " <<"\n";
//        }
        
        // if already hit, then take the earlier time
        if ((*CalCollection)[CellID[copyID]]->GetGlobalTime() > hitTime) {
            (*CalCollection)[CellID[copyID]]->SetGlobalTime(hitTime);
        }
        (*CalCollection)[CellID[copyID]]->AddEdep(edep);
        if (collect_PID_) {
            if(abs(pid) == 11)(*CalCollection)[CellID[copyID]]->AddElectronEdep(edep);
            if(abs(pid) == 13)(*CalCollection)[CellID[copyID]]->AddMuonEdep(edep);
            if(abs(pid) == 22)(*CalCollection)[CellID[copyID]]->AddPhotonEdep(edep);
        }
        if(verboseLevel>0){;}
    }
    
    return true;
}

void g4PSIScintillatorSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
    HCE->AddHitsCollection( HCID, CalCollection );
    
    *(SCHit_) = CalCollection->entries();
    for (int n = 0; n < CalCollection->entries(); n++) {
        SCCopyID_->push_back( (*CalCollection)[n]->GetCopyID() );
        SCParticle_->push_back( (*CalCollection)[n]->GetParticleID() );
        SCTime_->push_back( (*CalCollection)[n]->GetGlobalTime()/ns );
        SCW_->push_back( (*CalCollection)[n]->GetWeight() );
        SCEdep_->push_back( (*CalCollection)[n]->GetEdep() );
        SCPosX_->push_back( ((*CalCollection)[n]->GetPos()).x() );
        SCPosY_->push_back( ((*CalCollection)[n]->GetPos()).y() );
        SCPosZ_->push_back( ((*CalCollection)[n]->GetPos()).z() );
        SCPosHitX_->push_back( ((*CalCollection)[n]->GetPosHit()).x() );
        SCPosHitY_->push_back( ((*CalCollection)[n]->GetPosHit()).y() );
        SCPosHitZ_->push_back( ((*CalCollection)[n]->GetPosHit()).z() );
        SCDirHitX_->push_back( ((*CalCollection)[n]->GetDirHit()).x() );
        SCDirHitY_->push_back( ((*CalCollection)[n]->GetDirHit()).y() );
        SCDirHitZ_->push_back( ((*CalCollection)[n]->GetDirHit()).z() );
        SCPosOriginX_->push_back( ((*CalCollection)[n]->GetPosOrigin()).x() );
        SCPosOriginY_->push_back( ((*CalCollection)[n]->GetPosOrigin()).y() );
        SCPosOriginZ_->push_back( ((*CalCollection)[n]->GetPosOrigin()).z() );
        if (collect_PID_) {
            SCMuonEdep_->push_back((*CalCollection)[n]->GetMuonEdep());
            SCElectronEdep_->push_back((*CalCollection)[n]->GetElectronEdep());
            SCPhotonEdep_->push_back((*CalCollection)[n]->GetPhotonEdep());
        }
        //        SCProcess_->push_back( ((*CalCollection)[n]->GetCreatorProcess()) );
    };
    
}

void g4PSIScintillatorSD::clear()
{
}

void g4PSIScintillatorSD::DrawAll()
{
}

void g4PSIScintillatorSD::PrintAll()
{
}

G4int g4PSIScintillatorSD::GetProcessCode(G4String process_name)
{
    G4int process_code = -1;
    if (process_name == "None") process_code = 0;
    if (process_name == "eBrem" || process_name == "muBrem") process_code = 1;
    if (process_name == "eIoni" || process_name == "muIoni" || process_name == "hIoni") process_code = 2;
    if (process_name == "eMuls" || process_name == "muMuls" || process_name == "hMuls") process_code = 3;
    if (process_name == "Decay") process_code = 4;
    if (process_name == "MuPair") process_code = 5;
    if (process_name == "compt") process_code = 6;
    if (process_name == "conv") process_code = 7;
    if (process_name == "annihil") process_code = 8;
    if (process_name == "phot") process_code = 9;
    if (process_name == "muMinusCaptureAtRest") process_code = 10;
    if (process_name == "hadElastic") process_code = 11;
    if (process_name == "NeutronInelastic" || process_name == "neutronInelastic") process_code = 12;
    if (process_name == "nCapture") process_code = 13;
    if (process_name == "PhotonInelastic") process_code = 14;
    if (process_name == "ElectroNuclear") process_code = 15;
    if (process_name == "ProtonInelastic" || process_name == "protonInelastic") process_code = 16;
    if (process_name == "dInelastic") process_code = 17;
    if (process_name == "tInelastic") process_code = 18;
    if (process_name == "PionMinusInelastic" || process_name == "pi-Inelastic") process_code = 19;
    if (process_name == "CHIPSNuclearCaptureAtRest") process_code = 20;
    if (process_name == "PositronNuclear") process_code = 21;
    if (process_name == "PionPlusInelastic" || process_name == "pi+Inelastic") process_code = 22;
    if (process_name == "hBertiniCaptureAtRest") process_code = 23;
    if (process_name == "photonNuclear") process_code = 24;
    if (process_name == "electronNuclear") process_code = 25;
    if (process_name == "positronNuclear") process_code = 26;
    
    // If process_code is still -1, print out warning so it can be added to this list...
    if (process_code == -1) {
        G4cout << "*==*==* WARNING: Process name " << process_name << " not found in list." << G4endl;
    }
    
    return process_code;
}

void g4PSIScintillatorSD::InitTree(TTree *T) {
    
    SCHit_ = new int();
    SCW_ = new VectorDouble();
    SCCopyID_ = new VectorInt();
    SCParticle_ = new VectorInt();
    SCEdep_ = new VectorDouble();
    SCTime_ = new VectorDouble();
    SCPosX_ = new VectorDouble();
    SCPosY_ = new VectorDouble();
    SCPosZ_ = new VectorDouble();
    SCPosHitX_ = new VectorDouble();
    SCPosHitY_ = new VectorDouble();
    SCPosHitZ_ = new VectorDouble();
    SCDirHitX_ = new VectorDouble();
    SCDirHitY_ = new VectorDouble();
    SCDirHitZ_ = new VectorDouble();
    SCPosOriginX_ = new VectorDouble();
    SCPosOriginY_ = new VectorDouble();
    SCPosOriginZ_ = new VectorDouble();

    if (collect_PID_) {
        SCMuonEdep_ = new VectorDouble();
        SCElectronEdep_ = new VectorDouble();
        SCPhotonEdep_ = new VectorDouble();
    }
    //    SCProcess_ = new VectorDouble();
    
    G4String label = GetName();
    T->Branch((label+"_Hit").c_str(), SCHit_);
    T->Branch((label+"_CopyID").c_str(), SCCopyID_);
    T->Branch((label+"_W").c_str(), SCW_);
    T->Branch((label+"_ParticleID").c_str(), SCParticle_);
    T->Branch((label+"_Edep").c_str(), SCEdep_);
    T->Branch((label+"_Time").c_str(), SCTime_);
    T->Branch((label+"_PosX").c_str(), SCPosX_);
    T->Branch((label+"_PosY").c_str(), SCPosY_);
    T->Branch((label+"_PosZ").c_str(), SCPosZ_);
    T->Branch((label+"_PosHitX").c_str(), SCPosHitX_);
    T->Branch((label+"_PosHitY").c_str(), SCPosHitY_);
    T->Branch((label+"_PosHitZ").c_str(), SCPosHitZ_);
    T->Branch((label+"_DirHitX").c_str(), SCDirHitX_);
    T->Branch((label+"_DirHitY").c_str(), SCDirHitY_);
    T->Branch((label+"_DirHitZ").c_str(), SCDirHitZ_);
    T->Branch((label+"_PosOriX").c_str(), SCPosOriginX_);
    T->Branch((label+"_PosOriY").c_str(), SCPosOriginY_);
    T->Branch((label+"_PosOriZ").c_str(), SCPosOriginZ_);
    if (collect_PID_) {
        T->Branch((label+"_MuonEdep").c_str(), SCMuonEdep_);
        T->Branch((label+"_ElectronEdep").c_str(), SCElectronEdep_);
        T->Branch((label+"_PhotonEdep").c_str(), SCPhotonEdep_);
    }
    //    T->Branch((label+"_Process").c_str(), SCProcess_);
}

void g4PSIScintillatorSD::DeleteEventData() {
    
    delete SCHit_;
    delete SCCopyID_;
    delete SCW_;
    delete SCParticle_;
    delete SCEdep_;
    delete SCTime_;
    delete SCPosX_;
    delete SCPosY_;
    delete SCPosZ_;
    delete SCPosHitX_;
    delete SCPosHitY_;
    delete SCPosHitZ_;
    delete SCDirHitX_;
    delete SCDirHitY_;
    delete SCDirHitZ_;
    delete SCPosOriginX_;
    delete SCPosOriginY_;
    delete SCPosOriginZ_;

    if (collect_PID_) {
        delete SCMuonEdep_;
        delete SCElectronEdep_;
        delete SCPhotonEdep_;
    }
    //    delete SCProcess_;
}
