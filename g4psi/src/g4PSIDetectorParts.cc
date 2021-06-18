#include "g4PSIDetectorParts.hh"

#include <TVectorD.h>

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

g4PSIDetectorParts* g4PSIDetectorParts::fManager = 0;

g4PSIDetectorParts* g4PSIDetectorParts::getInstance() {
    if(!fManager) {
        fManager = new g4PSIDetectorParts();
    }
    return fManager;
}

void g4PSIDetectorParts::dispose() {
    delete fManager;
    fManager = NULL;
    
}

g4PSIDetectorParts::g4PSIDetectorParts() {
    detector_.clear();
    trigger_detectors_.clear();
    trigger_detectors_string_.clear();
    setup_ = kStandard2015;
    target_state_ = kFull;
    default_mother_ = NULL;
    target_ = NULL;
    target_wall_ = NULL;
    target_si_ = NULL;
    testMany_ = NULL;
    wiki_markup_filename_ = "";
    do_wiki_markup_ = false;
    userDouble1Par_ = 0;
    userDouble2Par_ = 0;
    userInt1Par_ = 0;
    userInt2Par_ = 0;
}

g4PSIDetectorParts::~g4PSIDetectorParts() {

    delete testMany_;
    
    /// \todo need to delete all detectors which are listed
}

g4PSIDetectorBase* g4PSIDetectorParts::GetDetector(G4String name) {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        
        /// \todo g4PSIDetectorBase needs to be written?
        
        if (detector_[cntr]->GetDetectorName().compare(name) == 0) {
            return detector_[cntr];
        }
    }
    G4Exception("g4PSIDetectorParts::GetDetector",
                ("Unknwon detector: " + name).c_str(),
                FatalException, "");
    return NULL;
}

void g4PSIDetectorParts::SetSD(G4SDManager *sdman) {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        detector_[cntr]->SetSD(sdman);
    };
}

void g4PSIDetectorParts::Info() {
    if (do_wiki_markup_) {
        int n = 0;
        wiki_markup_.open (wiki_markup_filename_);
        for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
            detector_[cntr]->Info(&wiki_markup_, ++n);
        };
        wiki_markup_.close();
    }
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        detector_[cntr]->Info();
    };
}

void g4PSIDetectorParts::Placement(G4LogicalVolume *m) {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        detector_[cntr]->SetMotherVolume(m);
        detector_[cntr]->Placement();
    };
}

void g4PSIDetectorParts::Placement() {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        if (!detector_[cntr]->GetMotherLVolume() ) detector_[cntr]->SetMotherVolume(default_mother_);
        detector_[cntr]->Placement();
    };
}

void g4PSIDetectorParts::AllOff() {
    for (int i = kShieldFloor; i != kTestPlanes; i++) {
        build[static_cast<BuildModes>(i)] = false;
    };
}

void g4PSIDetectorParts::RootInitTree(TTree *t) {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        detector_[cntr]->InitTree(t);
    }
}

void g4PSIDetectorParts::RootWrite() {
    
    TVectorD v1(kTestPlanes - kShieldFloor);
    for (int i = kShieldFloor; i != kTestPlanes; i++) {
        v1[i] = build[static_cast<BuildModes>(i)] ? 1.0 : 0.0;
    };
    v1.Write("DetectorParts");
    
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {

        detector_[cntr]->Write();
        detector_[cntr]->DeleteEventData();
    }
}

void g4PSIDetectorParts::DeleteEventData() {
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        //    std::cout << "g4PSIDetectorParts: DeleteEventData for: " << detector_[cntr]->GetName() << "\n";
        detector_[cntr]->DeleteEventData();
    }
}

void g4PSIDetectorParts::Print() {
    for (int i = kShieldFloor; i != kTestPlanes; i++) {
        G4cout << "Detector part: " << i << " is " << (build[static_cast<BuildModes>(i)] ? "on" : "off") << G4endl;
    };
    
    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    
    if (kill_primary_log_.size() > 0) {
        G4cout << "Warning: The primary particle will be killed if it enters any of the following logical volumes" << G4endl;
        for (size_t i = 0; i < kill_primary_log_.size(); i++) {
            G4LogicalVolume *l = kill_primary_log_[i];
            if (l) {
                G4cout << "(" << i << ") Logical volume : " << l->GetName() << G4endl;
            } else {
                G4cout << "(" << i << ") UNDEFINED VOLUME" << G4endl;
            }
        }
        G4cout << G4endl;
    }
}

void g4PSIDetectorParts::Add(G4LogicalVolume *m, g4PSIDetectorBase* d) {
    d->SetMotherVolume(m);
    Add(d);
}

void g4PSIDetectorParts::Add(g4PSIDetectorBase* d) {
    // Add a new detector to the list of detectors.
    // Signal an exception if a detector of that name already exists.
    for (unsigned cntr = 0; cntr < detector_.size(); cntr++) {
        if (detector_[cntr]->GetDetectorName() == d->GetDetectorName()) {
            G4Exception("g4PSIDetectorParts::Add()",
                        ("Error: multiple detector parts with the same name: " + d->GetDetectorName()).c_str(),
                        FatalException, "");
        };
    }
    detector_.push_back(d);
}

void g4PSIDetectorParts::BeginOfRun() {
    
    using namespace std;
    
    for (unsigned itrig = 0; itrig < trigger_detectors_string_.size(); itrig++) {
        
        detector_list set_of_trigger_detectors;
        
        istringstream iss(trigger_detectors_string_[itrig]);
        
        vector<G4String> tokens;
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter<vector<G4String> >(tokens));
        
        for (unsigned i = 0; i < tokens.size(); i++) {
            
            bool found_detector = false;
            for (unsigned j = 0; j < detector_.size(); j++) {
                if (detector_[j]->GetDetectorName() == tokens[i]) {
                    found_detector = true;
                    set_of_trigger_detectors.push_back(detector_[j]);
                }
            }
            if (!found_detector)
                G4Exception("g4PSIDetectorParts::BeginOfRun()",
                            ("Error: unknown setup: " + tokens[i]).c_str(),
                            FatalException, "");
            
        }
        trigger_detectors_.push_back(set_of_trigger_detectors);
    }
}

G4bool g4PSIDetectorParts::HasTrigger() {
    
    G4bool trig = trigger_detectors_.size() == 0 ? true : false;
    
    for (unsigned itrig = 0; itrig < trigger_detectors_.size(); itrig++) {
        
        G4bool this_trig = true;
        detector_list set_of_trigger_detectors = trigger_detectors_[itrig];
        for (unsigned cntr = 0; cntr < set_of_trigger_detectors.size(); cntr++) {
            this_trig = this_trig && set_of_trigger_detectors[cntr]->HasHit();
        }
        trig = trig || this_trig;
    }
    
    return trig;
}

void g4PSIDetectorParts::SetTriggerDetectors(G4String name) {
    std::cout << "New Trigger: `" << name << "'\n";
    trigger_detectors_string_.push_back(name);
}


void g4PSIDetectorParts::SetTargetState(TargetState state) {
    target_state_ = state;
}


g4PSIDetectorParts::TargetState g4PSIDetectorParts::GetTargetState() {
    return target_state_;
}


void g4PSIDetectorParts::AttachBiasingOperator(G4LogicalVolume *lv)
{
    if (lv) {
        if (!testMany_) {
            testMany_ = new GB01BOptrMultiParticleChangeCrossSection();
            testMany_->AddParticle("e-");
            testMany_->AddParticle("e+");
            testMany_->AddParticle("mu-");
            testMany_->AddParticle("mu+");   
        }
        testMany_->AttachTo(lv);
        G4cout << "*****  Attaching biasing operator " << testMany_->GetName()
        << " to logical volume \"" << lv->GetName() << "\""
        << G4endl;
    }
}


void g4PSIDetectorParts::SetWikiMarkupFile(G4String s) {
    do_wiki_markup_ = true;
    wiki_markup_filename_ = s;
}


void g4PSIDetectorParts::AddKillPrimaryLog(G4LogicalVolume *lvolume) {
    kill_primary_log_.push_back(lvolume);
}

void g4PSIDetectorParts::SetUserD1Par(G4double p) {
    userDouble1Par_ = p;
}

void g4PSIDetectorParts::SetUserD2Par(G4double p) {
    userDouble2Par_ = p;
}

void g4PSIDetectorParts::SetUserI1Par(G4int p) {
    userInt1Par_ = p;
}

void g4PSIDetectorParts::SetUserI2Par(G4int p){
    userInt2Par_ = p;
}

void g4PSIDetectorParts::SetUserT1Par(G4double p) {
    userDoubleT1Par_ = p;
}

void g4PSIDetectorParts::SetUserT2Par(G4double p) {
    userDoubleT2Par_ = p;
}

void g4PSIDetectorParts::SetUserT3Par(G4double p) {
    userDoubleT3Par_ = p;
}

void g4PSIDetectorParts::SetUserT4Par(G4double p) {
    userDoubleT4Par_ = p;
}

void g4PSIDetectorParts::SetUserT5Par(G4double p) {
    userDoubleT5Par_ = p;
}



