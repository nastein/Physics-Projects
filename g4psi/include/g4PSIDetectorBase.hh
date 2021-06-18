#ifndef g4PSIDetectorBase_h
#define g4PSIDetectorBase_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIDetectorBase
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include <TTree.h>
#include <iostream>
#include <fstream>

class g4PSIDetectorBase {
    
protected:
    g4PSIDetectorBase(G4String s);
    
public:
    virtual void Placement() = 0;
    virtual double GetZPos() {return 0;};
    virtual void SetSD(G4SDManager *SDman) = 0;
    virtual void Write() = 0;
    virtual void Info() = 0;
    virtual void InitTree(TTree *T) = 0;  // at the beginning of a run
    virtual void DeleteEventData() = 0;
    
    
    void Info(std::ofstream *file, int n);

    void SetMotherVolume(G4LogicalVolume *m) {mother_volume_ = m;};
    G4LogicalVolume* GetMotherLVolume() {return mother_volume_;};
    G4LogicalVolume* GetDetectorLVolume() {return detector_log_;};
    
    G4String GetDetectorName() {return label_;};
    G4int GetCollID() {return collID_;};
    void SetCollID(int ID) {collID_ = ID;};
    void SetHitRequired(G4bool t) {HitRequired_ = t;};
    G4bool HitRequired() {return HitRequired_;};
    
    virtual G4bool HasHit() {return true;}; /// No veto is the default for the trigger flag, overwrite in the derived classes, implement in the sensitive detectors
    
    void InfoTitle(G4String name);
    void InfoParDouble(G4String par, G4double val, G4String unit = "");
    void InfoPar3Double(G4String par, G4double val1, G4double val2, G4double val3, G4String unit = "");
    void InfoPar2Double(G4String par, G4double val1, G4double val2, G4String unit = "");
    void InfoParInt(G4String par, G4int val, G4String unit = "");
    void InfoParBool(G4String par, G4bool val);
    void InfoParString(G4String par, G4String val);
    void InfoEnd();
    
protected:
    G4String label_;
    G4LogicalVolume *log_;
    G4int collID_;
    
    G4LogicalVolume *mother_volume_;
    G4LogicalVolume *detector_log_;  // the detector log is the logical volume of the detector; this is the mother volume of subsequent detector parts.
    
    G4bool HitRequired_;
    G4bool info_wiki_;
    std::ofstream *info_wiki_file_;
    int info_wiki_n_;
};

#endif
