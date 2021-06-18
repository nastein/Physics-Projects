//
//  g4PSITestSphere.hh
//  g4PSI
//
//  Created by Steffen Strauch on 6/5/15.
//
//

#ifndef __g4PSI__g4PSITestSphere__
#define __g4PSI__g4PSITestSphere__

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSITestSphere
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSITrackerSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"

class g4PSITestSphere : public g4PSIDetectorBase {
    
public:
    g4PSITestSphere(G4String label,
                    G4double r,
                    G4double pth,
                    G4LogicalVolume *m = NULL);
    ~g4PSITestSphere();
    
public:
    void Placement();
    G4LogicalVolume* GetLog() {return testsphere_detector_log_;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    G4bool HasHit() {return SD_->SDHasHit();};

private:
    G4double r_;
    G4double p_threshold_;
    
    G4LogicalVolume *testsphere_detector_log_;
    G4LogicalVolume *testsphere_assembly_log_;
    g4PSITrackerSD* SD_;
    
};

#endif /* defined(__g4PSI__g4PSITestSphere__) */
