#ifndef g4PSITargetGrid_h
#define g4PSITargetGrid_h 1

#include "globals.hh"
#include "g4PSITargetBase.hh"
#include "g4PSITargetSD.hh"


class g4PSITargetGrid : public g4PSITargetBase
{
public:
    
    g4PSITargetGrid(G4String label, g4PSIChamberBase* chamber = NULL);
    
    void Placement();
    double GetZPos() {return 0;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:
    
    void init();
    
    G4LogicalVolume* target_log_;
    G4LogicalVolume* target_v_log_;
    G4LogicalVolume* target_h_log_;

    
    g4PSITargetSD* target_sd_;
};

#endif
