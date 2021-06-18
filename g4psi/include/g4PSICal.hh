#ifndef g4PSICal_h
#define g4PSICal_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSICal
//!
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "g4PSIDetectorBase.hh"
#include "g4PSIScintillatorSD.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"


class g4PSICal : public g4PSIDetectorBase {
    
public:
    
    enum CalType {
        kPbGlass,
        kPbSC
    };
    
    g4PSICal(G4String label, CalType ct, G4double z, G4double t_pb, G4double t_sc, G4int n);
    ~g4PSICal();
    
public:
    void Placement();
    double GetZPos();
    void SetSD(G4SDManager *SDman);
    G4String GetName() {return Cal_label_;};
    void Write();
    void Info();
    
    void InitTree(TTree *T);
    void DeleteEventData();
    void ReadEventData(G4HCofThisEvent *HCE);
    
private:
    CalType cal_type_;
    G4String Cal_label_;
    G4double Cal_z_;
    G4int cal_nx_;
    G4int cal_ny_;
    G4int cal_nz_;
    G4Material *cal_mat_;
    G4double cal_t_pb_;
    G4double cal_t_sc_;
    G4double cal_dx_;
    G4double cal_dy_;
    
    G4LogicalVolume *cal_sc_log_;
    G4LogicalVolume *cal_pb_log_;
    g4PSIScintillatorSD* SD_;
};

#endif
