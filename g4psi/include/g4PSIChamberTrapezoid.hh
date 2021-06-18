#ifndef g4PSIChamberTrapezoid_h
#define g4PSIChamberTrapezoid_h 1

#include "globals.hh"
#include "g4PSIChamberBase.hh"
#include "g4PSITargetSD.hh"
#include "g4PSIScintillatorSD.hh"


class g4PSIChamberTrapezoid : public g4PSIChamberBase
{
public:
    
    g4PSIChamberTrapezoid(G4String label);
    
    void Placement();
    double GetZPos() {return 0;};
    void SetSD(G4SDManager *SDman);
    void Write();
    void Info();
    void InitTree(TTree *T);
    void DeleteEventData();
    
private:

    void MakeWindow(G4String s,G4VSolid* &window, G4double w, G4double h, G4double t, G4double r);
    void MakeWindow(G4String s, G4LogicalVolume* &window_log, G4Material* mat, G4double w, G4double h, G4double t, G4double r);
    
    g4PSITargetSD* chamber_sd_;
    g4PSIScintillatorSD* veto_sd_;

    G4LogicalVolume* chamber_full_log_;
    G4LogicalVolume* chamber_full_vac_log_;
    G4LogicalVolume* chamber_entrance_window_log_;
    G4LogicalVolume* chamber_conflat_log_;
    G4LogicalVolume* chamber_conflat_vacuum_log_;
    G4LogicalVolume* chamber_exit_window_kapton_log_;
    G4LogicalVolume* chamber_side_exit_window_mylar_log_;
    G4LogicalVolume* chamber_side_exit_window_kevlar_log_;
    G4LogicalVolume* chamber_side_exit_window_poly_log_;
    G4LogicalVolume* chamber_veto_log_;
    G4LogicalVolume* chamber_side_veto_log_;
    G4LogicalVolume* chamber_downstream_window_reinforcement_log_;
    G4LogicalVolume* chamber_downstream_window_frame_log_;
    G4LogicalVolume* chamber_downstream_window_vacuum_log_;
    G4LogicalVolume* chamber_window_cover_log_;
    G4LogicalVolume* chamber_left_side_exit_frame_log_;
    G4LogicalVolume* chamber_right_side_exit_frame_log_;
    G4LogicalVolume* chamber_side_exit_frame_vacuum_log_;
};

#endif
