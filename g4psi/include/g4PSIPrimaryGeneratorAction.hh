#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "g4PSIPrimaryGeneratorMessenger.hh"
#include "globals.hh"

#include <fstream>
#include <iostream>

#include "G4VUserPrimaryGeneratorAction.hh"
#include <G4ParticleGun.hh>

class g4PSIPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    
    enum TBeamMixMode {
        DEFAULT_PARTICLE,
        BEAM_MIX_PLUS,
        BEAM_MIX_MINUS
    };
    
    enum TGunMode {
        DEFAULT_GUN,
        WALLUNIFORM,
        WALLUNIFORM_MOMDIST,
        WALLUNIFORM_FROM_CENTER,
        WALLCENTER,
        WALLCENTER_MOMDIST, /*!< shoot into the TOF-wall center with flat momentum distribution */
        BEAMLINE_2014,
        BEAMLINE_2015,
        BEAMLINE_2015_CENTER,
        BEAMLINE_PENCIL,
        BEAMLINE_PENCIL_FLAT,
        BEAMLINE_BOX,
        TARGET_MUON_DECAY,
        FULL4PI,   // 4pi from within target
        FULL4PI_FROM_CENTER,  // 4pi from (0,0,0)
        FULL4PI_DEFAULT,  // 4pi from gun position
        FULL4PI_MOMDIST,
        TURTLE,
        COSMIC /*!< comment */
    };
    
    g4PSIPrimaryGeneratorAction();
    ~g4PSIPrimaryGeneratorAction();
    
public:
    void GeneratePrimaries(G4Event* anEvent);
    void SetGunMode(TGunMode mode) {gunMode_ = mode;};
    void SetBeamMixMode(TBeamMixMode mode) {beamMixMode_ = mode;};
    void SetNPrimaryParticles(int n)
    {max_nPrimaryParticle_ = n;};
    
    void SetTurtleFile(G4String fn);
    void SetEventSeeds(G4long s1, G4long s2, bool ls, std::string fn);
    void SetBeamMomentumCenter(G4double p);
    void SetBeamMomentumSpread(G4double dpp);
    void SetPathLengthForRF(G4double path_length) {path0_ = path_length;};
    void SetBeamlineVertexZ(G4double vertexz) {z0_beamline = vertexz;};
    void SetRFJitter(G4double time) {proton_jitter = time;};
    G4double GetBeamMomentum();
    void SetTargetRadius(G4double r);
    
private:
    G4ParticleGun* particleGun_;
    g4PSIPrimaryGeneratorMessenger* gunMessenger;  // messenger of this class
    void GetBeamParticle(G4double &q, G4double &m, G4int &id);
    G4double GetRFTime(G4double m, G4double p, G4double z);
    G4double get_particle_momentum(G4double p0, G4double m, G4double theta);
    TGunMode gunMode_;
    TBeamMixMode beamMixMode_;
    G4int max_nPrimaryParticle_;
    G4double beamMomentumCenter_;
    G4double beamRelativeMomentumSpread_;
    G4double path0_;
    G4double z0_beamline;
    G4double proton_jitter;
    G4double beam_box_x0_;
    G4double beam_box_y0_;
    G4double beam_box_dx_;
    G4double beam_box_dy_;
    G4double target_radius_;
    
    long event_seed1_;
    long event_seed2_;
    bool load_seeds_;
    bool close_seed_file_;
    
    std::ifstream seed_file_;
    std::string seed_filename_;
    std::ifstream turtle_file_;
    std::string turtle_filename_;
};

#endif
