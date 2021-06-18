#ifndef g4PSITargetHit_h
#define g4PSITargetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class g4PSITargetHit : public G4VHit
{
public:
    
    g4PSITargetHit();
    //    g4PSITargetHit(G4LogicalVolume* logVol);
    ~g4PSITargetHit();
    g4PSITargetHit(const g4PSITargetHit &right);
    const g4PSITargetHit& operator=(const g4PSITargetHit &right);
    G4int operator==(const g4PSITargetHit &right) const;
    
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
    void operator delete(void *,void*){}
#endif
    
    void Draw();
    void Print();
    
public:
    
    inline void SetCopyID(G4int id) { copyID_ = id; };
    inline G4int GetCopyID() { return copyID_; };

    inline void SetParticleID(G4int id) { particleID_ = id; };
    inline G4int GetParticleID() { return particleID_; };
    
    inline void SetTrackID(G4int id) { trackID_ = id; };
    inline G4int GetTrackID() { return trackID_; };
    
    inline void SetProcessID(G4int id) { processID_ = id; };
    inline G4int GetProcessID() { return processID_; };

    inline void SetVertex(G4ThreeVector xyz) { vertex_ = xyz; };
    inline G4ThreeVector GetVertex() { return vertex_; };

    inline void SetWeight(G4double w) { weight_ = w; };
    inline G4double GetWeight() { return weight_; };

    inline void SetTheta(G4double x) { theta_ = x; };
    inline G4double GetTheta() { return theta_; };

    inline void SetMultipleScatterTheta(G4double x) { multiplescattertheta_ = x; };
    inline G4double GetMultipleScatterTheta() { return multiplescattertheta_; };

    inline void SetPIn(G4double x) { p_in_ = x; };
    inline G4double GetPIn() { return p_in_; };

    inline void SetPOut(G4double x) { p_out_ = x; };
    inline G4double GetPOut() { return p_out_; };

    inline void SetInScatter(G4bool x) { inscatter_= x; };
    inline G4bool GetInScatter() {return inscatter_; };
    
    
private:
    G4int copyID_;
    G4int particleID_;
    G4int trackID_;
    G4int processID_;
    G4ThreeVector vertex_;
    G4double weight_;
    G4double theta_;
    G4double multiplescattertheta_;
    G4double p_in_;
    G4double p_out_;
    G4bool inscatter_;
};

typedef G4THitsCollection<g4PSITargetHit> g4PSITargetHitsCollection;

extern G4Allocator<g4PSITargetHit> g4PSITargetHitAllocator;

inline void* g4PSITargetHit::operator new(size_t)
{

    void *aHit;
    aHit = (void *) g4PSITargetHitAllocator.MallocSingle();
    return aHit;
}

inline void g4PSITargetHit::operator delete(void *aHit)
{
    g4PSITargetHitAllocator.FreeSingle((g4PSITargetHit*) aHit);
}

#endif
