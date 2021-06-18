#ifndef g4PSITrackerHit_h
#define g4PSITrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class g4PSITrackerHit : public G4VHit
{
public:
    
    g4PSITrackerHit();
    g4PSITrackerHit(G4LogicalVolume* logVol);
    ~g4PSITrackerHit();
    g4PSITrackerHit(const g4PSITrackerHit &right);
    const g4PSITrackerHit& operator=(const g4PSITrackerHit &right);
    G4int operator==(const g4PSITrackerHit &right) const;
    
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
    void operator delete(void *,void*){}
#endif
    
    void Draw();
    void Print();
    
private:
    G4int copyID_;
    G4int particleID_;
    G4int trackID_;
    G4double edep_;
    G4ThreeVector pos_;
    G4ThreeVector dir_;
    G4double momentum_;
    G4double weight_;
    //  G4RotationMatrix rot;
    //  G4double globalTime;
    const G4LogicalVolume* pLogV_;
    
public:
    inline void SetEdep(G4double de) { edep_ = de; };
    inline G4double GetEdep() { return edep_; };
    
    inline void SetCopyID(G4int id) { copyID_ = id; };
    inline G4int GetCopyID() { return copyID_; };
    
    inline void SetTrackID(G4int id) { trackID_ = id; };
    inline G4int GetTrackID() { return trackID_; };
    
    inline void SetParticleID(G4int id) { particleID_ = id; };
    inline G4int GetParticleID() { return particleID_; };

    inline void SetPos(G4ThreeVector xyz) { pos_ = xyz; };
    inline G4ThreeVector GetPos() { return pos_; };
    
    inline void SetDir(G4ThreeVector xyz) { dir_ = xyz; };
    inline G4ThreeVector GetDir() { return dir_; };

    inline void SetMomentum(G4double p) { momentum_ = p; };
    inline G4double GetMomentum() { return momentum_; };

    inline void SetWeight(G4double w) { weight_ = w; };
    inline G4double GetWeight() { return weight_; };

    //  inline void SetRot(G4RotationMatrix rmat) { rot = rmat; };
    //  inline G4RotationMatrix GetRot() { return rot; };
    //  inline const G4LogicalVolume * GetLogV() { return pLogV; };
};

typedef G4THitsCollection<g4PSITrackerHit> g4PSITrackerHitsCollection;

extern G4Allocator<g4PSITrackerHit> g4PSITrackerHitAllocator;

inline void* g4PSITrackerHit::operator new(size_t)
{
    void *aHit;
    aHit = (void *) g4PSITrackerHitAllocator.MallocSingle();
    return aHit;
}

inline void g4PSITrackerHit::operator delete(void *aHit)
{
    g4PSITrackerHitAllocator.FreeSingle((g4PSITrackerHit*) aHit);
}

#endif
