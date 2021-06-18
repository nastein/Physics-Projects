#ifndef g4PSIScintillatorHit_h
#define g4PSIScintillatorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class g4PSIScintillatorHit : public G4VHit
{
public:
    
    //    g4PSIScintillatorHit();
    g4PSIScintillatorHit(G4LogicalVolume* logVol);
    ~g4PSIScintillatorHit();
    g4PSIScintillatorHit(const g4PSIScintillatorHit &right);
    const g4PSIScintillatorHit& operator=(const g4PSIScintillatorHit &right);
    G4int operator==(const g4PSIScintillatorHit &right) const;
    
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
    void operator delete(void *,void*){}
#endif
    
    void Draw();
    void Print();
    
private:
    G4int copyID;
    G4double weight;
    G4int particleID;
    G4double edep;
    G4double Muonedep;
    G4double Electronedep;
    G4double Photonedep;
    G4ThreeVector pos;
    G4ThreeVector pos_hit;
    G4ThreeVector dir_hit;
    G4ThreeVector pos_o;
    G4RotationMatrix rot;
    G4double globalTime;
    const G4LogicalVolume* pLogV;
    //    G4int process;
    
public:
    inline void SetEdep(G4double de) { edep = de; };
    inline void AddEdep(G4double de) { edep += de; };
    inline G4double GetEdep() { return edep; };

    inline void SetMuonEdep(G4double de) { Muonedep =de; };
    inline void AddMuonEdep(G4double de) { Muonedep += de; };
    inline G4double GetMuonEdep() { return Muonedep; };
    
    inline void SetElectronEdep(G4double de) { Electronedep =de; };
    inline void AddElectronEdep(G4double de) { Electronedep += de; };
    inline G4double GetElectronEdep() { return Electronedep; };
    
    inline void SetPhotonEdep(G4double de) { Photonedep =de; };
    inline void AddPhotonEdep(G4double de) { Photonedep += de; };
    inline G4double GetPhotonEdep() { return Photonedep; };
     
    inline void SetCopyID(G4int id) { copyID = id; };
    inline G4int GetCopyID() { return copyID; };

    inline void SetWeight(G4double w) { weight = w; };
    inline G4double GetWeight() { return weight; };

    inline void SetParticleID(G4int part) { particleID = part; };
    inline G4int GetParticleID() { return particleID; };
 
    inline void SetGlobalTime(G4double t) { globalTime = t; };
    inline G4double GetGlobalTime() { return globalTime; };
  
    inline void SetPos(G4ThreeVector xyz) { pos = xyz; };
    inline G4ThreeVector GetPos() { return pos; };
  
    inline void SetPosHit(G4ThreeVector xyz) { pos_hit = xyz; };
    inline G4ThreeVector GetPosHit() { return pos_hit; };

    inline void SetDirHit(G4ThreeVector xyz) { dir_hit = xyz; };
    inline G4ThreeVector GetDirHit() { return dir_hit; };
  
    inline void SetPosOrigin(G4ThreeVector xyz) { pos_o = xyz;};
    inline G4ThreeVector GetPosOrigin() { return pos_o; };
 
    inline void SetRot(G4RotationMatrix rmat) { rot = rmat; };
    inline G4RotationMatrix GetRot() { return rot; };
 
    inline const G4LogicalVolume * GetLogV() { return pLogV; };
    //    inline void SetCreatorProcess(G4int proc) { process = proc; };
    //    inline G4int GetCreatorProcess() { return process; };
    //}
};

typedef G4THitsCollection<g4PSIScintillatorHit> g4PSIScintillatorHitsCollection;

extern G4Allocator<g4PSIScintillatorHit> g4PSIScintillatorHitAllocator;

inline void* g4PSIScintillatorHit::operator new(size_t)
{

    void *aHit;
    aHit = (void *) g4PSIScintillatorHitAllocator.MallocSingle();
    return aHit;

}

inline void g4PSIScintillatorHit::operator delete(void *aHit)
{
    g4PSIScintillatorHitAllocator.FreeSingle((g4PSIScintillatorHit*) aHit);
}

#endif
