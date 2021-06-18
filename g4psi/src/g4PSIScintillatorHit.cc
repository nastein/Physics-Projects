#include "g4PSIScintillatorHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<g4PSIScintillatorHit> g4PSIScintillatorHitAllocator;

//g4PSIScintillatorHit::g4PSIScintillatorHit()
//{pLogV=NULL;}

g4PSIScintillatorHit::g4PSIScintillatorHit(G4LogicalVolume* logVol)
:copyID(0), weight(0), particleID(0), edep(0), Muonedep(0), Electronedep(0), Photonedep(0), pos(), pos_hit(),
dir_hit(), pos_o(), rot(), globalTime(), pLogV(logVol)
{;}

g4PSIScintillatorHit::~g4PSIScintillatorHit()
{;}

g4PSIScintillatorHit::g4PSIScintillatorHit(const g4PSIScintillatorHit &right)
: G4VHit()
{
    copyID = right.copyID;
    weight = right.weight;
    particleID = right.particleID;
    edep = right.edep;
    Muonedep = right.Muonedep;
    Electronedep = right.Electronedep;
    Photonedep = right.Photonedep;
    pos = right.pos;
    pos_hit = right.pos_hit;
    dir_hit = right.dir_hit;
    pos_o = right.pos_o;
    rot = right.rot;
    globalTime = right.globalTime;
    pLogV = right.pLogV;
}

const g4PSIScintillatorHit& g4PSIScintillatorHit::operator=(const g4PSIScintillatorHit &right)
{
    copyID = right.copyID;
    weight = right.weight;
    particleID = right.particleID;
    edep = right.edep;
    Muonedep = right.Muonedep;
    Electronedep = right.Electronedep;
    Photonedep = right.Photonedep;
    pos = right.pos;
    pos_hit = right.pos_hit;
    dir_hit = right.dir_hit;
    pos_o = right.pos_o;
    rot = right.rot;
    globalTime = right.globalTime;
    pLogV = right.pLogV;
    return *this;
}

G4int g4PSIScintillatorHit::operator==(const g4PSIScintillatorHit &right) const
{
    return (this==&right) ? 1 : 0;
}

void g4PSIScintillatorHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    
    if(pVVisManager)
    {
        G4Transform3D trans(rot,pos);
        G4VisAttributes attribs;
        const G4VisAttributes* pVA = pLogV->GetVisAttributes();
        if(pVA) attribs = *pVA;
        attribs.SetColour(G4Colour(1.,0.,0.));
        attribs.SetForceWireframe(false);
        attribs.SetForceSolid(true);
        pVVisManager->Draw(*pLogV,attribs,trans);
    }
}

void g4PSIScintillatorHit::Print()
{
    std::cout << "g4PSIScintillatorHit: edep = " << edep << " W = " << weight << "\n";
}


