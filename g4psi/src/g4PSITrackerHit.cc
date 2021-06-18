#include "g4PSITrackerHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<g4PSITrackerHit> g4PSITrackerHitAllocator;

g4PSITrackerHit::g4PSITrackerHit()
{pLogV_=NULL;}

g4PSITrackerHit::g4PSITrackerHit(G4LogicalVolume* logVol)
:pLogV_(logVol)
{;}

g4PSITrackerHit::~g4PSITrackerHit()
{;}

g4PSITrackerHit::g4PSITrackerHit(const g4PSITrackerHit &right)
: G4VHit()
{
    edep_ = right.edep_;
    copyID_ = right.copyID_;
    trackID_ = right.trackID_;
    weight_ = right.weight_;
    pos_ = right.pos_;
    dir_ = right.dir_;
    particleID_ = right.particleID_;
    momentum_ = right.momentum_;
    // rot = right.rot;
    //    pLogV = right.pLogV;
}

const g4PSITrackerHit& g4PSITrackerHit::operator=(const g4PSITrackerHit &right)
{
    edep_ = right.edep_;
    copyID_ = right.copyID_;
    trackID_ = right.trackID_;
    weight_ = right.weight_;
    pos_ = right.pos_;
    dir_ = right.dir_;
    particleID_ = right.particleID_;
    momentum_ = right.momentum_;
    // rot = right.rot;
    //    pLogV = right.pLogV;
    return *this;
}

G4int g4PSITrackerHit::operator==(const g4PSITrackerHit &right) const
{
    return (this==&right) ? 1 : 0;
}

void g4PSITrackerHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(pos_);
        circle.SetScreenSize(0.04);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,0.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

void g4PSITrackerHit::Print()
{
    std::cout << "g4PSITrackerHit: edep = " << edep_ << "\n";
}


