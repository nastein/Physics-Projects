#include "g4PSITargetHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<g4PSITargetHit> g4PSITargetHitAllocator;

g4PSITargetHit::g4PSITargetHit()
{;}

//g4PSITargetHit::g4PSITargetHit(G4LogicalVolume* logVol)
//:pLogV_(logVol)
//{;}

g4PSITargetHit::~g4PSITargetHit()
{;}

g4PSITargetHit::g4PSITargetHit(const g4PSITargetHit &right)
: G4VHit()
{
    copyID_ = right.copyID_;
    particleID_ = right.particleID_;
    trackID_ = right.trackID_;
    processID_ = right.processID_;
    vertex_ = right.vertex_;
    weight_ = right.weight_;
    theta_ = right.theta_;
    multiplescattertheta_ = right.multiplescattertheta_;
    p_in_ = right.p_in_;
    p_out_ = right.p_out_;
    inscatter_ = right.inscatter_;
}

const g4PSITargetHit& g4PSITargetHit::operator=(const g4PSITargetHit &right)
{
    copyID_ = right.copyID_;
    particleID_ = right.particleID_;
    trackID_ = right.trackID_;
    processID_ = right.processID_;
    vertex_ = right.vertex_;
    weight_ = right.weight_;
    theta_ = right.theta_;
    multiplescattertheta_ = right.multiplescattertheta_;
    p_in_ = right.p_in_;
    p_out_ = right.p_out_;
    inscatter_ = right.inscatter_;
    return *this;
}

G4int g4PSITargetHit::operator==(const g4PSITargetHit &right) const
{
    return (this==&right) ? 1 : 0;
}

void g4PSITargetHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(vertex_);
        circle.SetScreenSize(0.04);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(0.,0.,1.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

void g4PSITargetHit::Print()
{
    std::cout << "g4PSITargetHit: vertex = " << vertex_ << "\n";
}


