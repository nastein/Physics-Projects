//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4CoulombScatteringMUSE.cc 85246 2014-10-27 08:26:11Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CoulombScatteringMUSE
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 22.08.2004
//
// Modifications:
// 01.08.06 V.Ivanchenko add choice between G4eCoulombScatteringModel and
//          G4CoulombScatteringModel
// 15.06.15 S.Strauch copied for MUSE
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4CoulombScatteringMUSE.hh"
#include "G4SystemOfUnits.hh"
#include "G4eCoulombScatteringModelMUSE.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4Proton.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScatteringMUSE::G4CoulombScatteringMUSE(const G4String& name)
  : G4VEmProcess(name),q2Max(TeV*TeV),isInitialised(false)
{
  //  G4cout << "G4CoulombScatteringMUSE constructor "<< G4endl;
  SetBuildTableFlag(true);
  SetStartFromNullFlag(false);
  SetIntegral(true);
  SetSecondaryParticle(G4Proton::Proton());
  SetProcessSubType(fCoulombScattering);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScatteringMUSE::~G4CoulombScatteringMUSE()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4CoulombScatteringMUSE::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScatteringMUSE::InitialiseProcess(const G4ParticleDefinition* p)
{
  // second initialisation not allowed for the time being
  // this means that polar angle limit change will not be appled 
  // after first initialisation
  if(isInitialised) { return; }

  G4EmParameters* param = G4EmParameters::Instance();
  G4double a = param->FactorForAngleLimit()*CLHEP::hbarc/CLHEP::fermi;
  q2Max = 0.5*a*a;
  G4double theta = param->MscThetaLimit();

  // restricted or non-restricted cross section table
  G4bool yes = false;
  if(theta == CLHEP::pi) { yes = true; }
  SetStartFromNullFlag(yes);
  /*
  G4cout << "### G4CoulombScatteringMUSE::InitialiseProcess: "
  	 << p->GetParticleName()
	 << " Emin(MeV)= " << MinKinEnergy()/MeV
	 << " Emax(TeV)= " << MaxKinEnergy()/TeV
	 << " nbins= " << LambdaBinning()
	 << " theta= " << theta
	 << G4endl;
  */

  isInitialised = true;
  G4double mass = p->GetPDGMass();
  G4String name = p->GetParticleName();
  //G4cout << name << "  type: " << p->GetParticleType() 
  //<< " mass= " << mass << G4endl;
  yes = true;
  if (mass > GeV || p->GetParticleType() == "nucleus") {
    SetBuildTableFlag(false);
    yes = false;
    if(name != "GenericIon") { SetVerboseLevel(0); }
  } else {
    if(name != "e-" && name != "e+" &&
       name != "mu+" && name != "mu-" && name != "pi+" && 
       name != "kaon+" && name != "proton" ) { SetVerboseLevel(0); }
  }

  if(!EmModel(1)) { 
    if(yes) { SetEmModel(new G4eCoulombScatteringModelMUSE(), 1); }
    else    { SetEmModel(new G4IonCoulombScatteringModel(), 1); }
  }
  G4VEmModel* model = EmModel(1);
  G4double emin = std::max(param->MinKinEnergy(),model->LowEnergyLimit());
  G4double emax = std::min(param->MaxKinEnergy(),model->HighEnergyLimit());
  model->SetPolarAngleLimit(theta);
  model->SetLowEnergyLimit(emin);
  model->SetHighEnergyLimit(emax);
  AddEmModel(1, model);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoulombScatteringMUSE::MinPrimaryEnergy(const G4ParticleDefinition* part,
					       const G4Material* mat)
{
  // Pure Coulomb scattering
  G4double emin = 0.0;

  // Coulomb scattering combined with multiple or hadronic scattering
  G4double theta = G4EmParameters::Instance()->MscThetaLimit();

  if(0.0 < theta) {
    G4double p2 = q2Max*mat->GetIonisation()->GetInvA23()/(1.0 - cos(theta));
    G4double mass = part->GetPDGMass();
    emin = sqrt(p2 + mass*mass) - mass;
  }

  return emin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScatteringMUSE::PrintInfo()
{
  G4cout << "      " << G4EmParameters::Instance()->MscThetaLimit()/degree
	 << " < Theta(degree) < 180";

  if(q2Max < DBL_MAX) { G4cout << "; pLimit(GeV^1)= " << sqrt(q2Max)/GeV; }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
