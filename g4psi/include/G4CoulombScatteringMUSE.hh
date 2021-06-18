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
// $Id: G4CoulombScatteringMUSE.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4CoulombScatteringMUSE
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 12.03.2006
//
// Modifications:
//
// 15.06.2015 : S. Strauch, copied for MUSE
//
// Class Description:
//
// This class manages the process of Coulomb elastic scattering
//

// -------------------------------------------------------------------
//

#ifndef G4CoulombScatteringMUSE_h
#define G4CoulombScatteringMUSE_h 1

#include "G4VEmProcess.hh"
#include "G4VEmModel.hh"

class G4CoulombScatteringMUSE : public G4VEmProcess
{

public:

  G4CoulombScatteringMUSE(const G4String& name = "CoulombScat");

  virtual ~G4CoulombScatteringMUSE();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*);

private:

 // hide assignment operator
  G4CoulombScatteringMUSE & operator=(const G4CoulombScatteringMUSE &right);
  G4CoulombScatteringMUSE(const G4CoulombScatteringMUSE&);
  
  G4double q2Max;
  G4bool isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
