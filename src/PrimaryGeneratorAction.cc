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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  fParticleGun->SetParticleEnergy(0*eV); 
  
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  //fParticleGun->SetParticleMomentumDirection(G4RandomDirection());          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
    //G4int Z = 10, A = 24;
    G4int Z = 92, A = 235;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  } 
  const G4double crystal_side = 102.0*mm;
  const G4double crystal_length = 203.0*mm;
  
  G4double x = G4UniformRand(); 
  G4double y = G4UniformRand(); 
  G4double z = G4UniformRand();

  // The numbers generated below are exclusive of the boundaries; no particles will be generated on the surface of the crystal.
  G4double random_x = (crystal_side*x)*mm - crystal_side/2.0*mm;
  G4double random_y = (crystal_side*y)*mm - crystal_side/2.0*mm;
  G4double random_z = (crystal_length*z)*mm - crystal_length/2.0*mm;

  // In aluminum
//  G4double random_x = (crystal_side)*mm - crystal_side/2.0*mm;
//  G4double random_y = (crystal_side*y)*mm - crystal_side/2.0*mm;
//  G4double random_z = (crystal_length*z)*mm - crystal_length/2.0*mm;
  
  //fParticleGun->SetParticlePosition(G4ThreeVector(random_x*mm, random_y*mm, random_z*mm));
  fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 140 * mm));
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());  
  fParticleGun->SetParticleEnergy(0*eV); 

//  fParticleGun->SetParticlePosition(G4ThreeVector(100*mm, 100*mm, 100*mm));

  //create vertex
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
