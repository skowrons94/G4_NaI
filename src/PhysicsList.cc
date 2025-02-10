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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "G4PhysicsListHelper.hh"
#include "G4Radioactivation.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclideTable.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"

// Livermore
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4NuclideTable.hh"
//#include "BiasedRDPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VUserPhysicsList()
{
  //add new units for radioActive decays
  // 
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);

  // mandatory for G4NuclideTable
  //
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  //read new PhotonEvaporation data set 
  //
  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(false);
  deex->SetStoreAllLevels(true);
  deex->SetIsomerProduction(true);  
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                /std::log(2.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  
  // baryons
  G4Proton::ProtonDefinition();
  G4Neutron::NeutronDefinition();  

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();  
  
  ////
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  G4Radioactivation* radioactiveDecay = new G4Radioactivation();

  G4bool ARMflag = false;
  radioactiveDecay->SetARM(ARMflag);        //Atomic Rearangement

  // need to initialize atomic deexcitation
  //
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* deex = man->AtomDeexcitation();
  if (!deex) {
     ///G4EmParameters::Instance()->SetFluo(true);
     G4EmParameters::Instance()->SetAugerCascade(ARMflag);
     G4EmParameters::Instance()->SetDeexcitationIgnoreCut(ARMflag);    
     deex = new G4UAtomicDeexcitation();
     deex->InitialiseAtomicDeexcitation();
     man->SetAtomDeexcitation(deex);
  }

  // register radioactiveDecay
  //
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
  
  //printout
  //
  G4cout << "\n  Set atomic relaxation mode " << ARMflag << G4endl;
  
  // Modular physics lists copied from rdecay02
  
  G4int verb = 1;
  SetVerboseLevel(verb);
  
  // EM physics
  G4EmStandardPhysics* emList = new G4EmStandardPhysics();
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetAugerCascade(true);
  param->SetStepFunction(1, 0.11*CLHEP::mm);
  param->SetStepFunctionMuHad(1., 0.1*CLHEP::mm);
  emList->ConstructProcess();

//  G4EmPenelopePhysics* peList = new G4EmPenelopePhysics();
//  param->SetAugerCascade(true);
//  param->SetStepFunction(1., 0.1*CLHEP::mm);
//  param->SetStepFunctionMuHad(1., 0.1*CLHEP::mm);
//  peList->ConstructProcess();

  // Decay
  //G4DecayPhysics* decayList = new G4DecayPhysics();
  //decayList->ConstructProcess();

  // Hadron Elastic scattering
  G4HadronElasticPhysics* hadElasList = new G4HadronElasticPhysics(verb);
  hadElasList->ConstructProcess();
  
  // Hadron Inelastic physics
  G4HadronPhysicsFTFP_BERT* hadInelasList = new G4HadronPhysicsFTFP_BERT(verb);
//  G4HadronInelasticQBBC* hadInelasList = new G4HadronInelasticQBBC(verb);        
//  G4HadronPhysicsINCLXX* hadInelasList = new G4HadronPhysicsINCLXX(verb);
  hadInelasList->ConstructProcess();
  
   // Ion Elastic scattering
  G4IonElasticPhysics* ionElasList = new G4IonElasticPhysics(verb);
  ionElasList->ConstructProcess();
  
  // Ion Inelastic physics
  G4IonPhysics* ionInelasList = new G4IonPhysics(verb);
  ////G4IonPhysics* = new G4IonINCLXXPhysics(verb);
  ionInelasList->ConstructProcess();
    
  // Gamma-Nuclear Physics
  G4EmExtraPhysics* gnuc = new G4EmExtraPhysics(verb);
  gnuc->ElectroNuclear(true);
  gnuc->MuonNuclear(true);
  gnuc->ConstructProcess();

  // Livermore EM physics
//  G4EmLivermorePhysics* livermore = new G4EmLivermorePhysics();
//  livermore->ConstructProcess();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutsWithDefault();
  // Cuts for gamma
  SetCutValue(0.1*mm, "proton");
  SetCutValue(0.1*mm, "e-");
  SetCutValue(0.1*mm, "e+");
  SetCutValue(0.1*mm, "gamma");  
//  SetCutValue(1*mm, "alpha");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
