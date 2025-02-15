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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("rdecay01")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
 
  // Define histograms start values
  const G4int kMaxHisto = 15;
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"};
  const G4String title[] = 
          { "dummy",                                    //0
            "energy spectrum (%): e+ e-",               //1
            "energy spectrum (%): nu_e anti_nu_e",      //2
            "energy spectrum (%): gamma",               //3
            "energy spectrum (%): alpha",               //4
            "energy spectrum (%): ions",                //5
            "total kinetic energy per single decay (Q)",//6
            "momentum balance",                         //7
            "total time of life of decay chain",        //8
            "total visible energy in decay chain",      //9
            "energy spectrum (counts): e+ e-",          //10
            "energy spectrum (counts): nu_e anti_nu_e", //11
            "energy spectrum (counts): gamma",          //12
            "energy spectrum (counts): alpha",          //13
            "energy spectrum (counts): ions"            //14
          };

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  
  // Create directories
  analysisManager->SetNtupleDirectoryName("ntuple");
  
  // Create 1st ntuple (id = 0)
  analysisManager->CreateNtuple("AlphaEnergy", "Edep");
  analysisManager->CreateNtupleDColumn("MeV");    // column Id = 0
  analysisManager->CreateNtupleDColumn("posX");   // column Id = 1
  analysisManager->CreateNtupleDColumn("posY");   // column Id = 2
  analysisManager->CreateNtupleDColumn("posZ");   // column Id = 3
  analysisManager->CreateNtupleDColumn("time");   // column Id = 4
  analysisManager->FinishNtuple();
}

//void HistoManager::FillNtuple(G4double energyAbs)
//{
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //Fill 1st ntuple ( id = 0)
  //analysisManager->FillNtupleDColumn(0, 0, energyAbs);
  //analysisManager->AddNtupleRow(0);
//}

void HistoManager::FillNtuple(G4int ntuple_ID, G4int column_ID, G4double quantity_to_fill)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(ntuple_ID, column_ID, quantity_to_fill);
  analysisManager->AddNtupleRow(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
