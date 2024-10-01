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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "DetectorConstruction.hh"
#include "TrackingMessenger.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "Analysis.hh"

#include "G4PSEnergyDeposit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fWorldSize = .3*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const int n_crystals = 1;

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //
  // Define the dimensions of the NaI crystal and its surrounding parts.
  //
  const auto body_length = 254*mm;
  const auto body_side = 109*mm;
  const auto body_thickness = 1*mm;

  const auto body_inner_length = body_length - body_thickness;
  const auto body_inner_side = body_side - body_thickness - 5*mm;

  const auto crystal_side = 102*mm;
  const auto crystal_length = 203*mm;

  const auto support_side = crystal_side;
  const auto support_length = 1*mm;

  const auto base_side = body_inner_side;
  const auto base_length = 3.4*mm;

  const auto SiPM_side = 25*mm;
  const auto SiPM_length = 1.35*mm;

  const auto quartz_side = crystal_side;
  const auto quartz_length = 12.5*mm;
  
  //
  // Define the materials used in the NaI scintillator detector. 
  //   
  G4NistManager *nist = G4NistManager::Instance();
  
  G4Material* matAl = nist->FindOrBuildMaterial("G4_Al");
  G4Material* matNaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material* matQuartz = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* matMylar = nist->FindOrBuildMaterial("G4_MYLAR");
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* matGalactic = nist->FindOrBuildMaterial("G4_Galactic");
  
  //     
  // Define the 'World' volume. 
  //
  G4Box* solidWorld = new G4Box("World", fWorldSize/2, fWorldSize/2, fWorldSize/2);                
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, matAir, "logicWorld");
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0); 
  
  //
  // Define the NaI solid volumes. 
  //  
  //G4Box *solidBody_outer_total = new G4Box("Body_outer_total", body_side*0.5, body_side*0.5, body_length*0.5);
  //G4Box *solidBody_inner = new G4Box("Body_inner", body_inner_side*0.5, body_inner_side*0.5, body_inner_length*0.5);
  //G4Box *solidBody_outer = new G4Box("Body_outer", solidBody_outer_total, solidBody_inner);
  G4Box *solidCrystal = new G4Box("solidCrystal", crystal_side*0.5, crystal_side*0.5, crystal_length*0.5);
  G4Box *solidQuartz = new G4Box("solidQuartz", quartz_side*0.5, quartz_side*0.5, quartz_length*0.5);
  G4Box *solidSiPM = new G4Box("soliddSiPM", SiPM_side*0.5, SiPM_side*0.5, SiPM_length*0.5);
  G4Box *solidBase = new G4Box("solidBase", base_side*0.5, base_side*0.5, base_length*0.5);
  G4Box *solidSupport = new G4Box("solidSupport", support_side*0.5, support_side*0.5, support_length*0.5);
  
  //
  // Define the NaI logical volumes. 
  //
  G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal, matNaI, "logicCrystal");
  G4LogicalVolume* logicQuartz = new G4LogicalVolume(solidQuartz, matQuartz, "logicQuartz");
  G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, matAl, "logicSiPM");
  G4LogicalVolume* logicBase = new G4LogicalVolume(solidBase, matMylar, "logicBase");
  G4LogicalVolume* logicSupport = new G4LogicalVolume(solidSupport, matAl, "logicSupport");
  
  	//
  	// Define visualization attributes of NaI detector components. 
	//
	//logicBody_outer->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
	//logicBody_inner->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
	logicCrystal->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));
	logicQuartz->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
	logicSiPM->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
	logicSupport->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
	logicBase->SetVisAttributes(G4VisAttributes(G4Colour::Magenta()));
  
  //
  // Define the NaI physical volumes. 
  //
  
  	//
  	// Before placing, define some spatial constants. 
  	//
  	const auto sep = 1*um;
      	const auto body_z = -body_length*0.5;
      	const auto base_z = -1*(-body_inner_length*0.5 + base_length*0.5 + sep);
      	const auto support_z = -1*(-body_inner_length*0.5 + base_length + support_length*0.5 + sep);
      	const auto crystal_z = -1*(-body_inner_length*0.5 + base_length + support_length + crystal_length*0.5 + sep);
      	const auto quartz_z = -1*(-body_inner_length*0.5 + base_length + support_length + crystal_length + quartz_length*0.5 + sep);
      	const auto SiPM_z = -1*(-body_inner_length*0.5 + base_length + support_length + crystal_length + quartz_length + SiPM_length*0.5 					+ sep);
  
  G4VPhysicalVolume* physiCrystal = new G4PVPlacement(0, G4ThreeVector(), logicCrystal, "physiWorld", logicWorld, false, 0, true); 
  G4VPhysicalVolume* physiQuartz = new G4PVPlacement(0, G4ThreeVector(0., 0., -(crystal_length*0.5 + quartz_length*0.5)), logicQuartz, 						"physiWorld", logicWorld, false, 0, true); 
  G4VPhysicalVolume* physiSiPM = new G4PVPlacement(0, G4ThreeVector(0., 0., -(crystal_length*0.5 + quartz_length + SiPM_length*0.5)), 					logicSiPM, "physiSiPM", logicWorld, false, 0, true); 
  G4VPhysicalVolume* physiBase = new G4PVPlacement(0, G4ThreeVector(0., 0., crystal_length*0.5 + support_length + base_length*0.5), 						logicBase, "physiWorld", logicWorld, false, 0, true); 
  G4VPhysicalVolume* physiSupport = new G4PVPlacement(0, G4ThreeVector(0., 0., crystal_length*0.5 + support_length*0.5), logicSupport, 						"physiWorld", logicWorld, false, 0, true); 
                
  //
  // Always return the physical World.
  //  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "rdecay01/TrackerChamberSD";
  TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name of "Chamber_LV".
  SetSensitiveDetector("logicCrystal", aTrackerSD, true);
}

