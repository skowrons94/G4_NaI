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
#include "G4SubtractionSolid.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "Analysis.hh"

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fWorldSize = .3*m;
  fShieldThickness = 15*cm;  // Default shield thickness
}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define the dimensions of the NaI crystal and its surrounding parts
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
  
  // Define materials
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* matPb = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* matAl = nist->FindOrBuildMaterial("G4_Al");
  G4Material* matNaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material* matQuartz = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* matMylar = nist->FindOrBuildMaterial("G4_MYLAR");
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");
  
  // Calculate shield dimensions
  const auto maxDetectorDim = std::max(crystal_side, crystal_length);
  const auto shieldSide = maxDetectorDim + 2*fShieldThickness;
  const auto shieldLength = crystal_length + quartz_length + SiPM_length + 
                           base_length + support_length + 2*fShieldThickness;

  // Create world volume
  const auto worldSize = std::max(fWorldSize, 2.2 * std::max(shieldSide, shieldLength));
  G4Box* solidWorld = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, matAir, "World");
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0, G4ThreeVector(), 
                                                   logicWorld, "World", 0, false, 0);

  // Create lead shield
  G4Box* solidShieldOuter = new G4Box("shieldOuter", 
                                     shieldSide/2, shieldSide/2, shieldLength/2);
  G4Box* solidShieldInner = new G4Box("shieldInner", 
                                     (shieldSide - 2*fShieldThickness)/2,
                                     (shieldSide - 2*fShieldThickness)/2,
                                     (shieldLength - 2*fShieldThickness)/2);
  
  G4SubtractionSolid* solidShield = new G4SubtractionSolid("shield",
                                                          solidShieldOuter,
                                                          solidShieldInner);
  
  G4LogicalVolume* logicShield = new G4LogicalVolume(solidShield, matPb, "Shield");

  // Create detector components
  G4Box* solidCrystal = new G4Box("Crystal", 
                                 crystal_side/2, crystal_side/2, crystal_length/2);
  G4Box* solidQuartz = new G4Box("Quartz", 
                                quartz_side/2, quartz_side/2, quartz_length/2);
  G4Box* solidSiPM = new G4Box("SiPM", 
                              SiPM_side/2, SiPM_side/2, SiPM_length/2);
  G4Box* solidBase = new G4Box("Base", 
                              base_side/2, base_side/2, base_length/2);
  G4Box* solidSupport = new G4Box("Support", 
                                 support_side/2, support_side/2, support_length/2);

  // Create logical volumes
  G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal, matNaI, "Crystal");
  G4LogicalVolume* logicQuartz = new G4LogicalVolume(solidQuartz, matQuartz, "Quartz");
  G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, matAl, "SiPM");
  G4LogicalVolume* logicBase = new G4LogicalVolume(solidBase, matMylar, "Base");
  G4LogicalVolume* logicSupport = new G4LogicalVolume(solidSupport, matAl, "Support");

  // Set visualization attributes
  logicShield->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));
  logicCrystal->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));
  logicQuartz->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  logicSiPM->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
  logicSupport->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
  logicBase->SetVisAttributes(G4VisAttributes(G4Colour::Magenta()));

  // Place the shield
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), 
                    logicShield, "Shield", logicWorld, false, 0, true);

  // Define component positions relative to shield center
  const auto sep = 1*um;  // Small separation to avoid overlap
  const auto crystal_z = 0.0;  // Center crystal in shield
  const auto quartz_z = -(crystal_length/2 + quartz_length/2 + sep);
  const auto SiPM_z = -(crystal_length/2 + quartz_length + SiPM_length/2 + sep);
  const auto base_z = (crystal_length/2 + base_length/2 + sep);
  const auto support_z = (crystal_length/2 + base_length + support_length/2 + sep);

  // Place detector components
  new G4PVPlacement(0, G4ThreeVector(0, 0, crystal_z), 
                    logicCrystal, "Crystal", logicWorld, false, 0, true);
  
  new G4PVPlacement(0, G4ThreeVector(0, 0, quartz_z), 
                    logicQuartz, "Quartz", logicWorld, false, 0, true);
  
  new G4PVPlacement(0, G4ThreeVector(0, 0, SiPM_z), 
                    logicSiPM, "SiPM", logicWorld, false, 0, true);
  
  new G4PVPlacement(0, G4ThreeVector(0, 0, base_z), 
                    logicBase, "Base", logicWorld, false, 0, true);
  
  new G4PVPlacement(0, G4ThreeVector(0, 0, support_z), 
                    logicSupport, "Support", logicWorld, false, 0, true);

  return physiWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  G4String trackerChamberSDname = "rdecay01/TrackerChamberSD";
  TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector("Crystal", aTrackerSD, true);
}