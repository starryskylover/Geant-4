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
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

namespace B1
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    G4VPhysicalVolume* DetectorConstruction::Construct()
    {
        // Get nist material manager
        G4NistManager* nist = G4NistManager::Instance();

        // Envelope parameters
        //
        G4double env_sizeXY = 100 * cm, env_sizeZ = 100 * cm;
        G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

        G4bool checkOverlaps = true;

        G4double world_sizeXY = 1.2 * env_sizeXY;
        G4double world_sizeZ = 1.2 * env_sizeZ;
        G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

        auto solidWorld = new G4Box("World",                           // its name
            0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

        auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
            world_mat,                                       // its material
            "World");                                        // its name


        auto physWorld = new G4PVPlacement(nullptr,  // no rotation
            G4ThreeVector(),                           // at (0,0,0)
            logicWorld,                                // its logical volume
            "World",                                   // its name
            nullptr,                                   // its mother  volume
            false,                                     // no boolean operation
            0,                                         // copy number
            checkOverlaps);                            // overlaps checking

        //
        // Envelope
        //

          auto solidEnv = new G4Box("Envelope",                    // its name
            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

        auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
            env_mat,                                     // its material
            "Envelope");                                 // its name



        new G4PVPlacement(nullptr,  // no rotation
            G4ThreeVector(),          // at (0,0,0)
            logicEnv,                 // its logical volume
            "Envelope",               // its name
            logicWorld,               // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking

        G4Material* DetectorMat_1 = nist->FindOrBuildMaterial("G4_Ge");

        G4double DetectorRadius = 66.3 / 2 * mm;
        G4double DetectorZ = 51.3 / 2 * mm;

        G4ThreeVector pos1 = G4ThreeVector(0, 0, -DetectorZ);
        auto solidShape1 = new G4Tubs(
            "TubeDetector", 0., DetectorRadius, DetectorZ, 0. * deg, 360. * deg); 

        auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
            DetectorMat_1,                                        // its material
            "Shape1");                                         // its name

        new G4PVPlacement(nullptr,  // no rotation
            pos1,                     // at position
            logicShape1,              // its logical volume
            "Shape1",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking


      // Al-container
      //
        G4Material* MaterialAlCorpus = nist->FindOrBuildMaterial("G4_Al");


        G4double AlRadiusCorpus = (66.3 + 6) / 2 * mm;
        G4double AlCorpusZ = (51.3 + 3) / 2 * mm;

        G4ThreeVector positionAlCorpus = G4ThreeVector(0, 0 * cm, -DetectorZ);

        auto AlCorpusSolid = new G4Tubs(
            "TubeDetector", AlRadiusCorpus, AlRadiusCorpus + 1 * mm, AlCorpusZ, 0. * deg, 360. * deg);


        auto logicalAlCorpus = new G4LogicalVolume(AlCorpusSolid,  // its solid
            MaterialAlCorpus,                                 // its material
            "Shape1");                                         // its name

        new G4PVPlacement(nullptr,  // no rotation
            positionAlCorpus,
            logicalAlCorpus,              // its logical volume
            "AlCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking

        auto Kryshka1Solid = new G4Tubs(
            "TubeDetector", 0, AlRadiusCorpus + 1 * mm, 0.5 * mm, 0. * deg, 360. * deg);

        auto logicKryshka1 = new G4LogicalVolume(Kryshka1Solid,  // its solid
            MaterialAlCorpus,                                        // its material
            "Shape1");                                         // its name

        G4double ZKryshka1 = (51.3 / 2 + 3) * mm;

        G4ThreeVector posKryska1 = G4ThreeVector(0, 0 * cm, ZKryshka1);

        new G4PVPlacement(nullptr,  // no rotation
            posKryska1,                     // at position
            logicKryshka1,              // its logical volume
            "AlCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking


        auto Kryshka2Solid = new G4Tubs(
            "TubeDetector", 0, AlRadiusCorpus + 1 * mm, 0.5 * mm, 0. * deg, 360. * deg);

        auto logicKryska2 = new G4LogicalVolume(Kryshka2Solid,  // its solid
            MaterialAlCorpus,                                        // its material
            "Shape1");                                         // its name

        G4double Zkryshka2 = -DetectorZ - 1 * mm;

        G4ThreeVector posKryska2 = G4ThreeVector(0, 0 * cm, Zkryshka2);

        new G4PVPlacement(nullptr,  // no rotation
            posKryska2,                     // at position //pos1
            logicKryska2,              // its logical volume
            "AlCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking


        G4Material* MaterialCup = nist->FindOrBuildMaterial("G4_Cu");

        G4double CuprumRadius = (66.3 + 6 + 2) / 2 * mm;
        G4double CuprumZ = (51.3 + 3 + 50 + 1) / 2 * mm;

        auto CuprumSolid = new G4Tubs(
            "TubeDetector", CuprumRadius, CuprumRadius + 50 * mm, 28.15 * mm, 0. * deg, 360. * deg);
        // 1mm+51.3mm+3mm+1mm

        auto logicalCuprum = new G4LogicalVolume(CuprumSolid,  // its solid
            MaterialCup,                                        // its material
            "Shape1");                                         // its name

        G4ThreeVector posCuprum = G4ThreeVector(0, 0 * cm, -DetectorZ - 1 * mm);

        new G4PVPlacement(nullptr,  // no rotation
            posCuprum,                     // at position //pos1
            logicalCuprum,              // its logical volume
            "CuCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking

        auto CupKryshkaSolid = new G4Tubs(
            "TubeDetector", 0, CuprumRadius + 50 * mm, 25 * mm, 0. * deg, 360. * deg);

        auto CupKryshka = new G4LogicalVolume(CupKryshkaSolid,  // its solid
            MaterialCup,                                        // its material
            "Shape1");                                         // its name

        G4double ZCuprumKryshka = (51.3 / 2 + 3 + 1) * mm;

        G4ThreeVector posCuprumKryshka = G4ThreeVector(0, 0 * cm, ZCuprumKryshka);

        new G4PVPlacement(nullptr,  // no rotation
            posCuprumKryshka,                     // at position
            CupKryshka,              // its logical volume
            "CuCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking



        G4Material* MaterialPb = nist->FindOrBuildMaterial("G4_Pb");

        G4double PbRadius = (66.3 + 6 + 2 + 100) / 2 * mm;

        auto PbSolid = new G4Tubs(
            "TubeDetector", PbRadius, PbRadius + 100 * mm, 53.15 * mm, 0. * deg, 360. * deg);
        // 1mm+51.3mm+3mm+1mm+50

        auto logicalPb = new G4LogicalVolume(PbSolid,  // its solid
            MaterialPb,                                        // its material
            "Shape1");                                         // its name

        G4ThreeVector posPb = G4ThreeVector(0, 0 * cm, -DetectorZ - 1 * mm);

        new G4PVPlacement(nullptr,  // no rotation
            posPb,                     // at position //pos1
            logicalPb,              // its logical volume
            "PbCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking

        auto PbKryshkaSolid = new G4Tubs(
            "TubeDetector", 0, PbRadius + 100 * mm, 50 * mm, 0. * deg, 360. * deg);

        auto PbKryshka = new G4LogicalVolume(PbKryshkaSolid,  // its solid
            MaterialPb,                                        // its material
            "Shape1");                                         // its name

        G4double ZPbKryshka = (51.3 / 2 + 3 + 1 + 50) * mm;

        G4ThreeVector posPbKryshka = G4ThreeVector(0, 0 * cm, ZPbKryshka);

        new G4PVPlacement(nullptr,  // no rotation
            posPbKryshka,                     // at position
            PbKryshka,              // its logical volume
            "PbCorpus",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking




        //
        // Shape 2
        //

        double a;
        //Set Carbon element as "elH"
        a = 1.01 * g / mole;
        G4Element* elH = new G4Element("Hydrogen", "H", 1, a);

        //Set Carbon element as "elC"
        a = 12.01 * g / mole;
        G4Element* elC = new G4Element("Carbon", "C", 6, a);

        G4double density = 1.032 * g / cm3;

        G4Material* plastic = new G4Material("PlasticSc", density, 2);
        plastic->AddElement(elH, 10);
        plastic->AddElement(elC, 11);

        G4ThreeVector pos2 = G4ThreeVector(0, 0 * cm, -50 * cm);

        // Trapezoid shape
        G4double shape2_dxa = 100 * cm, shape2_dxb = 100 * cm;
        G4double shape2_dya = 100 * cm, shape2_dyb = 100 * cm;
        G4double shape2_dz = 5 * cm;
        auto solidShape2 = new G4Trd("Shape2",  // its name
            0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
            0.5 * shape2_dz);  // its size

        auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
            plastic,                                        // its material
            "Shape2");                                         // its name

        new G4PVPlacement(nullptr,  // no rotation
            pos2,                     // at position
            logicShape2,              // its logical volume
            "Shape2",                 // its name
            logicEnv,                 // its mother  volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // overlaps checking

        // Set Shape2 as scoring volume
        //
        fScoringVolume = logicShape1;

        //
        //always return the physical World
        //
        return physWorld;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
