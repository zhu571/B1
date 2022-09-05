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
/// \file DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SDManager.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VisAttributes;

/// Detector construction class to define materials and geometry.

namespace B1
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

   

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  //  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  
  
  private:
  // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    void BuildSensitiveDetector(G4LogicalVolume* lv);//

  private: 
    G4bool checkOverlaps; // Option to switch on/off checking of volumes overlaps

    std::vector<G4VisAttributes*> fVisAttributes;//可视化界面几何体颜色设置用

    G4LogicalVolume* logicWorld;
    G4VPhysicalVolume* physWorld;

    G4LogicalVolume* logicDSSD1[10000]; 
    G4LogicalVolume* logicDSSD2[10000]; 
    G4LogicalVolume* logicDSSD3[10000];
    G4VPhysicalVolume* physDSSD1[10000];
    G4VPhysicalVolume* physDSSD2[10000];
    G4VPhysicalVolume* physDSSD3[10000];
    


 // protected:
  //  G4LogicalVolume* fScoringVolume = nullptr;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
