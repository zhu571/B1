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
/// \file SteppingAction.hh
/// \brief Definition of the B1::SteppingAction class

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

class G4LogicalVolume;
class G4VAnalysisManager;

/// Stepping action class
///

namespace B1
{

class EventAction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    virtual ~SteppingAction();

    // method from the base class
  public:
    virtual void UserSteppingAction(const G4Step*);

  private:
  
  G4int    EventID;//
  G4int    ParentID;//
  G4int    TrackID;//
  G4int    CurrentStepNumber;//       // Total steps number up to now
  G4String PName;//                           //particle name
  G4double TrackWeight;//                        // Track Weight
  G4String CreatorProcess;//
  G4double EDep;//
  G4double TrackLength;//          // Accumulated track length
  G4double StepLength;//
  G4int    TrackStatus;//

  G4double Mass;//                 // Dynamical mass of the particle静质量
  G4double Charge;//             // Dynamical Charge of the particle
  G4double MagneticMoment;//    // Dynamical MagneticMoment of the particle


  G4String VolNamePre;//                   //G4VPhysicalVolume name
  G4String VolNamePost;//
  G4double GlobalTimePre;//                   //Time since event is created
  G4double GlobalTimePost;//
  G4double LocalTimePre;//                  // Time since track is created
  G4double LocalTimePost;//
  G4double ProperTimePre;//            // Time since track is created (in rest frame of particle)
  G4double ProperTimePost;//
  G4int    StepStatusPre;//
  G4int    StepStatusPost;//
  G4double EkPre;//                    
  G4double EkPost;//
  G4ThreeVector PosPre;//
  G4ThreeVector PosPost;//
  G4ThreeVector MomentumDirectionPre;//
  G4double VelocityPre;//
  G4ThreeVector MomentumDirectionPost;//
  G4double VelocityPost;//
  G4int Info;//

  G4LogicalVolume * preStepVolume;
  G4LogicalVolume * postStepVolume;
  
  G4VAnalysisManager* analysisManager;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
