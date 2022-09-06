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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4TrackVector.hh"
#include "G4SteppingManager.hh"

#include "G4CsvAnalysisManager.hh"
#include "G4XmlAnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{
  analysisManager = G4RootAnalysisManager::Instance();
  // analysisManager = G4CsvAnalysisManager::Instance();
  // analysisManager = G4XmlAnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Collect data step by step

  // Information in G4Step includes:
  // Pointers to PreStep and PostStepPoint
  // Geometrical step length (step length before the correction of multiple scattering)
  // True step length (step length after the correction of multiple scattering)
  // Increment of position and time between PreStepPoint and PostStepPoint
  // Increment of momentum and energy between PreStepPoint and PostStepPoint. (Note: to get the en-
  // ergy deposited in the step, you cannot use this 'Delta energy'. You have to use 'Total energy deposit' as below.)
  // Pointer to G4Track
  // Total energy deposited during the step - this is the sum of
  // • the energy deposited by the energy loss process, and
  // • the energy lost by secondaries which have NOT been generated because each of their energies was below
  // the cut threshold
  // Energy deposited not by ionization during the step
  // Secondary tracks created during tracking for the current track.
  // • NOTE: all secondaries are included. NOT only secondaries created in the CURRENT step.
  
  // step->GetTrack();//G4Track*
  // step->GetPreStepPoint();//G4StepPoint*
  // step->GetPostStepPoint();//G4StepPoint*
  // // Before the end of the AlongStepDoIt loop,StepLength keeps the initial value which is determined by the shortest geometrical Step proposed by a physics process. After finishing the AlongStepDoIt, it will be set equal to 'StepLength' in G4Step. 
  // step->GetStepLength();//G4double
  // step->GetTotalEnergyDeposit();//G4double    total energy deposit 
  // step->GetNonIonizingEnergyDeposit();//G4double    total non-ionizing energy deposit 
  // step->GetControlFlag();//G4SteppingControl    cotrole flag for stepping
  // step->GetDeltaPosition();//G4ThreeVector
  // step->GetDeltaTime();//G4double


  // Information in G4StepPoint (PreStepPoint and PostStepPoint) includes:
  // (x, y, z, t)
  // (px, py, pz, Ek)
  // Momentum direction (unit vector)
  // Pointers to physical volumes
  // Safety
  // Beta, gamma
  // Polarization
  // Step status
  // Pointer to the physics process which defined the current step and its DoIt type
  // Pointer to the physics process which defined the previous step and its DoIt type
  // Total track length
  // Global time (time since the current event began)
  // Local time (time since the current track began)
  // Proper time
  
  // steppoint->GetPosition();//G4ThreeVector&
  // steppoint->GetGlobalTime();//G4double    Time since the event in which the track belongs is created.
  // steppoint->GetLocalTime();//G4double    Time since the track is created.
  // steppoint->GetProperTime();//G4double    Proper time of the particle.
  // steppoint->GetMomentumDirection();//G4ThreeVector&
  // steppoint->GetMomentum();//G4ThreeVector    Total momentum of the track
  // steppoint->GetTotalEnergy();//G4double    Total energy of the track
  // steppoint->GetKineticEnergy();//G4double    Kinetic Energy of the track
  // steppoint->GetVelocity();//G4double
  // steppoint->GetBeta();//G4double    Velocity of the track in unit of c(light velocity)
  // steppoint->GetGamma();//G4double    Gamma factor (1/sqrt[1-beta*beta]) of the track
  // steppoint->GetPhysicalVolume();//G4VPhysicalVolume*
  // steppoint->GetTouchable();//G4VTouchable*
  // steppoint->GetTouchableHandle();//G4TouchableHandle&
  // steppoint->GetMaterial();//G4Material*
  // steppoint->GetMaterialCutsCouple();//G4MaterialCutsCouple*
  // steppoint->GetSensitiveDetector();//G4VSensitiveDetector*
  // steppoint->GetSafety();//G4double
  // steppoint->GetPolarization();//G4ThreeVector&
  // steppoint->GetStepStatus();//G4StepStatus
  // steppoint->GetProcessDefinedStep();//G4VProcess*    If the pointer is 0, this means the Step is defined by the user defined limit in the current volume.
  // steppoint->GetMass();//G4double
  // steppoint->GetCharge();//G4double
  // steppoint->GetMagneticMoment();//G4double
  // steppoint->GetWeight();//G4double


  // G4SteppingManager 类中有很多函数
  // G4SteppingManager* steppingManager = fpSteppingManager;
  // G4TrackVector* fSecondary = steppingManager->GetfSecondary();
  // for(size_t lp1 = 0;lp1 < (*fSecondary).size(); lp1++)
  //   {
  //     std::cout<<lp1<<std::endl;
  //   }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  EDep = step->GetTotalEnergyDeposit();

  G4Track* aTrack = step->GetTrack();
  G4ParticleDefinition* theparticle = aTrack->GetDefinition();
  PName = theparticle->GetParticleName();
  EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  //TrackID = aTrack->GetTrackID();
  //ParentID = aTrack->GetParentID();
  
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();


  VolNamePre = preStepPoint->GetPhysicalVolume()->GetName();
 
 
 // G4LogicalVolume
  preStepVolume = preStepPoint->GetPhysicalVolume()->GetLogicalVolume();
  if(postStepPoint->GetPhysicalVolume())//判断是否在world外
    postStepVolume = postStepPoint->GetPhysicalVolume()->GetLogicalVolume();
    
    
  PosPre = preStepPoint->GetPosition();
  PosPost = postStepPoint->GetPosition();


  analysisManager->FillNtupleIColumn(0, EventID);
  //analysisManager->FillNtupleIColumn(1, ParentID);
  //analysisManager->FillNtupleIColumn(2, TrackID);
  analysisManager->FillNtupleSColumn(1, PName);
  analysisManager->FillNtupleDColumn(2, EDep);
  analysisManager->FillNtupleSColumn(3, VolNamePre);
  analysisManager->FillNtupleDColumn(4, PosPre.x());
  analysisManager->FillNtupleDColumn(5, PosPre.y());
  analysisManager->FillNtupleDColumn(6, PosPre.z());
  analysisManager->FillNtupleDColumn(7, PosPost.x());
  analysisManager->FillNtupleDColumn(8, PosPost.y());
  analysisManager->FillNtupleDColumn(9, PosPost.z());

  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
