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
/// \file ActionInitialization.cc
/// \brief Implementation of the B1::ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1030
#include "G4MultiRunAction.hh"
#include "G4MultiEventAction.hh"
#include "G4MultiTrackingAction.hh"
#include "G4MultiSteppingAction.hh"
#endif



namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);

#if G4VERSION_NUMBER >= 1030
  
  G4MultiRunAction* actsRun = new G4MultiRunAction;
  G4MultiEventAction* actsEvent = new G4MultiEventAction;
  //G4MultiTrackingAction* actsTrack = new G4MultiTrackingAction;
  G4MultiSteppingAction* actsStep = new G4MultiSteppingAction;
  
  //...0123456789876543210...0123456789876543210...

  actsRun->push_back(G4UserRunActionUPtr(new RunAction));
  actsEvent->push_back(G4UserEventActionUPtr(new EventAction));
  //actsTrack->push_back(G4UserTrackingActionUPtr(new TrackingAction));
  actsStep->push_back(G4UserSteppingActionUPtr(new SteppingAction));
  
  //SetUserAction(new StackingAction);


  //...0123456789876543210...0123456789876543210...
  
  SetUserAction(actsRun);
  //SetUserAction(actsTrack);
  SetUserAction(actsEvent);
  SetUserAction(actsStep);
  
#else

  G4cout<<"It need G4VERSION_NUMBER >= 1030"<<G4endl;
  
#endif 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
