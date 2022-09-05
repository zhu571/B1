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
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4EventManager.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

#include "G4DigiManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    // anEvent->GetEventID();//  Returns the event ID
  
  // // This method sets a new primary vertex. This method must be invoked exclusively by G4VPrimaryGenerator concrete class.
  // anEvent->AddPrimaryVertex(G4PrimaryVertex* aPrimaryVertex);

  // //  Returns i-th primary vertex of the event.
  // G4PrimaryVertex *primaryvertex = anEvent->GetPrimaryVertex(/*G4int i=0*/);// return G4PrimaryVertex*
  // primaryvertex->GetPosition();//G4ThreeVector
  // primaryvertex->GetX0();//G4double
  // primaryvertex->GetY0();//G4double
  // primaryvertex->GetZ0();//G4double
  // primaryvertex->GetT0();//G4double
  // primaryvertex->GetNumberOfParticle();//G4int
  // primaryvertex->GetPrimary(/*G4int i=0*/);//G4PrimaryParticle*
  // primaryvertex->GetWeight();//G4double
  // primaryvertex->GetUserInformation();//G4VUserPrimaryVertexInformation*
  // primaryvertex->SetUserInformation(/*G4VUserPrimaryVertexInformation* anInfo*/);

  // G4PrimaryParticle* primaryparticle = primaryvertex->GetPrimary(/*G4int i=0*/);
  // primaryparticle->GetPDGcode();//G4int
  // primaryparticle->GetG4code();//G4ParticleDefinition*
  // primaryparticle->GetMass();//G4double
  // primaryparticle->GetCharge();//G4double
  // primaryparticle->GetKineticEnergy();//G4double
  // primaryparticle->GetMomentumDirection();//G4ThreeVector&
  // primaryparticle->GetTotalMomentum();//G4double
  // primaryparticle->GetTotalEnergy();//G4double
  // primaryparticle->GetMomentum();//G4ThreeVector
  // primaryparticle->GetPx();//G4ThreeVector
  // primaryparticle->GetPy();//G4ThreeVector
  // primaryparticle->GetPz();//G4ThreeVector
  // primaryparticle->GetTrackID();//G4int     "trackID" will be set if this particle is sent to G4EventManager and converted to G4Track. Otherwise = -1.
  // primaryparticle->GetPolarization();//G4ThreeVector
  // primaryparticle->GetPolX();//G4double
  // primaryparticle->GetPolY();//G4double
  // primaryparticle->GetPolZ();//G4double
  // primaryparticle->GetWeight();//G4double
  // primaryparticle->GetProperTime();//G4double   
  // primaryparticle->GetUserInformation();//G4VUserPrimaryParticleInformation*
  // primaryparticle->SetUserInformation(/*G4VUserPrimaryParticleInformation* anInfo*/);

  
  // anEvent->SetUserInformation(/*G4VUserEventInformation* anInfo*/);
  // anEvent->GetUserInformation();//return G4VUserEventInformation*

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  // // 读取该事件SD中记录的数据
  // G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();//SD  (hits collections of this event)
  // if(HCE)
  //   {
  //     size_t nHCinHCE = HCE->GetCapacity();//SD个数
  //     // 
  //     for(size_t i = 0; i < nHCinHCE; i++)//遍历SD
  // 	{
  // 	  // 拿到HitCollection

  // 	}
  //   }

  // G4DCofThisEvent* *DCE = anEvent->GetDCofThisEvent();// (digi collections of this event)

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
