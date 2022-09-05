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
/// \file RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4CsvAnalysisManager.hh"
#include "G4XmlAnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RunAction::BeginOfRunAction(const G4Run*)
{
  analysisManager = G4RootAnalysisManager::Instance();
  // analysisManager = G4CsvAnalysisManager::Instance();
  // analysisManager = G4XmlAnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  analysisManager->CreateNtuple("t", "Geant4 data !");

  analysisManager->CreateNtupleIColumn("EventID");
  //analysisManager->CreateNtupleIColumn("ParentID");
  //analysisManager->CreateNtupleIColumn("TrackID");
  //analysisManager->CreateNtupleIColumn("CurrentStepNumber");
  analysisManager->CreateNtupleSColumn("PName");
  //analysisManager->CreateNtupleDColumn("TrackWeight");
  //analysisManager->CreateNtupleSColumn("CreatorProcess");
  analysisManager->CreateNtupleDColumn("EDep");
  //analysisManager->CreateNtupleDColumn("TrackLength");
  //analysisManager->CreateNtupleDColumn("StepLength");
  //analysisManager->CreateNtupleIColumn("TrackStatus");
  //analysisManager->CreateNtupleDColumn("Mass");
  //analysisManager->CreateNtupleDColumn("Charge");
  //analysisManager->CreateNtupleDColumn("MagneticMoment");
  analysisManager->CreateNtupleSColumn("VolNamePre");
  //analysisManager->CreateNtupleSColumn("VolNamePost");
  //analysisManager->CreateNtupleDColumn("GlobalTimePre");
  //analysisManager->CreateNtupleDColumn("GlobalTimePost");
  //analysisManager->CreateNtupleDColumn("LocalTimePre");
  //analysisManager->CreateNtupleDColumn("LocalTimePost");
  //analysisManager->CreateNtupleDColumn("ProperTimePre");
  //analysisManager->CreateNtupleDColumn("ProperTimePost");
  //analysisManager->CreateNtupleIColumn("StepStatusPre");
  //analysisManager->CreateNtupleIColumn("StepStatusPost");
  //analysisManager->CreateNtupleDColumn("EkPre");
  //analysisManager->CreateNtupleDColumn("EkPost");
  //analysisManager->CreateNtupleDColumn("xPre");
  //analysisManager->CreateNtupleDColumn("yPre");
  //analysisManager->CreateNtupleDColumn("zPre");
  //analysisManager->CreateNtupleDColumn("xPost");
  //analysisManager->CreateNtupleDColumn("yPost");
  //analysisManager->CreateNtupleDColumn("zPost");
  //analysisManager->CreateNtupleDColumn("xMomentumDirectionPre");
  //analysisManager->CreateNtupleDColumn("yMomentumDirectionPre");
  //analysisManager->CreateNtupleDColumn("zMomentumDirectionPre");
  //analysisManager->CreateNtupleDColumn("Velocitypre");
  //analysisManager->CreateNtupleDColumn("xMomentumDirectionPost");
  //analysisManager->CreateNtupleDColumn("yMomentumDirectionPost");
  //analysisManager->CreateNtupleDColumn("zMomentumDirectionPost");
  //analysisManager->CreateNtupleDColumn("Velocitypost");

  
  analysisManager->FinishNtuple();
  analysisManager->OpenFile("outputdata");//输出文件名
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
