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
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "TROOT.h"

#include "FTF_BIC.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysic G4HadronPhysicsFTF_BIC G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "FTFP_BERT_ATL.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysic G4HadronPhysicsFTFP_BERT_ATL G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut 
#include "FTFP_BERT.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysic G4HadronPhysicsFTFP_BERT G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut        
#include "FTFP_BERT_HP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsHP G4HadronPhysicsFTFP_BERT_HP G4StoppingPhysics G4IonPhysics            
#include "FTFP_BERT_TRV.hh"// G4EmStandardPhysicsGS G4EmExtraPhysics G4DecayPhysics G4HadronHElasticPhysics G4HadronPhysicsFTFP_BERT_TRV G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "FTFP_INCLXX.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsINCLXX G4StoppingPhysics G4IonINCLXXPhysics G4NeutronTrackingCut
#include "FTFP_INCLXX_HP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsHP G4HadronPhysicsINCLXX G4StoppingPhysics G4IonINCLXXPhysics
#include "G4GenericPhysicsList.hh"                    
#include "LBE.hh"// 这个比较复杂 -_-
#include "NuBeam.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsNuBeam G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "QBBC.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsXS G4StoppingPhysics G4IonPhysics G4HadronInelasticQBBC G4NeutronTrackingCut
#include "QGS_BIC.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsQGS_BIC G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "QGSP_BERT.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsQGSP_BERT G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "QGSP_BERT_HP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsHP G4HadronPhysicsQGSP_BERT_HP G4StoppingPhysics G4IonPhysics
#include "QGSP_BIC_AllHP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsPHP G4HadronPhysicsQGSP_BIC_AllHP G4StoppingPhysics G4IonPhysicsPHP
#include "QGSP_BIC.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsQGSP_BIC G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "QGSP_BIC_HP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsHP G4HadronPhysicsQGSP_BIC_HP G4StoppingPhysics G4IonPhysics
#include "QGSP_FTFP_BERT.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsQGSP_FTFP_BERT G4StoppingPhysics G4IonPhysics G4NeutronTrackingCut
#include "QGSP_INCLXX.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysics G4HadronPhysicsINCLXX G4StoppingPhysics G4IonINCLXXPhysics G4NeutronTrackingCut
#include "QGSP_INCLXX_HP.hh"// G4EmStandardPhysics G4EmExtraPhysics G4DecayPhysics G4HadronElasticPhysicsHP G4HadronPhysicsINCLXX G4StoppingPhysics G4IonINCLXXPhysics
#include "Shielding.hh"// 这个比较复杂,分好几种情况 -_-

#include "TROOT.h"

#include "G4ParticleHPManager.hh"

// 关于图形界面与交互接口
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"



#include "Randomize.hh"//随机数这里产生
#include <ctime>// initialize random seed

using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4int seconds =  time(NULL);
  G4Random::setTheSeed(seconds);

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }




  // Construct the default run manager
  G4RunManager* mtrunManager = NULL;

  #ifdef SINGLETHREAD
    mtrunManager = new G4RunManager;
  #else
    mtrunManager = new G4MTRunManager;
    mtrunManager->SetNumberOfThreads(1);
  #endif
    // mtrunManager->SetUserInitialization(new wuWorkerInitialization);


    // Set mandatory initialization classes ，410版本的框架是这样的，ActionInitialization来管理。
    //
    // Detector construction
    mtrunManager->SetUserInitialization(new DetectorConstruction());

    //physics
    mtrunManager->SetUserInitialization(new QGSP_BIC());

    // User action initialization
    mtrunManager->SetUserInitialization(new ActionInitialization());

    // Initialize G4 kernel
    mtrunManager->Initialize();


  
  // Print   Data source of this Partile HP calculation
  // G4ParticleHPManager::GetInstance()->DumpDataSource();
<<<<<<< HEAD
  ROOT::EnableThreadSafety(); // 这一行非常重要，少了程序就崩溃
=======


  ROOT::EnableThreadSafety(); // 这一行非常重要，少了程序就崩溃

>>>>>>> a5ec991fd00c2f592556d739c4ce09eb32c03836
  

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete mtrunManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
