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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4PhysicalConstants.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),particleGun(NULL)
{
  particleGun = new G4ParticleGun(1);///*G4int n_particle*/
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if(particleGun)
    {
      delete particleGun;
      particleGun =NULL;
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pp = 0;

  // -------------------------
  
  //自定义带电粒子
  // G4int IonZ = 1;
  // G4int IonA = 1;
  // G4double IonEstar = 0.0; //exitition energy
  // G4double IonQ = 1;
  // G4cout<<"ion:Z-A-Q-E*"<<IonZ<<" "<<IonA<<" "<<IonQ<<" "<<IonEstar<<G4endl;
  // pp = particleTable->GetIonTable()->GetIon(IonZ, IonA, IonEstar);//4.10.01版本强制 G4IonTable.hh
  // particleGun->SetParticleCharge(IonQ);

  //Geant4已经定义的粒子
  pp = particleTable->FindParticle("e-");

  // -------------------------
  
  if(pp)
    particleGun->SetParticleDefinition(pp);
  else
    G4cout<<"##Null pp in PrimaryGeneratorAction::SetParticleGun()"<<G4endl;



  //particle mass
  G4double mp = 1.6735238e-27;
  G4double m21Mg = 34.890758e-27;
  G4double m20Na = 33.222993e-27;

  G4double Ed =  4.484*MeV;
  G4double Ep1 = 2.364*MeV;
  G4double Ep2 = 2.12*MeV;

  G4double E2V_Mg = std::sqrt(2*mp/(mp*mp+mp*m21Mg));
  G4double E2V_p2 = std::sqrt(2*m20Na/(mp*mp+mp*m20Na));

  G4double v23 = std::sqrt(Ep1)*E2V_Mg;
  G4double v2c = std::sqrt(Ep2)*E2V_p2;

  
  //primary particle kinetic energy
  particleGun->SetParticleEnergy(Ep1);



  // primary particle position
  G4double x1 = 0.;
  G4double y1 = 0.;
  G4double z1 = 0.;
  particleGun->SetParticlePosition(G4ThreeVector(x1, y1, z1));


  // primary particle moving direction
  G4double theta1= acos((G4UniformRand()-0.5)*2);
  G4double phi1 = G4UniformRand()*2.0*pi;
  G4double cosPX1 = sin(theta1)*cos(phi1);
  G4double cosPY1 = sin(theta1)*sin(phi1);
  G4double cosPZ1 = cos(theta1);
  G4ThreeVector directPri1(cosPX1, cosPY1, cosPZ1);    
  particleGun->SetParticleMomentumDirection(directPri1);

  //这个调用一次设置一次粒子，一次模拟要同时发射多个不同粒子就得多次调用它
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleEnergy(Ep2);


  // primary particle position
  G4double x2 = 0.;
  G4double y2 = 0.;
  G4double z2 = 0.;
  particleGun->SetParticlePosition(G4ThreeVector(x2, y2, z2));


  // primary particle moving direction
  G4double theta2= acos((G4UniformRand()-0.5)*2);
  G4double phi2 = G4UniformRand()*2.0*pi;
  G4double cosPX2 = sin(theta2)*cos(phi2);
  G4double cosPY2 = sin(theta2)*sin(phi2);
  G4double cosPZ2 = cos(theta2);
  //G4ThreeVector directPri2(cosPX2, cosPY2, cosPZ2);


  //速度合成
 
  G4double v23x = v23*cosPX1;
  G4double v23y = v23*cosPY1;
  G4double v23z = v23*cosPZ1;

  G4double v2x = v2*cosPX2;
  G4double v2y = v2*cosPY2;
  G4double v2z = v2*cosPZ2;

  G4double v2xl = v23x +v2x;
  G4double v2yl = v23y +v2y;
  G4double v2zl = v23z +v2z;


  G4double cosPX2l = v2xl/(v2xl*v2xl+v2yl*v2yl+v2zl*v2zl);
  G4double cosPY2l = v2yl/(v2xl*v2xl+v2yl*v2yl+v2zl*v2zl);
  G4double cosPZ2l = v2zl/(v2xl*v2xl+v2yl*v2yl+v2zl*v2zl);

  G4ThreeVector directPri2(cosPX2l,cosPY2l,cosPZ2l);
  ///particle
  particleGun->SetParticleMomentumDirection(directPri2);

  //这个调用一次设置一次粒子，一次模拟要同时发射多个不同粒子就得多次调用它
  particleGun->GeneratePrimaryVertex(anEvent);


  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();// return the pointer to G4ParticleTable object. G4ParticleTable is a "singleton" and can get its pointer by this function. At the first time of calling this function, the G4ParticleTable object is instantiated 
  // G4ParticleDefinition* pp = 0;
  // // ion 有以下几个设置方法，选择一个来设置即可
  // pp = particleTable->GetIonTable()->GetIon(/*G4int Z*/, /*G4int A*/, /*G4int lvl=0*/);//Z: Atomic Number   A: Atomic Mass (nn + np +nlambda)
  // pp = particleTable->GetIonTable()->GetIon(/*G4int Z*/, /*G4int A*/, /*G4int L*/, /*G4int lvl*/);//L: Number of Lmabda  E: Excitaion energy
  // pp = particleTable->GetIonTable()->GetIon(/*G4int Z*/, /*G4int A*/, /*G4double E*/, /*G4int J=0*/);//lvl:  Isomer Level 0: ground state)
  // pp = particleTable->GetIonTable()->GetIon(/*G4int Z*/, /*G4int A*/, /*G4int L*/, /*G4double E*/, /*G4int J=0*/);//J: Total Angular momentum (in unit of 1/2) : not used
  // // G4 内部定义好的粒子 如alpha nutron proton 等等,选择一个来设置即可
  // pp = particleTable->FindParticle(/*G4int  PDGEncoding*/ );// returns a pointer to the particle (0 if not contained)
  // pp = particleTable->FindParticle(/*const G4String &particle_name*/);// returns a pointer to the particle (0 if not contained)
  // pp = particleTable->FindParticle(/*const G4ParticleDefinition *particle*/);// returns a pointer to the particle (0 if not contained)
  // pp = particleTable->FindAntiParticle(/*G4int  PDGEncoding*/ );// returns a pointer to its anti-particle (0 if not contained)
  // pp = particleTable->FindAntiParticle(/*const G4String &particle_name*/);// returns a pointer to its anti-particle (0 if not contained)
  // pp = particleTable->FindAntiParticle(/*const G4ParticleDefinition *particle*/);// returns a pointer to its anti-particle (0 if not contained)
  
  // particleGun->SetParticleDefinition(/*G4ParticleDefinition * aParticleDefinition*/);
  // particleGun->SetParticleEnergy(/*G4double aKineticEnergy*/);
  // particleGun->SetParticleMomentum(/*G4double aMomentum*/);
  // particleGun->SetParticleMomentum(/*G4ParticleMomentum aMomentum*/);
  // particleGun->SetParticleMomentumDirection(/*G4ParticleMomentum aMomentumDirection*/);
  // particleGun->SetParticleCharge(/*G4double aCharge*/);
  // particleGun->SetParticlePolarization(/*G4ThreeVector aVal*/);
  // particleGun->SetNumberOfParticles(/*G4int i*/);
  // particleGun->SetParticlePosition(/*G4ThreeVector aPosition*/);
  // particleGun->SetParticleTime(/*G4double aTime*/);

  // particleGun->GeneratePrimaryVertex(anEvent);//这个调用一次设置一次粒子，一次模拟要同时发射多个不同粒子就得多次调用它

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


