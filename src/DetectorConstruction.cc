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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"//继承G4VPhysicalVolume
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4PSEnergyDeposit3D.hh"
#include "G4PSNofStep3D.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSFlatSurfaceCurrent3D.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"
#include "G4SDNeutralFilter.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4CutTubs.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),checkOverlaps(true),fVisAttributes(),logicWorld(0)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  //释放可视化界面的颜色设置变量
  for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  //快速添加SD
  // BuildSensitiveDetector(logicWorld);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //在这里先定义所有可能用到的材料
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  //调用G4自身定义好的材料
  nist->FindOrBuildMaterial("G4_Galactic");//
  //nist->FindOrBuildElement("G4_Si");
  //nist->FindOrBuildMaterial("G4_Ge");//
  // nist->FindOrBuildMaterial("G4_Pu");
  // nist->FindOrBuildMaterial("G4_H");
  // nist->FindOrBuildMaterial("G4_Al");
  // nist->FindOrBuildMaterial("G4_lH2");//
  // nist->FindOrBuildMaterial("G4_lN2");//
  // nist->FindOrBuildMaterial("G4_lO2");//
  // nist->FindOrBuildMaterial("G4_lAr");//Liquid argon
  // nist->FindOrBuildMaterial("G4_Be");
  // nist->FindOrBuildMaterial("G4_WATER");
  // nist->FindOrBuildMaterial("G4_WATER_VAPOR");//水蒸气
  // nist->FindOrBuildMaterial("G4_POLYETHYLENE");//聚乙烯
  // nist->FindOrBuildMaterial("G4_BGO");//
  // nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");//二氧化碳
  // nist->FindOrBuildMaterial("G4_LEAD_OXIDE");//氧化铅
  // nist->FindOrBuildMaterial("G4_MYLAR");//mylar膜
  // nist->FindOrBuildMaterial("G4_PLEXIGLASS");//有机玻璃
  // nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");//不锈钢
  // nist->FindOrBuildMaterial("G4_LUCITE");//透明合成树脂(有机玻璃)
  // nist->FindOrBuildMaterial("G4_CONCRETE");//混凝土
  // nist->FindOrBuildMaterial("G4_GRAPHITE");//石墨
  // nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");//二氧化硅
  // nist->FindOrBuildMaterial("G4_RUBBER_NATURAL");//天然橡胶
  // nist->FindOrBuildMaterial("G4_PbWO4");//
  // nist->FindOrBuildMaterial("G4_URANIUM_OXIDE");//氧化铀 
  // nist->FindOrBuildMaterial("G4_URANIUM_MONOCARBIDE");//碳化铀
  // nist->FindOrBuildMaterial("G4_URANIUM_DICARBIDE");//二碳化铀

  
  

  // // elements
  // G4Element* H = nist->FindOrBuildElement("H",false);//1
  // G4Element* Li = nist->FindOrBuildElement("Li",false);//3
  // G4Element* C = nist->FindOrBuildElement("C",false);//6
  // G4Element* N = nist->FindOrBuildElement("N",false);//7
  // G4Element* O = nist->FindOrBuildElement("O",false);//8
  // G4Element* Mg = nist->FindOrBuildElement("Mg",false);//12
  // G4Element* Al = nist->FindOrBuildElement("Al",false);//13
   G4Element* Si = nist->FindOrBuildElement("Si",false);//14
  // G4Element* P = nist->FindOrBuildElement("P",false);//15
  // G4Element* S =  nist->FindOrBuildElement("S",false);//16
  // G4Element* Cr = nist->FindOrBuildElement("Cr",false);//24
  // G4Element* Mn = nist->FindOrBuildElement("Mn",false);//25
  // G4Element* Fe = nist->FindOrBuildElement("Fe",false);//26
  // G4Element* Ni = nist->FindOrBuildElement("Ni",false);//28
  // G4Element* I = nist->FindOrBuildElement("I",false);//53
  // G4Element* Cs = nist->FindOrBuildElement("Cs",false);//55
  // G4Element* Ce = nist->FindOrBuildElement("Ce",false);//58
  
  //elements material
  //G4double fractionmass;
  G4Material* dssdmat = new G4Material("dssdmat",2.32*g/cm3,1);
  dssdmat->AddElement(Si, 1.00);


  // G4Isotope* U4 = new G4Isotope("U234",92,234,234.02*g/mole);
  // G4Isotope* U5 = new G4Isotope("U235",92,235,235.01*g/mole);
  // G4Isotope* U6 = new G4Isotope("U236",92,236,236.04*g/mole);
  // G4Isotope* U8 = new G4Isotope("U238",92,238,238.03*g/mole);
  // G4Element* HEU58 = new G4Element("Highly-enriched Uranium 58", "HEU58", 2);
  // HEU58->AddIsotope(U5, 0.93);
  // HEU58->AddIsotope(U8, 0.07);

  // G4Element* HEU4568 = new G4Element("Highly-enriched Uranium 4568","HEU4568",4);
  // HEU4568->AddIsotope(U4,0.0097);
  // HEU4568->AddIsotope(U5,0.9315);
  // HEU4568->AddIsotope(U6,0.0024);
  // HEU4568->AddIsotope(U8,0.0564);
  
  //---------------------------------------------------------------------------------

  // //Scintillator(BC408) 塑闪 
  // G4Material* BC408 = new G4Material("BC408", 1.032*g/cm3, 2);
  // BC408->AddElement(H, 11);BC408->AddElement(C, 10);
  // BC408->AddElement(H, 10);BC408->AddElement(C, 9);

  // // LiquidScint(NE213) 液闪 
  // G4Material* NE213 = new G4Material("NE213",0.874*g/cm3,2);
  // NE213->AddElement(H,1212);
  // NE213->AddElement(C,1000);
  
  // // He-3 detector materials
  // G4Material* matHe3  = new G4Material("He3",  2., 3.*g/mole, 0.00049*g/cm3, kStateGas);
  
  // // Uranium material
  // G4Material* matHEU58 = new G4Material("HEU58", 19.1*g/cm3, 1, kStateSolid);
  // matHEU58->AddElement(HEU58, 1.00);

  // G4Material* matHEU4568 = new G4Material("HEU4568",18.75*g/cm3,1);
  // matHEU4568->AddElement(HEU4568, 1.0);

  
  // G4Material* matSteel = new G4Material("Steel",7.788*g/cm3,9);
  // matSteel->AddElement(Fe,62.1805*perCent);
  // matSteel->AddElement(Cr,20.26*perCent);
  // matSteel->AddElement(Mn,9.37*perCent);
  // matSteel->AddElement(Ni,7.5*perCent);
  // matSteel->AddElement(Si,0.34*perCent);
  // matSteel->AddElement(N,0.29*perCent);
  // matSteel->AddElement(C,0.04*perCent);
  // matSteel->AddElement(P,0.018*perCent);
  // matSteel->AddElement(S,0.0015*perCent);
  
  
  // Print materials，运行时在终端输出材料信息
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  //通过G4Material::GetMaterial()获取DefineMaterials()中定义的材料！
  G4Material* world_mat =  G4Material::GetMaterial("G4_Galactic");
  G4Material* dssd_mat =  G4Material::GetMaterial("dssdmat");

  //     
  // World
  //
  G4double sizeXYZ = 1.0*m;

  /*G4double dssdX = 5.0*cm;
  G4double dssdY = 5.0*cm;*/
  G4double dssd1Z = 142.0*um;
  G4double dssd2Z = 40.0*um;
  G4double dssd3Z = 304.0*um;
  G4double dssddx = 3.125*mm;
  G4double dssddy = 3.125*mm;

  G4Box* solidWord = 
    new G4Box("World",
    0.5*sizeXYZ,0.5*sizeXYZ,0.5*sizeXYZ);
  
  logicWorld =
    new G4LogicalVolume(solidWord,
                        world_mat,
                        "World");
  
  physWorld = 
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicWorld,
                      "World",
                      0,
                      false,
                      0,
                      checkOverlaps);

  char dssd1name[200];
  char dssd2name[200];
  char dssd3name[200];
/*
  G4Box* solidDSSD1 =    
    new G4Box("DSSD1",                    //its name
  	      0.5*dssddx, 0.5*dssddy, 0.5*dssd1Z); //its size

*/
  /*
  G4Box* solidDSSD2 =    
    new G4Box("DSSD2",                    //its name
  	      0.5*dssddx, 0.5*dssddy, 0.5*dssd2Z); //its size
*/

  G4Box* solidDSSD3 =    
    new G4Box("DSSD3",                    //its name
  	      0.5*dssddx, 0.5*dssddy, 0.5*dssd3Z); //its size

 /* logicDSSD =
    new G4LogicalVolume(solidDSSD,            //its solid
			dssd_mat,             //its material
			"DSSD");         //its name

  physDSSD =
    new G4PVPlacement(0,                       //no rotation set 0
		      G4ThreeVector(), //at (0,0,0)
		      logicDSSD,               //its logical volume
		      "DSSD",              //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      checkOverlaps);          //overlaps checking
*/

    /*
  for (size_t i = 0; i < 16; i++)
{
  for (size_t j = 0; j < 16; j++)
  {
    int num1 = 16*i+j;
    sprintf(dssd2name,"dssd2|%d",num1);

    logicDSSD2[num1] =
     new G4LogicalVolume(solidDSSD2,
                         dssd_mat,
                         dssd2name);
                      
    physDSSD2[num1] = 
    new G4PVPlacement(0,
                      G4ThreeVector((24.84375-3.125*j)*mm,(24.84375-3.125*i)*mm,0),
                      logicDSSD2[num1],
                      dssd2name,
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
  }
}
  for (size_t l = 0; l < 16; l++)
{
  for (size_t m = 0; m < 16; m++)
  {
    int num2 = 16*l+m;
    sprintf(dssd1name,"dssd1|%d",num2);

    logicDSSD1[num2] =
     new G4LogicalVolume(solidDSSD1,
                         dssd_mat,
                         dssd1name);
                      
    physDSSD1[num2] = 
    new G4PVPlacement(0,
                      G4ThreeVector((24.84375-3.125*m)*mm,(24.84375-3.125*l)*mm,19*mm),
                      logicDSSD1[num2],
                      dssd1name,
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
  }
}*/
  for (size_t u = 0; u < 16; u++)
{
  for (size_t o = 0; o < 16; o++)
  {
    int num3 = 16*u+o;
    sprintf(dssd3name,"dssd3|%d",num3);

    logicDSSD3[num3] =
     new G4LogicalVolume(solidDSSD3,
                         dssd_mat,
                         dssd3name);
                      
    physDSSD3[num3] = 
    new G4PVPlacement(0,
                      G4ThreeVector((24.84375-3.125*o)*mm,(24.84375-3.125*u)*mm,152*um),
                      logicDSSD3[num3],
                      dssd3name,
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
  }
  
}



  //===============================


    // visualization attributes ------------------------------------------------
    //可视化界面几何体颜色设置（可有可无）
  /* 
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    //visAttributes->SetVisibility(false);//不显示边框s
    logicWorld->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,0.0)); 
    for (size_t i = 0; i < 16; i++)
    {
      for (size_t j = 0; j < 16; j++)
      {
        int numc = 10*i+j;

        logicDSSD1[numc]->SetVisAttributes(visAttributes);
        fVisAttributes.push_back(visAttributes);

      }
      
    }

    for (size_t i = 0; i < 16; i++)
    {
      for (size_t j = 0; j < 16; j++)
      {
        int numc = 10*i+j;

        logicDSSD2[numc]->SetVisAttributes(visAttributes);
        fVisAttributes.push_back(visAttributes);

      }
      
    }

    for (size_t i = 0; i < 16; i++)
    {
      for (size_t j = 0; j < 16; j++)
      {
        int numc = 10*i+j;

        logicDSSD3[numc]->SetVisAttributes(visAttributes);
        fVisAttributes.push_back(visAttributes);

      }
      
    }
  
   */
//    logicDSSD->SetVisAttributes(visAttributes);
//    fVisAttributes.push_back(visAttributes);

    // G4Colour  white   (1.0, 1.0, 1.0) ;
    // G4Colour  grey    (0.5, 0.5, 0.5) ;
    // G4Colour  lgrey   (.85, .85, .85) ;
    // G4Colour  red     (1.0, 0.0, 0.0) ;
    // G4Colour  blue    (0.0, 0.0, 1.0) ;
    // G4Colour  cyan    (0.0, 1.0, 1.0) ;
    // G4Colour  magenta (1.0, 0.0, 1.0) ; 
    // G4Colour  yellow  (1.0, 1.0, 0.0) ;
    // G4Colour  orange  (.75, .55, 0.0) ;
    // G4Colour  lblue   (0.0, 0.0, .75) ;
    // G4Colour  lgreen  (0.0, .75, 0.0) ;
    // G4Colour  green   (0.0, 1.0, 0.0) ;
    // G4Colour  brown   (0.7, 0.4, 0.1) ;

  //===================================================================
  //
  //always return the physical World
  //
  return physWorld;  

}

void DetectorConstruction::BuildSensitiveDetector(G4LogicalVolume* lv)
{
  // if(!lv){
  //   G4cout<<"$$No given logical volume, no SD built"<<G4endl;
  //   return;
  // }
  // G4SDManager* SDman=G4SDManager::GetSDMpointer();
  // wuSensitiveDetector* sd1=new wuSensitiveDetector(lv->GetName()); //
  // SDman->AddNewDetector(sd1);
  // lv->SetSensitiveDetector(sd1);
}

}
