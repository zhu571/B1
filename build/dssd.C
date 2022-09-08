#define dssd_cxx
#include "dssd.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void dssd::Loop()
{



//   In a ROOT session, you can do:
//      root> .L dssd.C
//      root> dssd t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branchopf

   TFile* opf = new TFile("dssd.root","recreate");
   TTree *tree=new TTree("tree","dssd");

   Int_t ID = -1;
   TH1F *h1 = new TH1F("h1","",1000,0,10);
   Double_t ene =0;
   Int_t dssdxid =0;
   Int_t dssdyid =0;



   //tree->Branch("ene",  &ene,  "ene/D");
   tree->Branch("dssdxid",  &dssdxid,  "dssdxid/I");
   tree->Branch("dssdyid",  &dssdyid,  "dssdyid/I");





   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //user code
   /*
      if (VolNamePre[0]!='D') continue;
      if (PName[0]!='e') continue;

      if (EventID!=ID)
      {
         if (ID!=-1)
         {
            if (ene>0)
            {
               h1->Fill(ene);
               tree->Fill();
            }
         
         }
            ID = EventID;
            ene = 0;
         
      }
      ene += EDep;
      
      
   */
  if (xPre >0)
  {
      dssdxid = xPre/0.3125 + 1;
  }
   else if (xPre < 0)
   {
      dssdxid = xPre/0.3125 - 1;
   }
   else if (xPre = 0)
   {
      dssdxid = 0;
   }

   
    if (yPre >0)
  {
      dssdxid = yPre/0.3125 + 1;
  }
   else if (yPre < 0)
   {
      dssdxid = yPre/0.3125 - 1;
   }
   else if (y2Pre = 0)
   {
      dssdxid = 0;
   }

  
  
   



   }
   TCanvas* c1 = new TCanvas("c1");
   h1->Draw();
   h1->Write();
   tree->Write();
   opf->Close();
}
