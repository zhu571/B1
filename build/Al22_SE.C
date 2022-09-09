#define Al22_SE_cxx
#include "Al22_SE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Al22_SE::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Al22_SE.C
//      root> Al22_SE t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current TreeC
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
//by  b_branchname->GetEntry(ientry); //read only this branch


   TFile* opf = new TFile("Al22_SE02.root","recreate");
   TTree *tree=new TTree("tree","dssd");


   Int_t ID = -1;
//   Double_t ene =0;
   Int_t dssdxid =0;
   Int_t dssdyid =0;
   Int_t dssdxidpre = 0;
   Int_t dssdyidpre = 0;
   Double_t dssdzpos =0;
   Int_t pid = 0;
   Int_t Pid = 1;
//   Char_t VolNameAPre[10];
//   Char_t VolNameAPost[1000];

//   VolNameAPre[0] = 'g';

   //tree->Branch("ene",  &ene,  "ene/D");
   tree->Branch("ID", &ID, "ID/I");
   tree->Branch("dssdxid",  &dssdxid,  "dssdxid/I");
   tree->Branch("dssdyid",  &dssdyid,  "dssdyid/I");
   tree->Branch("dssdzpos",  &dssdzpos,  "dssdzpos/D");
   tree->Branch("pid",&pid,"pid/I");
//   tree->Branch("VolNameAPre", &VolNameAPre, "VolNameAPre/C");
//   tree->Branch("VolNameAPost", &VolNameAPost, "VolNameAPost/C");




   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
//////user code///////////////////////////////////////////////////////////////////////////////////////////////////////////


  if (xPre >0)
   {
      dssdxid = xPre/3.125 + 1;
   }
   else if (xPre < 0)
   {
      dssdxid = xPre/3.125 - 1;
   }
   else if (xPre == 0)
   {
      dssdxid = 0;
   }

   
   if (yPre >0)
   {
      dssdyid = yPre/3.125 + 1;
   }
   else if (yPre < 0)
   {
      dssdyid = yPre/3.125 - 1;
   }     
   else if (yPre == 0)
   {
      dssdyid = 0;
   }



   dssdzpos = zPre;


   //   *VolNameAPost = *VolNamePost;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//if (VolNamePre[0]!='d') continue;

   if (abs(zPre)>15)
   {

      if (dssdxid!=dssdxidpre||dssdyid!=dssdyidpre)
      {

      if (ID==EventID)
      {
         pid = 2;
      }
      else if(ID!=EventID)
      {
         pid = 1;
      }

      if (pid==2||Pid==2)
      {
         tree->Fill();

      }
  
      
      Pid =pid;

      ID=EventID;
   }
   dssdyidpre = dssdyid;
   dssdxidpre = dssdxid;

   }







  
   }
   tree->Write();
   opf->Close();

}
