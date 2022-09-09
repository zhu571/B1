#define Angle_cxx
#include "Angle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Angle::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Angle.C
//      root> Angle t
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
//by  b_branchname->GetEntry(ientry); //read only this branch

   TFile* opf = new TFile("Angle.root","recreate");
   TTree *tree=new TTree("tree","dssd");


   Int_t pid = -100;
   Int_t xsidid = -100;
   Int_t ysidid = -100;
   Int_t ID = -100;

//   Int_t pxsidid=-100;
//   Int_t pysidid=-100;
//   Int_t pEventID=-100;


   tree->Branch("ID",&ID,"ID/I");
   tree->Branch("xsidid",&xsidid,"xsidid/I");
   tree->Branch("ysidid",&ysidid,"ysidid/I");
   tree->Branch("pid",&pid,"pid/I");



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   if(jentry%1000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
////////////////////////////usercode/////////////////////////////////////////////////////////////////////////

   if (ID!=EventEID)
   {
      pid = 1;
   }
   else if(ID==EventEID)
   {
      pid = 2;
   }

   ID = EventEID;
   ysidid = yid;
   xsidid = xid;

   tree->Fill();



   }



   tree->Write();
   opf->Close();
}
