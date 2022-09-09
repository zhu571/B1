#define Al22_Angle_cxx
#include "Al22_Angle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Al22_Angle::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Al22_Angle.C
//      root> Al22_Angle t
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


   TFile* opf = new TFile("Al22_Angle.root","recreate");
   TTree *tree=new TTree("tree","dssd");

   Int_t xid =-100;
   Int_t yid =-100;

   Int_t EventEID = 0;


   tree->Branch("EventEID",&EventEID,"EventEID/I");
   tree->Branch("xid",&xid,"xid/I");
   tree->Branch("yid",&yid,"yid/I");







   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   if(jentry%1000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
/////////////////user code////////////////////////////////////////////////////////////////////////////////////

      if (xid == dssdxid && yid == dssdyid) continue;

      xid = dssdxid;
      yid = dssdyid;
      EventEID = EventAID;

      tree->Fill();

   }
   tree->Write();
   opf->Close();
}
