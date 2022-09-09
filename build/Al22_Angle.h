//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  8 21:05:56 2022 by ROOT version 6.26/06
// from TTree tree/dssd
// found on file: Al22_SE01.root
//////////////////////////////////////////////////////////

#ifndef Al22_Angle_h
#define Al22_Angle_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Al22_Angle {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventAID;
   Int_t           dssdxid;
   Int_t           dssdyid;
   Double_t        dssdzpos;

   // List of branches
   TBranch        *b_EventAID;   //!
   TBranch        *b_dssdxid;   //!
   TBranch        *b_dssdyid;   //!
   TBranch        *b_dssdzpos;   //!

   Al22_Angle(TTree *tree=0);
   virtual ~Al22_Angle();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Al22_Angle_cxx
Al22_Angle::Al22_Angle(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Al22_SE01.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Al22_SE01.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

Al22_Angle::~Al22_Angle()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Al22_Angle::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Al22_Angle::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Al22_Angle::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventAID", &EventAID, &b_EventAID);
   fChain->SetBranchAddress("dssdxid", &dssdxid, &b_dssdxid);
   fChain->SetBranchAddress("dssdyid", &dssdyid, &b_dssdyid);
   fChain->SetBranchAddress("dssdzpos", &dssdzpos, &b_dssdzpos);
   Notify();
}

Bool_t Al22_Angle::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Al22_Angle::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Al22_Angle::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Al22_Angle_cxx
