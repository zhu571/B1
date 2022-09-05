//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 10 14:56:08 2022 by ROOT version 6.26/06
// from TChain t/
//////////////////////////////////////////////////////////

#ifndef DSSD_h
#define DSSD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DSSD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventID;
   Char_t          PName[128];
   Double_t        EDep;
   Char_t          VolNamePre[128];

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_PName;   //!
   TBranch        *b_EDep;   //!
   TBranch        *b_VolNamePre;   //!

   DSSD(TTree *tree=0);
   virtual ~DSSD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DSSD_cxx
DSSD::DSSD(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outputdata_t1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outputdata_t1.root");
      }
      f->GetObject("t",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("t","");
      chain->Add("outputdata_t0.root/t");
      chain->Add("outputdata_t1.root/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

DSSD::~DSSD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DSSD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DSSD::LoadTree(Long64_t entry)
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

void DSSD::Init(TTree *tree)
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

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("PName", PName, &b_PName);
   fChain->SetBranchAddress("EDep", &EDep, &b_EDep);
   fChain->SetBranchAddress("VolNamePre", VolNamePre, &b_VolNamePre);
   Notify();
}

Bool_t DSSD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DSSD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DSSD::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DSSD_cxx
