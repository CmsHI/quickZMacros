//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 16 11:49:15 2015 by ROOT version 6.02/10
// from TTree ztree/Z boson candidate events
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef ztree_h
#define ztree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ztree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   Int_t           hiBin;
   Int_t           Ztype;
   Float_t         Zmass;
   Float_t         weight;
   Float_t         Zpt;
   Float_t         Zeta;
   Float_t         Zphi;
   Float_t         Zrapidity;
   Float_t         Zlepton1Pt;
   Float_t         Zlepton2Pt;
   Float_t         Zlepton1Eta;
   Float_t         Zlepton2Eta;
   Float_t         Zlepton1Phi;
   Float_t         Zlepton2Phi;
   Int_t           Zcharge;
   Int_t           njet;
   Float_t         jetpt[11];   //[njet]
   Float_t         jeteta[11];   //[njet]
   Float_t         jetphi[11];   //[njet]
   Int_t           jetID[11];   //[njet]
   Int_t           nTrk;
   Float_t         trkPt[5000];   //[nTrk]
   Float_t         trkPtError[5000];   //[nTrk]
   UChar_t         trkNHit[5000];   //[nTrk]
   UChar_t         trkNlayer[5000];   //[nTrk]
   Float_t         trkEta[5000];   //[nTrk]
   Float_t         trkPhi[5000];   //[nTrk]
   Int_t           trkCharge[5000];   //[nTrk]
   UChar_t         trkNVtx[5000];   //[nTrk]
   Bool_t          highPurity[5000];   //[nTrk]
   Bool_t          tight[5000];   //[nTrk]
   Bool_t          loose[5000];   //[nTrk]
   Float_t         trkChi2[5000];   //[nTrk]
   UChar_t         trkNdof[5000];   //[nTrk]
   Float_t         trkDxy1[5000];   //[nTrk]
   Float_t         trkDxyError1[5000];   //[nTrk]
   Float_t         trkDz1[5000];   //[nTrk]
   Float_t         trkDzError1[5000];   //[nTrk]
   Bool_t          trkFake[5000];   //[nTrk]
   UChar_t         trkAlgo[5000];   //[nTrk]
   UChar_t         trkOriginalAlgo[5000];   //[nTrk]
   Float_t         trkMVA[5000];   //[nTrk]
   Int_t           pfType[5000];   //[nTrk]
   Float_t         pfCandPt[5000];   //[nTrk]
   Float_t         pfEcal[5000];   //[nTrk]
   Float_t         pfHcal[5000];   //[nTrk]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_Ztype;   //!
   TBranch        *b_Zmass;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_Zpt;   //!
   TBranch        *b_Zeta;   //!
   TBranch        *b_Zphi;   //!
   TBranch        *b_Zrapidity;   //!
   TBranch        *b_Zlepton1Pt;   //!
   TBranch        *b_Zlepton2Pt;   //!
   TBranch        *b_Zlepton1Eta;   //!
   TBranch        *b_Zlepton2Eta;   //!
   TBranch        *b_Zlepton1Phi;   //!
   TBranch        *b_Zlepton2Phi;   //!
   TBranch        *b_Zcharge;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jetpt;   //!
   TBranch        *b_jeteta;   //!
   TBranch        *b_jetphi;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkNVtx;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_tight;   //!
   TBranch        *b_loose;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_trkOriginalAlgo;   //!
   TBranch        *b_trkMVA;   //!
   TBranch        *b_pfType;   //!
   TBranch        *b_pfCandPt;   //!
   TBranch        *b_pfEcal;   //!
   TBranch        *b_pfHcal;   //!

   ztree(std::string thisfilename);
   virtual ~ztree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outfname, std::string tag="");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ztree_cxx
ztree::ztree(std::string thisfilename) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile *f = TFile::Open(thisfilename.data());
  TTree * tree = (TTree*) f->Get("ztree");
  Init(tree);
}

ztree::~ztree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ztree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ztree::LoadTree(Long64_t entry)
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

void ztree::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("Ztype", &Ztype, &b_Ztype);
   fChain->SetBranchAddress("Zmass", &Zmass, &b_Zmass);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("Zpt", &Zpt, &b_Zpt);
   fChain->SetBranchAddress("Zeta", &Zeta, &b_Zeta);
   fChain->SetBranchAddress("Zphi", &Zphi, &b_Zphi);
   fChain->SetBranchAddress("Zrapidity", &Zrapidity, &b_Zrapidity);
   fChain->SetBranchAddress("Zlepton1Pt", &Zlepton1Pt, &b_Zlepton1Pt);
   fChain->SetBranchAddress("Zlepton2Pt", &Zlepton2Pt, &b_Zlepton2Pt);
   fChain->SetBranchAddress("Zlepton1Eta", &Zlepton1Eta, &b_Zlepton1Eta);
   fChain->SetBranchAddress("Zlepton2Eta", &Zlepton2Eta, &b_Zlepton2Eta);
   fChain->SetBranchAddress("Zlepton1Phi", &Zlepton1Phi, &b_Zlepton1Phi);
   fChain->SetBranchAddress("Zlepton2Phi", &Zlepton2Phi, &b_Zlepton2Phi);
   fChain->SetBranchAddress("Zcharge", &Zcharge, &b_Zcharge);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jetpt", jetpt, &b_jetpt);
   fChain->SetBranchAddress("jeteta", jeteta, &b_jeteta);
   fChain->SetBranchAddress("jetphi", jetphi, &b_jetphi);
   fChain->SetBranchAddress("jetID", jetID, &b_jetID);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   fChain->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkNVtx", trkNVtx, &b_trkNVtx);
   fChain->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   fChain->SetBranchAddress("tight", tight, &b_tight);
   fChain->SetBranchAddress("loose", loose, &b_loose);
   fChain->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   fChain->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   fChain->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   fChain->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   fChain->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   fChain->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   fChain->SetBranchAddress("trkFake", trkFake, &b_trkFake);
   fChain->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
   fChain->SetBranchAddress("trkOriginalAlgo", trkOriginalAlgo, &b_trkOriginalAlgo);
   fChain->SetBranchAddress("trkMVA", trkMVA, &b_trkMVA);
   fChain->SetBranchAddress("pfType", pfType, &b_pfType);
   fChain->SetBranchAddress("pfCandPt", pfCandPt, &b_pfCandPt);
   fChain->SetBranchAddress("pfEcal", pfEcal, &b_pfEcal);
   fChain->SetBranchAddress("pfHcal", pfHcal, &b_pfHcal);
   Notify();
}

Bool_t ztree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ztree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ztree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ztree_cxx
