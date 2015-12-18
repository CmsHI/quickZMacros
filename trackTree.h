   const static int NTRK = 100000 ; // 53122
   const static int NVTX = 100 ; // 53122
   // Declaration of leaf types

      //if((trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>15 && (trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>pfHcal[i]+pfEcal[i]) continue;} //Calo Matching
   
   Int_t           nEv_;
   Int_t           nLumi_;
   Int_t           nBX_;
   Int_t           nRun_;
   Int_t           N_;
   Int_t           nVtx_;
   Int_t           nTrk_;
   Int_t           maxPtVtx_;
   Int_t           maxMultVtx_;
   Int_t           nTrkVtx_[NVTX];   //[nVtx]
   Float_t         normChi2Vtx_[NVTX];   //[nVtx]
   Float_t         sumPtVtx_[NVTX];   //[nVtx]
   Float_t         xVtx_[NVTX];   //[nVtx]
   Float_t         yVtx_[NVTX];   //[nVtx]
   Float_t         zVtx_[NVTX];   //[nVtx]
   Float_t         xVtxErr_[NVTX];   //[nVtx]
   Float_t         yVtxErr_[NVTX];   //[nVtx]
   Float_t         zVtxErr_[NVTX];   //[nVtx]
   Float_t         vtxDist2DErr_[NVTX];   //[nVtx]
   Float_t         vtxDist2DSig_[NVTX];   //[nVtx]
   Float_t         vtxDist2D_[NVTX];   //[nVtx]
   Float_t         vtxDist3DErr_[NVTX];   //[nVtx]
   Float_t         vtxDist3DSig_[NVTX];   //[nVtx]
   Int_t           nVtxSim_;
   Float_t         xVtxSim_[NVTX];   //[nVtxSim]
   Float_t         yVtxSim_[NVTX];   //[nVtxSim]
   Float_t         zVtxSim_[NVTX];   //[nVtxSim]
   Float_t         trkPt_[NTRK];   //[nTrk]
   Float_t         trkPtError_[NTRK];   //[nTrk]
   UChar_t         trkNHit_[NTRK];   //[nTrk]
   UChar_t         trkNlayer_[NTRK];   //[nTrk]
   Float_t         trkEta_[NTRK];   //[nTrk]
   Float_t         trkPhi_[NTRK];   //[nTrk]
   Int_t           trkCharge_[NTRK];   //[nTrk]
   UChar_t         trkNVtx_[NTRK];   //[nTrk]
   Int_t           nTrkTimesnVtx_;
   Bool_t          trkAssocVtx_[NTRK];   //[nTrkTimesnVtx]
   Bool_t          highPurity_[NTRK];   //[nTrk]
   Bool_t          tight_[NTRK];   //[nTrk]
   Bool_t          loose_[NTRK];   //[nTrk]
   Float_t         trkChi2_[NTRK];   //[nTrk]
   UChar_t         trkNdof_[NTRK];   //[nTrk]
   Float_t         trkDxy1_[NTRK];   //[nTrk]
   Float_t         trkDxyError1_[NTRK];   //[nTrk]
   Float_t         trkDz1_[NTRK];   //[nTrk]
   Float_t         trkDzError1_[NTRK];   //[nTrk]
   Bool_t          trkFake_[NTRK];   //[nTrk]
   UChar_t         trkAlgo_[NTRK];   //[nTrk]
   UChar_t         trkOriginalAlgo_[NTRK];   //[nTrk]
   Float_t         trkMVA_[NTRK];   //[nTrk]
   Int_t           pfType_[NTRK];   //[nTrk]
   Float_t         pfCandPt_[NTRK];   //[nTrk]
   Float_t         pfEcal_[NTRK];   //[nTrk]
   Float_t         pfHcal_[NTRK];   //[nTrk]


   
void initTrackTree(TTree *tree)
{
   tree->SetBranchAddress("nEv", &nEv_);
   tree->SetBranchAddress("nLumi", &nLumi_);
   tree->SetBranchAddress("nBX", &nBX_);
   tree->SetBranchAddress("nRun", &nRun_);
   tree->SetBranchAddress("N", &N_);
   tree->SetBranchAddress("nVtx", &nVtx_);
   tree->SetBranchAddress("nTrk", &nTrk_);
   tree->SetBranchAddress("maxPtVtx", &maxPtVtx_);
   tree->SetBranchAddress("maxMultVtx", &maxMultVtx_);
   tree->SetBranchAddress("nTrkVtx", nTrkVtx_);
   tree->SetBranchAddress("normChi2Vtx", normChi2Vtx_);
   tree->SetBranchAddress("sumPtVtx", sumPtVtx_);
   tree->SetBranchAddress("xVtx", xVtx_);
   tree->SetBranchAddress("yVtx", yVtx_);
   tree->SetBranchAddress("zVtx", zVtx_);
   tree->SetBranchAddress("xVtxErr", xVtxErr_);
   tree->SetBranchAddress("yVtxErr", yVtxErr_);
   tree->SetBranchAddress("zVtxErr", zVtxErr_);
   tree->SetBranchAddress("vtxDist2D", vtxDist2D_);
   tree->SetBranchAddress("vtxDist2DErr", vtxDist2DErr_);
   tree->SetBranchAddress("vtxDist2DSig", vtxDist2DSig_);
   tree->SetBranchAddress("vtxDist3DErr", vtxDist3DErr_);
   tree->SetBranchAddress("vtxDist3DSig", vtxDist3DSig_);
   tree->SetBranchAddress("nVtxSim", &nVtxSim_);
   tree->SetBranchAddress("xVtxSim", &xVtxSim_);
   tree->SetBranchAddress("yVtxSim", &yVtxSim_);
   tree->SetBranchAddress("zVtxSim", &zVtxSim_);
   tree->SetBranchAddress("trkPt", trkPt_);
   tree->SetBranchAddress("trkPtError", trkPtError_);
   tree->SetBranchAddress("trkNHit", trkNHit_);
   tree->SetBranchAddress("trkNlayer", trkNlayer_);
   tree->SetBranchAddress("trkEta", trkEta_);
   tree->SetBranchAddress("trkPhi", trkPhi_);
   tree->SetBranchAddress("trkCharge", trkCharge_);
   tree->SetBranchAddress("trkNVtx", trkNVtx_);
   tree->SetBranchAddress("nTrkTimesnVtx", &nTrkTimesnVtx_);
   tree->SetBranchAddress("trkAssocVtx", trkAssocVtx_);
   tree->SetBranchAddress("highPurity", highPurity_);
   tree->SetBranchAddress("tight", tight_);
   tree->SetBranchAddress("loose", loose_);
   tree->SetBranchAddress("trkChi2", trkChi2_);
   tree->SetBranchAddress("trkNdof", trkNdof_);
   tree->SetBranchAddress("trkDxy1", trkDxy1_);
   tree->SetBranchAddress("trkDxyError1", trkDxyError1_);
   tree->SetBranchAddress("trkDz1", trkDz1_);
   tree->SetBranchAddress("trkDzError1", trkDzError1_);
   tree->SetBranchAddress("trkFake", trkFake_);
   tree->SetBranchAddress("trkAlgo", trkAlgo_);
   tree->SetBranchAddress("trkOriginalAlgo", trkOriginalAlgo_);
   tree->SetBranchAddress("trkMVA", trkMVA_);
   tree->SetBranchAddress("pfType", pfType_);
   tree->SetBranchAddress("pfCandPt", pfCandPt_);
   tree->SetBranchAddress("pfEcal", pfEcal_);
   tree->SetBranchAddress("pfHcal", pfHcal_);
}

