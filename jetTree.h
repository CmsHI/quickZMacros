   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[300];   //[nref]
   Float_t         jtpt[300];   //[nref]
   Float_t         jteta[300];   //[nref]
   Float_t         jty[300];   //[nref]
   Float_t         jtphi[300];   //[nref]
   Float_t         jtpu[300];   //[nref]
   Float_t         jtm[300];   //[nref]
   Float_t         discr_fr01[300];   //[nref]
   Float_t         trackMax[300];   //[nref]
   Float_t         trackSum[300];   //[nref]
   Int_t           trackN[300];   //[nref]
   Float_t         trackHardSum[300];   //[nref]
   Int_t           trackHardN[300];   //[nref]
   Float_t         chargedMax[300];   //[nref]
   Float_t         chargedSum[300];   //[nref]
   Int_t           chargedN[300];   //[nref]
   Float_t         chargedHardSum[300];   //[nref]
   Int_t           chargedHardN[300];   //[nref]
   Float_t         photonMax[300];   //[nref]
   Float_t         photonSum[300];   //[nref]
   Int_t           photonN[300];   //[nref]
   Float_t         photonHardSum[300];   //[nref]
   Int_t           photonHardN[300];   //[nref]
   Float_t         neutralMax[300];   //[nref]
   Float_t         neutralSum[300];   //[nref]
   Int_t           neutralN[300];   //[nref]
   Float_t         eMax[300];   //[nref]
   Float_t         eSum[300];   //[nref]
   Int_t           eN[300];   //[nref]
   Float_t         muMax[300];   //[nref]
   Float_t         muSum[300];   //[nref]
   Int_t           muN[300];   //[nref]
   Float_t         discr_ssvHighEff[300];   //[nref]
   Float_t         discr_ssvHighPur[300];   //[nref]
   Float_t         discr_csvMva[300];   //[nref]
   Float_t         discr_csvSimple[300];   //[nref]
   Float_t         discr_muByIp3[300];   //[nref]
   Float_t         discr_muByPt[300];   //[nref]
   Float_t         discr_prob[300];   //[nref]
   Float_t         discr_probb[300];   //[nref]
   Float_t         discr_tcHighEff[300];   //[nref]
   Float_t         discr_tcHighPur[300];   //[nref]
   Float_t         ndiscr_ssvHighEff[300];   //[nref]
   Float_t         ndiscr_ssvHighPur[300];   //[nref]
   Float_t         ndiscr_csvSimple[300];   //[nref]
   Float_t         ndiscr_muByPt[300];   //[nref]
   Float_t         ndiscr_prob[300];   //[nref]
   Float_t         ndiscr_probb[300];   //[nref]
   Float_t         ndiscr_tcHighEff[300];   //[nref]
   Float_t         ndiscr_tcHighPur[300];   //[nref]
   Float_t         pdiscr_csvSimple[300];   //[nref]
   Float_t         pdiscr_prob[300];   //[nref]
   Float_t         pdiscr_probb[300];   //[nref]
   Int_t           nsvtx[300];   //[nref]
   Int_t           svtxntrk[300];   //[nref]
   Float_t         svtxdl[300];   //[nref]
   Float_t         svtxdls[300];   //[nref]
   Float_t         svtxdl2d[300];   //[nref]
   Float_t         svtxdls2d[300];   //[nref]
   Float_t         svtxm[300];   //[nref]
   Float_t         svtxpt[300];   //[nref]
   Int_t           nIPtrk[300];   //[nref]
   Int_t           nselIPtrk[300];   //[nref]
   /*Int_t           nIP;
   Int_t           ipJetIndex[30000];   //[nIP]
   Float_t         ipPt[30000];   //[nIP]
   Float_t         ipProb0[30000];   //[nIP]
   Float_t         ipProb1[30000];   //[nIP]
   Float_t         ip2d[30000];   //[nIP]
   Float_t         ip2dSig[30000];   //[nIP]
   Float_t         ip3d[30000];   //[nIP]
   Float_t         ip3dSig[30000];   //[nIP]
   Float_t         ipDist2Jet[30000];   //[nIP]
   Float_t         ipDist2JetSig[30000];   //[nIP]
   Float_t         ipClosest2Jet[30000];   //[nIP]*/
   Float_t         mue[300];   //[nref]
   Float_t         mupt[300];   //[nref]
   Float_t         mueta[300];   //[nref]
   Float_t         muphi[300];   //[nref]
   Float_t         mudr[300];   //[nref]
   Float_t         muptrel[300];   //[nref]
   Int_t           muchg[300];   //[nref]

void initjetTree(TTree *tree) {

   tree->SetBranchAddress("evt", &evt);
   tree->SetBranchAddress("b", &b);
   tree->SetBranchAddress("nref", &nref);
   tree->SetBranchAddress("rawpt", rawpt);
   tree->SetBranchAddress("jtpt", jtpt);
   tree->SetBranchAddress("jteta", jteta);
   tree->SetBranchAddress("jty", jty);
   tree->SetBranchAddress("jtphi", jtphi);
   tree->SetBranchAddress("jtpu", jtpu);
   tree->SetBranchAddress("jtm", jtm);
   tree->SetBranchAddress("discr_fr01", discr_fr01);
   tree->SetBranchAddress("trackMax", trackMax);
   tree->SetBranchAddress("trackSum", trackSum);
   tree->SetBranchAddress("trackN", trackN);
   tree->SetBranchAddress("trackHardSum", trackHardSum);
   tree->SetBranchAddress("trackHardN", trackHardN);
   tree->SetBranchAddress("chargedMax", chargedMax);
   tree->SetBranchAddress("chargedSum", chargedSum);
   tree->SetBranchAddress("chargedN", chargedN);
   tree->SetBranchAddress("chargedHardSum", chargedHardSum);
   tree->SetBranchAddress("chargedHardN", chargedHardN);
   tree->SetBranchAddress("photonMax", photonMax);
   tree->SetBranchAddress("photonSum", photonSum);
   tree->SetBranchAddress("photonN", photonN);
   tree->SetBranchAddress("photonHardSum", photonHardSum);
   tree->SetBranchAddress("photonHardN", photonHardN);
   tree->SetBranchAddress("neutralMax", neutralMax);
   tree->SetBranchAddress("neutralSum", neutralSum);
   tree->SetBranchAddress("neutralN", neutralN);
   tree->SetBranchAddress("eMax", eMax);
   tree->SetBranchAddress("eSum", eSum);
   tree->SetBranchAddress("eN", eN);
   tree->SetBranchAddress("muMax", muMax);
   tree->SetBranchAddress("muSum", muSum);
   tree->SetBranchAddress("muN", muN);
   tree->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff);
   tree->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur);
   tree->SetBranchAddress("discr_csvMva", discr_csvMva);
   tree->SetBranchAddress("discr_csvSimple", discr_csvSimple);
   tree->SetBranchAddress("discr_muByIp3", discr_muByIp3);
   tree->SetBranchAddress("discr_muByPt", discr_muByPt);
   tree->SetBranchAddress("discr_prob", discr_prob);
   tree->SetBranchAddress("discr_probb", discr_probb);
   tree->SetBranchAddress("discr_tcHighEff", discr_tcHighEff);
   tree->SetBranchAddress("discr_tcHighPur", discr_tcHighPur);
   tree->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff);
   tree->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur);
   tree->SetBranchAddress("ndiscr_csvSimple", ndiscr_csvSimple);
   tree->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt);
   tree->SetBranchAddress("ndiscr_prob", ndiscr_prob);
   tree->SetBranchAddress("ndiscr_probb", ndiscr_probb);
   tree->SetBranchAddress("ndiscr_tcHighEff", ndiscr_tcHighEff);
   tree->SetBranchAddress("ndiscr_tcHighPur", ndiscr_tcHighPur);
   tree->SetBranchAddress("pdiscr_csvSimple", pdiscr_csvSimple);
   tree->SetBranchAddress("pdiscr_prob", pdiscr_prob);
   tree->SetBranchAddress("pdiscr_probb", pdiscr_probb);
   tree->SetBranchAddress("nsvtx", nsvtx);
   tree->SetBranchAddress("svtxntrk", svtxntrk);
   tree->SetBranchAddress("svtxdl", svtxdl);
   tree->SetBranchAddress("svtxdls", svtxdls);
   tree->SetBranchAddress("svtxdl2d", svtxdl2d);
   tree->SetBranchAddress("svtxdls2d", svtxdls2d);
   tree->SetBranchAddress("svtxm", svtxm);
   tree->SetBranchAddress("svtxpt", svtxpt);
   tree->SetBranchAddress("nIPtrk", nIPtrk);
   tree->SetBranchAddress("nselIPtrk", nselIPtrk);
   /*tree->SetBranchAddress("nIP", &nIP);
   tree->SetBranchAddress("ipJetIndex", ipJetIndex);
   tree->SetBranchAddress("ipPt", ipPt);
   tree->SetBranchAddress("ipProb0", ipProb0);
   tree->SetBranchAddress("ipProb1", ipProb1);
   tree->SetBranchAddress("ip2d", ip2d);
   tree->SetBranchAddress("ip2dSig", ip2dSig);
   tree->SetBranchAddress("ip3d", ip3d);
   tree->SetBranchAddress("ip3dSig", ip3dSig);
   tree->SetBranchAddress("ipDist2Jet", ipDist2Jet);
   tree->SetBranchAddress("ipDist2JetSig", ipDist2JetSig);
   tree->SetBranchAddress("ipClosest2Jet", ipClosest2Jet);*/
   tree->SetBranchAddress("mue", mue);
   tree->SetBranchAddress("mupt", mupt);
   tree->SetBranchAddress("mueta", mueta);
   tree->SetBranchAddress("muphi", muphi);
   tree->SetBranchAddress("mudr", mudr);
   tree->SetBranchAddress("muptrel", muptrel);
   tree->SetBranchAddress("muchg", muchg);

}
