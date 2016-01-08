#include "ggTree.h"

int hiBin;
int Zcharge, Ztype; //type 1 muon, type 2 electron
float Zmass, Zpt, Zeta, Zrapidity, Zphi, weight, weightall;
int njet;
float jetpt[50], jeteta[50], jetphi[50]; 
int jetID[50];

void initZtree(TTree *ztree) {
 ztree->SetBranchAddress("hiBin",	&hiBin);
 ztree->SetBranchAddress("weight",	&weight);
 ztree->SetBranchAddress("weightall",	&weightall);
 ztree->SetBranchAddress("Ztype",	&Ztype);
 ztree->SetBranchAddress("Zmass",	&Zmass);
 ztree->SetBranchAddress("Zpt",		&Zpt);
 ztree->SetBranchAddress("Zeta",	&Zeta);
 ztree->SetBranchAddress("Zphi",	&Zphi);
 ztree->SetBranchAddress("Zrapidity",	&Zrapidity);
 ztree->SetBranchAddress("Zcharge",	&Zcharge);
 ztree->SetBranchAddress("njet",	&njet);
 ztree->SetBranchAddress("jetpt",	jetpt);
 ztree->SetBranchAddress("jeteta",	jeteta);
 ztree->SetBranchAddress("jetphi",	jetphi);
 ztree->SetBranchAddress("jetID",	jetID);
}

void ggZjetPlots(bool pp = 1, int mu = 1, float ptcutZ = 40) {

 int color = kRed;
 if(mu == 2) color = kBlue;

 gStyle->SetOptStat(0);
 TH1::SetDefaultSumw2();

 TLatex *tx = new TLatex();
 tx->SetTextFont(42);
 tx->SetTextSize(0.04);
 tx->SetNDC(kTRUE);

 TLatex *tx2 = new TLatex();
 tx2->SetTextFont(62);
 tx2->SetTextSize(0.04);
 tx2->SetNDC(kTRUE);

 TLatex *tx3 = new TLatex();
 tx3->SetTextFont(52);
 tx3->SetTextSize(0.04);
 tx3->SetNDC(kTRUE);

 TFile *fin, *fMC;
 if(pp) {
   if(mu == 1) {
     fin = TFile::Open("Zevents_pp_data_dimuon_6Jan.root");
     fMC = TFile::Open("Zevents_pp_MC_ZmumuJet_6Jan.root");
   }
   else {
     fin = TFile::Open("Zevents_pp_data_dielectron_6Jan.root");
     fMC = TFile::Open("Zevents_pp_MC_ZeeJet_6Jan.root");
   }
 }
 else {
   if(mu == 1) {
     fin = TFile::Open("Zevents_PbPb_data_ZMM_7Jan.root");
     fMC = TFile::Open("Zevents_PbPb_MC_ZmumuJet_8Jan.root");
   }
   else {
     fin = TFile::Open("Zevents_PbPb_data_ZEE_7Jan.root");
     fMC = TFile::Open("Zevents_PbPb_MC_ZeeJet_8Jan.root");
   }
 }



 int nData = 0; float nMC = 0;

 TTree *ztree = (TTree*)fin->Get("ztree");
 initZtree(ztree);

 const int nCentBins = 6;
 const float centBins[nCentBins+1] = {0,20,40,60,80,100,200};
 TH1F *hZcent = new TH1F("hZcent",";centrality bin;Counts",nCentBins,centBins);

 TH1F *hZmass = new TH1F("hZmass",";M^{#mu#mu} (GeV);Counts",30,60,120);
 TH1F *hZpt = new TH1F("hZpt",";p_{T}^{#mu#mu} (GeV);Counts",30,30,180);
 TH1F *hZy = new TH1F("hZy",";y^{#mu#mu};Counts",18,-2.4,2.4);
 TH1F *hZphi = new TH1F("hZphi",";#phi^{#mu#mu};Counts",18,-pi,pi);
 TH1F *dPhiZjet = new TH1F("dPhiZjet",";#Delta#phi_{Zjet};Counts",30,0,pi);
 TH1F *dEtaZjet = new TH1F("dEtaZjet",";|#Delta#eta_{Zjet}|;Counts",25,0,5);
 TH2F *dPhidEta = new TH2F("dPhidEta",";|#Delta#eta_{Zjet}|;#Delta#phi_{Zjet}",25,0,5,30,0,pi);
 TH1F *xZjet = new TH1F("xZjet",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 //TH1F *xZjetNear = new TH1F("xZjetNear",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 TH1F *nJetZcut = new TH1F("nJetZcut",";number of jets;Counts",10,0,10);
 TH1F *hJetPt = new TH1F("hJetPt",";jet p_{T} (GeV);Counts",30,30,180);

 if(mu != 1) {
   hZmass->GetXaxis()->SetTitle("M^{ee} (GeV)");
   hZpt->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
   hZy->GetXaxis()->SetTitle("y^{ee}");
 }

 for(int j=0; j<ztree->GetEntries(); j++) {

   ztree->GetEntry(j);
   if(Zcharge != 0 || Ztype != mu) continue;
   if(!pp) hZcent->Fill(hiBin);

   if(Zmass < 80 || Zmass > 110) continue;
   if(Zpt < ptcutZ) continue;

   hZmass->Fill(Zmass);
   hZpt->Fill(Zpt);
   hZy->Fill(Zrapidity);
   hZphi->Fill(Zphi);
   nData++;
   int nJets = 0;

   for(int ij=0; ij<njet; ij++) {

     if(fabs(jeteta[ij])>2.0) continue;
     //if(pp && jetID[ij] != 1) continue;
     float dphi = getDphi(Zphi,jetphi[ij]);
     dPhiZjet->Fill(fabs(dphi));
     float deta = Zeta - jeteta[ij];
     dEtaZjet->Fill(fabs(deta));
     dPhidEta->Fill(fabs(deta),fabs(dphi));
     float xj = jetpt[ij]/Zpt;
     if(dphi>2*pi/3) {
       xZjet->Fill(xj);
     }
     //else xZjetNear->Fill(xj);
     hJetPt->Fill(jetpt[ij]);
     nJets++;

   }

   nJetZcut->Fill(nJets);

 }

 TTree *ztreeMC = (TTree*)fMC->Get("ztree");
 initZtree(ztreeMC);

 TH1F *hZmassMC = new TH1F("hZmassMC",";M^{#mu#mu} (GeV);Counts",30,60,120);
 TH1F *hZptMC = new TH1F("hZptMC",";p_{T}^{#mu#mu} (GeV);Counts",30,30,180);
 TH1F *hZyMC = new TH1F("hZyMC",";y^{#mu#mu};Counts",18,-2.4,2.4);
 TH1F *hZphiMC = new TH1F("hZphiMC",";#phi^{#mu#mu};Counts",18,-pi,pi);
 TH1F *dPhiZjetMC = new TH1F("dPhiZjetMC",";#Delta#phi_{Zjet};Counts",30,0,pi);
 TH1F *dEtaZjetMC = new TH1F("dEtaZjetMC",";|#Delta#eta_{Zjet}|;Counts",25,0,5);
 TH2F *dPhidEtaMC = new TH2F("dPhidEtaMC",";|#Delta#eta_{Zjet}|;#Delta#phi_{Zjet}",25,0,5,30,0,pi);
 TH1F *xZjetMC = new TH1F("xZjetMC",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 //TH1F *xZjetNearMC = new TH1F("xZjetNearMC",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 TH1F *nJetZcutMC = new TH1F("nJetZcutMC",";number of jets;Counts",10,0,10);
 TH1F *hJetPtMC = new TH1F("hJetPtMC",";jet p_{T} (GeV);Counts",30,30,180);

 hZmassMC->SetLineColor(color);
 hZptMC->SetLineColor(color);
 hZyMC->SetLineColor(color);
 hZphiMC->SetLineColor(color);
 dPhiZjetMC->SetLineColor(color);
 dEtaZjetMC->SetLineColor(color);
 xZjetMC->SetLineColor(color);
 //xZjetNearMC->SetLineColor(kCyan);
 nJetZcutMC->SetLineColor(color);
 hJetPtMC->SetLineColor(color);

 if(mu != 1) {
   hZmassMC->GetXaxis()->SetTitle("M^{ee} (GeV)");
   hZptMC->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
   hZyMC->GetXaxis()->SetTitle("y^{ee}");
   hZphiMC->GetXaxis()->SetTitle("#phi^{ee}");
 }

 for(int j=0; j<ztreeMC->GetEntries(); j++) {

   ztreeMC->GetEntry(j);
   float w = 1;
   if(pp) w = weight;
   else w = weightall;

   if(Zcharge != 0 || Ztype != mu) continue;

   if(Zmass < 80 || Zmass > 110) continue;
   if(Zpt < ptcutZ) continue;

   hZmassMC->Fill(Zmass,w);
   hZptMC->Fill(Zpt,w);
   hZyMC->Fill(Zrapidity,w);
   hZphiMC->Fill(Zphi,w);
   nMC+= w;
   int nJets = 0;

   for(int ij=0; ij<njet; ij++) {

     if(fabs(jeteta[ij])>2.0) continue;
     //if(pp && jetID[ij] != 1) continue;
     float dphi = getDphi(Zphi,jetphi[ij]);
     dPhiZjetMC->Fill(fabs(dphi),w);
     float deta = Zeta - jeteta[ij];
     dEtaZjetMC->Fill(fabs(deta),w);
     dPhidEtaMC->Fill(fabs(deta),fabs(dphi),w);
     float xj = jetpt[ij]/Zpt;
     if(dphi>2*pi/3) {
       xZjetMC->Fill(xj,w);
     }
     //else xZjetNearMC->Fill(xj,w);
     hJetPtMC->Fill(jetpt[ij],w);
     nJets++;

   }

   nJetZcutMC->Fill(nJets,w);

 }

 hZmassMC->Scale(float(nData)/float(nMC));
 hZptMC->Scale(float(nData)/float(nMC));
 hZyMC->Scale(float(nData)/float(nMC));
 hZphiMC->Scale(float(nData)/float(nMC));
 dPhiZjetMC->Scale(float(nData)/float(nMC));
 dEtaZjetMC->Scale(float(nData)/float(nMC));
 xZjetMC->Scale(float(nData)/float(nMC));
 //xZjetNearMC->Scale(float(nData)/float(nMC));
 nJetZcutMC->Scale(float(nData)/float(nMC));
 hJetPtMC->Scale(float(nData)/float(nMC));

 if(!pp) {
   float Npart[nCentBins] = {358.8,264.3,189.2,131.4,86.94,21.85};
   float Ncoll[nCentBins] = {1626.0,1005.0,606.4,348.2,186.2,30.74};
   float x[nCentBins], xerr[nCentBins], y[nCentBins], yerr[nCentBins];
   for(int i=0; i<nCentBins; i++) {
     x[i] = Npart[nCentBins-1-i];
     xerr[i] = 0;
     float binwidth = centBins[nCentBins-i]-centBins[nCentBins-1-i];
     y[i] = hZcent->GetBinContent(nCentBins-i)/Ncoll[nCentBins-1-i]/binwidth;
     yerr[i] = hZcent->GetBinError(nCentBins-i)/Ncoll[nCentBins-1-i]/binwidth;
   }
   TCanvas *c0 = new TCanvas();
   TGraphErrors *g = new TGraphErrors(nCentBins,x,y,xerr,yerr);
   g->SetMarkerColor(color);
   g->SetLineColor(color);
   g->GetXaxis()->SetTitle("N_{part}");
   g->GetYaxis()->SetTitle("dN^{Z} / N_{coll}");
   g->Draw("aep");
   if(mu == 1) {
     tx->DrawLatex(0.60,0.86,"p_{T}^{#mu} > 20 GeV");
     tx->DrawLatex(0.60,0.80,"|#eta^{#mu}| < 2.4");
     //tx->DrawLatex(0.60,0.74,"60 < m^{#mu#mu} < 120 GeV");
     tx->DrawLatex(0.60,0.74,"80 < m^{#mu#mu} < 110 GeV");
   }
   else {
     tx->DrawLatex(0.60,0.86,"p_{T}^{e} > 20 GeV");
     tx->DrawLatex(0.60,0.80,"|#eta^{e}| < 2.5");
     tx->DrawLatex(0.60,0.74,"60 < m^{ee} < 120 GeV");
   }
   tx->DrawLatex(0.13,0.94,"180 #mub^{-1} (PbPb 5 TeV)");
 }

 TCanvas *c1 = new TCanvas();
 hZmassMC->Draw("hist");
 hZmass->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c2 = new TCanvas();
 hZptMC->Draw("hist");
 hZpt->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c3 = new TCanvas();
 hZyMC->Draw("hist");
 hZy->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c3b = new TCanvas();
 hZphiMC->Draw("hist");
 hZphi->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c4 = new TCanvas();
 nJetZcutMC->Draw("hist");
 nJetZcut->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.61,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.61,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.61,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.61,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.61,0.74,"akPu4PF jets");
 tx->DrawLatex(0.61,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.61,0.62,"|#eta^{jet}| < 2.0");

 TCanvas *c4b = new TCanvas();
 hJetPtMC->Draw("hist");
 hJetPt->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.61,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.61,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.61,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.61,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.61,0.74,"akPu4PF jets");
 tx->DrawLatex(0.61,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.61,0.62,"|#eta^{jet}| < 2.0");

 TCanvas *c5 = new TCanvas();
 dPhiZjetMC->Draw("hist");
 dPhiZjet->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.17,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.17,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.17,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.17,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.17,0.74,"akPu4PF jets");
 tx->DrawLatex(0.17,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.17,0.62,"|#eta^{jet}| < 2.0");
 //tx->DrawLatex(0.17,0.56,"jet ID");

 TCanvas *c6 = new TCanvas();
 dEtaZjetMC->Draw("hist");
 dEtaZjet->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.63,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.63,0.74,"akPu4PF jets");
 tx->DrawLatex(0.63,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.63,0.62,"|#eta^{jet}| < 2.0");

 TCanvas *c7a = new TCanvas();
 dPhidEta->Draw("colz");

 TCanvas *c7b = new TCanvas();
 dPhidEtaMC->Draw("colz");

 TArrow *arrow = new TArrow();
 arrow->SetLineColor(kBlack);
 arrow->SetFillColor(kBlack);
 arrow->SetLineWidth(2);
 arrow->SetArrowSize(0.03);

 TArrow *arrowMC = (TArrow*)arrow->Clone();
 arrowMC->SetLineColor(color);
 arrowMC->SetFillColor(color);

 TArrow *arrow2 = (TArrow*)arrow->Clone();
 arrow2->SetLineColor(kBlue);
 arrow2->SetFillColor(kBlue);

 TCanvas *c8 = new TCanvas();
 xZjetMC->Draw("hist");
 xZjet->Draw("ep same");
 //xZjetNear->SetMarkerColor(kBlue);
 //xZjetNear->SetLineColor(kBlue);
 //xZjetNear->SetMarkerStyle(24);
 //xZjetNear->Draw("same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"317 #mub^{-1} (PbPb 5 TeV)");
 //tx->DrawLatex(0.63,0.86,"60 < m^{Z} < 120 GeV");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.63,0.74,"ak4PF jets");
 else tx->DrawLatex(0.63,0.74,"akPu4PF jets");
 tx->DrawLatex(0.63,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.63,0.62,"#Delta#phi_{Zjet} > 2#pi/3");
 arrow->DrawArrow(xZjet->GetMean(),xZjet->GetMaximum()*0.2,xZjet->GetMean(),0);
 arrowMC->DrawArrow(xZjetMC->GetMean(),xZjet->GetMaximum()*0.2,xZjetMC->GetMean(),0);
 //arrow2->DrawArrow(xZjetNear->GetMean(),xZjet->GetMaximum()*0.2,xZjetNear->GetMean(),0);

}
