#include "ggTree.h"

int hiBin;
int Zcharge, Ztype; //type 1 muon, type 2 electron
float Zmass, Zpt, Zeta, Zrapidity, Zphi;
int njet;
float jetpt[50], jeteta[50], jetphi[50]; 
int jetID[50];

void initZtree(TTree *ztree) {
 ztree->SetBranchAddress("hiBin",	&hiBin);
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
     fin = TFile::Open("Zevents_pp_data_dimuon_14Dec.root");
     fMC = TFile::Open("Zevents_pp_MC_ZmumuJet_14Dec.root");
   }
   else {
     fin = TFile::Open("Zevents_pp_data_diphoton_15Dec.root");
     fMC = TFile::Open("Zevents_pp_MC_ZeeJet_14Dec.root");
   }
 }
 else {
   if(mu == 1) {
     fin = TFile::Open("Zevents_PbPb_Express_14Dec.root");
     fMC = TFile::Open("Zevents_PbPb_MC_ZmumuJet_14Dec.root");
   }
   else {
     fin = TFile::Open("Zevents_PbPb_Express_14Dec.root");
     fMC = TFile::Open("Zevents_PbPb_MC_ZeeJet_14Dec.root");
   }
 }



 int nData = 0; int nMC = 0;

 TTree *ztree = (TTree*)fin->Get("ztree");
 initZtree(ztree);

 const int nCentBins = 6;
 const float centBins[nCentBins+1] = {0,20,40,60,80,100,200};
 TH1F *hZcent = new TH1F("hZcent",";centrality bin;Counts",nCentBins,centBins);

 TH1F *hZmass = new TH1F("hZmass",";M^{#mu#mu} (GeV);Counts",30,60,120);
 TH1F *hZpt = new TH1F("hZpt",";p_{T}^{#mu#mu} (GeV);Counts",30,30,180);
 TH1F *hZy = new TH1F("hZy",";y^{#mu#mu};Counts",25,-2.5,2.5);
 TH1F *dPhiZjet = new TH1F("dPhiZjet",";#Delta#phi_{Zjet};Counts",30,0,pi);
 TH1F *dEtaZjet = new TH1F("dEtaZjet",";|#Delta#eta_{Zjet}|;Counts",25,0,5);
 TH2F *dPhidEta = new TH2F("dPhidEta",";|#Delta#eta_{Zjet}|;#Delta#phi_{Zjet}",25,0,5,30,0,pi);
 TH1F *xZjet = new TH1F("xZjet",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 //TH1F *xZjetNear = new TH1F("xZjetNear",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);

 if(mu != 1) {
   hZmass->GetXaxis()->SetTitle("M^{ee} (GeV)");
   hZpt->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
   hZy->GetXaxis()->SetTitle("y^{ee}");
 }

 for(int j=0; j<ztree->GetEntries(); j++) {

   ztree->GetEntry(j);
   if(Zcharge != 0 || Ztype != mu) continue;
   if(!pp) hZcent->Fill(hiBin);

   hZmass->Fill(Zmass);

   if(Zpt < ptcutZ || Zmass < 80 || Zmass > 110) continue;


   hZpt->Fill(Zpt);
   hZy->Fill(Zrapidity);
   nData++;

   for(int ij=0; ij<njet; ij++) {

     if(fabs(jeteta[ij])>2.0) continue;
     //if(pp && jetID[ij] != 1) continue;
     float dphi = getDphi(Zphi,jetphi[ij]);
     dPhiZjet->Fill(fabs(dphi));
     float deta = Zeta - jeteta[ij];
     dEtaZjet->Fill(fabs(deta));
     dPhidEta->Fill(fabs(deta),fabs(dphi));
     float xj = jetpt[ij]/Zpt;
     if(dphi>2*pi/3) xZjet->Fill(xj);
     //else xZjetNear->Fill(xj);

   }

 }

 TTree *ztreeMC = (TTree*)fMC->Get("ztree");
 initZtree(ztreeMC);

 TH1F *hZmassMC = new TH1F("hZmassMC",";M^{#mu#mu} (GeV);Counts",30,60,120);
 TH1F *hZptMC = new TH1F("hZptMC",";p_{T}^{#mu#mu} (GeV);Counts",30,30,180);
 TH1F *hZyMC = new TH1F("hZyMC",";y^{#mu#mu};Counts",18,-2.4,2.4);
 TH1F *dPhiZjetMC = new TH1F("dPhiZjetMC",";#Delta#phi_{Zjet};Counts",30,0,pi);
 TH1F *dEtaZjetMC = new TH1F("dEtaZjetMC",";|#Delta#eta_{Zjet}|;Counts",25,0,5);
 TH2F *dPhidEtaMC = new TH2F("dPhidEtaMC",";|#Delta#eta_{Zjet}|;#Delta#phi_{Zjet}",25,0,5,30,0,pi);
 TH1F *xZjetMC = new TH1F("xZjetMC",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);
 //TH1F *xZjetNearMC = new TH1F("xZjetNearMC",";p_{T}^{jet} / p_{T}^{Z};Counts",10,0,2);

 hZmassMC->SetLineColor(color);
 hZptMC->SetLineColor(color);
 hZyMC->SetLineColor(color);
 dPhiZjetMC->SetLineColor(color);
 dEtaZjetMC->SetLineColor(color);
 xZjetMC->SetLineColor(color);
 //xZjetNearMC->SetLineColor(kCyan);

 if(mu != 1) {
   hZmassMC->GetXaxis()->SetTitle("M^{ee} (GeV)");
   hZptMC->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
   hZyMC->GetXaxis()->SetTitle("y^{ee}");
 }

 for(int j=0; j<ztreeMC->GetEntries(); j++) {

   ztreeMC->GetEntry(j);
   if(Zcharge != 0 || Ztype != mu) continue;

   hZmassMC->Fill(Zmass);

   if(Zpt < ptcutZ || Zmass < 80 || Zmass > 110) continue;


   hZptMC->Fill(Zpt);
   hZyMC->Fill(Zrapidity);
   nMC++;

   for(int ij=0; ij<njet; ij++) {

     if(fabs(jeteta[ij])>2.0) continue;
     //if(pp && jetID[ij] != 1) continue;
     float dphi = getDphi(Zphi,jetphi[ij]);
     dPhiZjetMC->Fill(fabs(dphi));
     float deta = Zeta - jeteta[ij];
     dEtaZjetMC->Fill(fabs(deta));
     dPhidEtaMC->Fill(fabs(deta),fabs(dphi));
     float xj = jetpt[ij]/Zpt;
     if(dphi>2*pi/3) xZjetMC->Fill(xj);
     //else xZjetNearMC->Fill(xj);

   }

 }

 hZmassMC->Scale(double(nData)/double(nMC));
 hZptMC->Scale(double(nData)/double(nMC));
 hZyMC->Scale(double(nData)/double(nMC));
 dPhiZjetMC->Scale(double(nData)/double(nMC));
 dEtaZjetMC->Scale(double(nData)/double(nMC));
 xZjetMC->Scale(double(nData)/double(nMC));
 //xZjetNearMC->Scale(double(nData)/double(nMC));

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
   tx->DrawLatex(0.60,0.86,"p_{T}^{#mu} > 20 GeV");
   tx->DrawLatex(0.60,0.80,"|#eta^{#mu}| < 2.4");
   tx->DrawLatex(0.60,0.74,"60 < m^{#mu#mu} < 120 GeV");
   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 }

 TCanvas *c1 = new TCanvas();
 hZmassMC->Draw("hist");
 hZmass->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c1b = new TCanvas();
 hZmassMC->DrawNormalized("hist");
 hZmass->DrawNormalized("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");

 TCanvas *c2 = new TCanvas();
 hZptMC->Draw("hist");
 hZpt->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c3 = new TCanvas();
 hZyMC->Draw("hist");
 hZy->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));

 TCanvas *c4 = new TCanvas();
 dPhiZjetMC->Draw("hist");
 dPhiZjet->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 tx->DrawLatex(0.17,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.17,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.17,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.17,0.74,"akPu4PF jets");
 tx->DrawLatex(0.17,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.17,0.62,"|#eta^{jet}| < 2.0");
 //tx->DrawLatex(0.17,0.56,"jet ID");

 TCanvas *c5 = new TCanvas();
 dEtaZjetMC->Draw("hist");
 dEtaZjet->Draw("ep same");
 if(pp) tx->DrawLatex(0.13,0.94,"13 pb^{-1} (pp 5 TeV)");
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
 tx->DrawLatex(0.63,0.86,"80 < m^{Z} < 110 GeV");
 tx->DrawLatex(0.63,0.80,Form("p_{T}^{Z} > %g GeV",ptcutZ));
 if(pp) tx->DrawLatex(0.63,0.74,"ak4PF jets");
 else   tx->DrawLatex(0.63,0.74,"akPu4PF jets");
 tx->DrawLatex(0.63,0.68,"p_{T}^{jet} > 30 GeV");
 tx->DrawLatex(0.63,0.62,"|#eta^{jet}| < 2.0");

 TCanvas *c6 = new TCanvas();
 dPhidEta->Draw("colz");

 TCanvas *c7 = new TCanvas();
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
 else   tx->DrawLatex(0.13,0.94,"Express (PbPb 5 TeV)");
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
