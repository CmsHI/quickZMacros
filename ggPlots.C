#include "commonFunctions.h"
#include "ggTree.h"

void setMuHistoStyle(TH1* h) {
 h->SetMarkerColor(kRed);
 h->SetLineColor(kRed);
}

void setElHistoStyle(TH1* h) {
 h->SetMarkerColor(kBlue);
 h->SetLineColor(kBlue);
}

void ggPlots() {

 gStyle->SetOptStat(0);

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

 TFile *fin = TFile::Open("Zevents_pp_ZMM_5Dec.root");

 //dimuon mass
 TCanvas *cmm = new TCanvas();

 TH1F *mumass = (TH1F*)fin->Get("mmMass"); setMuHistoStyle(mumass);
 TH1F *mumassSS = (TH1F*)fin->Get("mmMassSS"); setMuHistoStyle(mumassSS);
 mumassSS->SetMarkerStyle(25);

 mumass->Draw("ep");
 mumassSS->Draw("ep same");

 TLegend *lm = new TLegend(0.16,0.60,0.37,0.72);
 lm->SetBorderSize(0);
 lm->SetTextFont(42);
 lm->SetTextSize(0.04);
 lm->AddEntry(mumass,"Opposite charge","lp");
 lm->AddEntry(mumassSS,"Same charge","lp");
 lm->Draw();

 tx->DrawLatex(0.17,0.86,"Tight muons");
 tx->DrawLatex(0.17,0.80,"p_{T}^{#mu} > 20 GeV");
 tx->DrawLatex(0.17,0.74,"|#eta^{#mu}| < 2.4");

 tx->DrawLatex(0.13,0.94,"Express 262548-263261 (PbPb 5 TeV)");

 tx2->DrawLatex(0.79,0.85,"CMS");
 tx3->DrawLatex(0.70,0.80,"Preliminary");

 //dimuon pt
 TCanvas *cpm = new TCanvas();
 TH1F *mmPt = (TH1F*)fin->Get("mmPt"); setMuHistoStyle(mmPt);
 mmPt->Draw("ep");

 //dimuon rapidity
 TCanvas *cym = new TCanvas();
 TH1F *mmY = (TH1F*)fin->Get("mmY"); setMuHistoStyle(mmY);
 mmY->Draw("ep");

 //dimuon phi
 TCanvas *cfm = new TCanvas();
 TH1F *mmPhi = (TH1F*)fin->Get("mmPhi"); setMuHistoStyle(mmPhi);
 mmPhi->Draw("ep");

 //dielectron mass
 TCanvas *cme = new TCanvas();

 TH1F *eemass = (TH1F*)fin->Get("eeMass"); setElHistoStyle(eemass);
 TH1F *eemassSS = (TH1F*)fin->Get("eeMassSS"); setElHistoStyle(eemassSS);
 eemassSS->SetMarkerStyle(25);

 eemass->Draw("ep");
 eemassSS->Draw("ep same");

 TLegend *le = new TLegend(0.16,0.60,0.37,0.72);
 le->SetBorderSize(0);
 le->SetTextFont(42);
 le->SetTextSize(0.04);
 le->AddEntry(eemass,"Opposite charge","lp");
 le->AddEntry(eemassSS,"Same charge","lp");
 le->Draw();

 tx->DrawLatex(0.17,0.86,"Medium electrons");
 tx->DrawLatex(0.17,0.80,"p_{T}^{e} > 20 GeV");
 tx->DrawLatex(0.17,0.74,"|#eta^{e}| < 2.5");

 tx->DrawLatex(0.13,0.94,"Express 262163 - 262274 (pp 5 TeV)");

 tx2->DrawLatex(0.79,0.85,"CMS");
 tx3->DrawLatex(0.70,0.80,"Preliminary");

 //dielectron pt
 TCanvas *cpe = new TCanvas();
 TH1F *eePt = (TH1F*)fin->Get("eePt"); setElHistoStyle(eePt);
 eePt->Draw("ep");

 //dielectron rapidity
 TCanvas *cye = new TCanvas();
 TH1F *eeY = (TH1F*)fin->Get("eeY"); setElHistoStyle(eeY);
 eeY->Draw("ep");

 //dielectron phi
 TCanvas *cfe = new TCanvas();
 TH1F *eePhi = (TH1F*)fin->Get("eePhi"); setElHistoStyle(eePhi);
 eePhi->Draw("ep");

}
