#include "ggTree.h"
#include "jetTree.h"
#include "trackTree.h"

bool goodMuon(int imu) {
  if(muChi2NDF->at(imu)<10
      && muMuonHits->at(imu)>0
      && muStations->at(imu)>1
      && muTrkLayers->at(imu)>5
      && muPixelHits->at(imu)>0
      && fabs(muInnerD0->at(imu))<0.2
      && fabs(muInnerDz->at(imu))<0.5
      ) return true;
  else return false;
}

bool goodElectron(int i) {
  if(fabs(eleSCEta->at(i))<1.479) {
    if(  fabs(eledEtaAtVtx->at(i))<0.0094
        && fabs(eledPhiAtVtx->at(i))<0.0296
        && eleSigmaIEtaIEta->at(i)<0.0101
        && eleHoverE->at(i)<0.0372
        && fabs(eleD0->at(i))<0.0151
        && fabs(eleDz->at(i))<0.238
        && fabs(eleEoverPInv->at(i))<0.118
        && eleMissHits->at(i) <= 2
        ) return true;
    else return false;
  }
  if(fabs(eleSCEta->at(i))>1.479) {
    if(  fabs(eledEtaAtVtx->at(i))<0.00773
        && fabs(eledPhiAtVtx->at(i))<0.148
        && eleSigmaIEtaIEta->at(i)<0.0287
        && eleHoverE->at(i)<0.0546
        && fabs(eleD0->at(i))<0.0535
        && fabs(eleDz->at(i))<0.572
        && fabs(eleEoverPInv->at(i))<0.104
        && eleMissHits->at(i) <= 1
        ) return true;
    else return false;
  }
  else return false;
}

bool goodJet(int i) {
  if(	neutralSum[i]/rawpt[i] < 0.9
      && chargedSum[i]/rawpt[i] > 0.01
      && chargedN[i]+photonN[i]+neutralN[i]+eN[i]+muN[i] > 0
      && chargedMax[i]/rawpt[i] < 0.99
      && photonSum[i]/rawpt[i] < 0.99
      && eSum[i]/rawpt[i] < 0.99
      ) return true;
  else return false;
}

void ggTrack(TString infilename="HiForest.root", TString outfilename="Zevents.root") {

  float leptonptcut = 20;

  TH1::SetDefaultSumw2();

  int hiBin;

  int Ztype; //type 1 muon, type 2 electron
  float Zmass, Zpt, Zeta, Zrapidity, Zphi;
  float Zlepton1Pt, Zlepton2Pt, Zlepton1Eta, Zlepton2Eta, Zlepton1Phi, Zlepton2Phi;
  float weight = 0;
  int Zcharge;


  int njet;
  float jetpt[200], jeteta[200], jetphi[200]; 
  int jetID[200];
  
  Int_t           nTrk;
  Float_t         trkPt[100000];   //[nTrk]
  Float_t         trkPtError[100000];   //[nTrk]
  UChar_t         trkNHit[100000];   //[nTrk]
  UChar_t         trkNlayer[100000];   //[nTrk]
  Float_t         trkEta[100000];   //[nTrk]
  Float_t         trkPhi[100000];   //[nTrk]
  Int_t           trkCharge[100000];   //[nTrk]
  UChar_t         trkNVtx[100000];   //[nTrk]
  Bool_t          highPurity[100000];   //[nTrk]
  Bool_t          tight[100000];   //[nTrk]
  Bool_t          loose[100000];   //[nTrk]
  Float_t         trkChi2[100000];   //[nTrk]
  UChar_t         trkNdof[100000];   //[nTrk]
  Float_t         trkDxy1[100000];   //[nTrk]
  Float_t         trkDxyError1[100000];   //[nTrk]
  Float_t         trkDz1[100000];   //[nTrk]
  Float_t         trkDzError1[100000];   //[nTrk]
  Bool_t          trkFake[100000];   //[nTrk]
  UChar_t         trkAlgo[100000];   //[nTrk]
  UChar_t         trkOriginalAlgo[100000];   //[nTrk]
  Float_t         trkMVA[100000];   //[nTrk]
  Int_t           pfType[100000];   //[nTrk]
  Float_t         pfCandPt[100000];   //[nTrk]
  Float_t         pfEcal[100000];   //[nTrk]
  Float_t         pfHcal[100000];   //[nTrk]

  TTree *ztree = new TTree("ztree","Z boson candidate events");
  ztree->Branch("run",	&run,	"run/I");
  ztree->Branch("event",	&event,	"event/I");
  ztree->Branch("lumis",	&lumis,	"lumis/I");
  ztree->Branch("hiBin", &hiBin, "hiBin/I");
  ztree->Branch("Ztype",	&Ztype,	"Ztype/I");
  ztree->Branch("Zmass",	&Zmass,	"Zmass/F");
  ztree->Branch("Zpt",	&Zpt,	"Zpt/F");
  ztree->Branch("Zeta",	&Zeta,	"Zeta/F");
  ztree->Branch("Zphi",	&Zphi,	"Zphi/F");
  ztree->Branch("Zrapidity",	&Zrapidity,	"Zrapidity/F");
  ztree->Branch("Zlepton1Pt",	&Zlepton1Pt,	"Zlepton1Pt/F");
  ztree->Branch("Zlepton2Pt",	&Zlepton2Pt,	"Zlepton2Pt/F");
  ztree->Branch("Zlepton1Eta",	&Zlepton1Eta,	"Zlepton1Eta/F");
  ztree->Branch("Zlepton2Eta",	&Zlepton2Eta,	"Zlepton2Eta/F");
  ztree->Branch("Zlepton1Phi",	&Zlepton1Phi,	"Zlepton1Phi/F");
  ztree->Branch("Zlepton2Phi",	&Zlepton2Phi,	"Zlepton2Phi/F");
  ztree->Branch("Zcharge",	&Zcharge,	"Zcharge/I");
  ztree->Branch("njet",	&njet,	"njet/I");
  ztree->Branch("jetpt",	&jetpt,	"jetpt[njet]/F");
  ztree->Branch("jeteta",	&jeteta,	"jeteta[njet]/F");
  ztree->Branch("jetphi",	&jetphi,	"jetphi[njet]/F");
  ztree->Branch("jetID",	&jetID,	"jetID[njet]/I");
  ztree->Branch("nTrk",	&nTrk,	"nTrk/I");
  ztree->Branch("trkPt",	&trkPt,	"trkPt[nTrk]/F");
  ztree->Branch("trkPtError", &trkPtError,"trkPtError[nTrk]/F");
  ztree->Branch("trkNHit", &trkNHit,"trkNHit[nTrk]/b");
  ztree->Branch("trkNlayer", &trkNlayer,"trkNlayer[nTrk]/b");
  ztree->Branch("trkEta", &trkEta,"trkEta[nTrk]/F");
  ztree->Branch("trkPhi", &trkPhi,"trkPhi[nTrk]/F");
  ztree->Branch("trkCharge", &trkCharge,"trkCharge[nTrk]/I");
  ztree->Branch("trkNVtx", &trkNVtx,"trkNVtx[nTrk]/b");
  ztree->Branch("highPurity", &highPurity,"highPurity[nTrk]/O");
  ztree->Branch("tight", &tight,"tight[nTrk]/O");
  ztree->Branch("loose", &loose,"loose[nTrk]/O");
  ztree->Branch("trkChi2", &trkChi2,"trkChi2[nTrk]/F");
  ztree->Branch("trkNdof", &trkNdof,"trkNdof[nTrk]/b");
  ztree->Branch("trkDxy1", &trkDxy1,"trkDxy1[nTrk]/F");
  ztree->Branch("trkDxyError1", &trkDxyError1,"trkDxyError1[nTrk]/F");
  ztree->Branch("trkDz1", &trkDz1,"trkDz1[nTrk]/F");
  ztree->Branch("trkDzError1", &trkDzError1,"trkDzError1[nTrk]/F");
  ztree->Branch("trkFake", &trkFake,"trkFake[nTrk]/O");
  ztree->Branch("trkAlgo", &trkAlgo,"trkAlgo[nTrk]/b");
  ztree->Branch("trkOriginalAlgo", &trkOriginalAlgo,"trkOriginalAlgo[nTrk]/b");
  ztree->Branch("trkMVA", &trkMVA,"trkMVA[nTrk]/F");
  ztree->Branch("pfType", &pfType,"pfType[nTrk]/I");
  ztree->Branch("pfCandPt", &pfCandPt,"pfCandPt[nTrk]/F");
  ztree->Branch("pfEcal", &pfEcal,"pfEcal[nTrk]/F");
  ztree->Branch("pfHcal", &pfHcal,"pfHcal[nTrk]/F");
  ztree->Branch("weight", &weight,"weight/F");



  const int nPtBins = 13;
  const double PtBins[nPtBins+1]={0,2.5,5.0,7.5,10.0,12.5,15.0,20,30,40,50,70,100,150};
  const int nYBins = 13;

  TH1F *mmMass = new TH1F("mmMass",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);
  TH1F *mmMassSS = new TH1F("mmMassSS",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);

  TH1F *mmPt = new TH1F("mmPt",";p_{T}^{#mu#mu} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  TH1F *mmY = new TH1F("mmY",";y^{#mu#mu};Events",nYBins,-2.6,2.6);
  TH1F *mmPhi = new TH1F("mmPhi",";#phi^{#mu#mu};Events",20,-pi,pi);

  TH1F *eeMass = new TH1F("eeMass",";m^{ee} (GeV/c^{2});Events",30,60,120);
  TH1F *eeMassSS = new TH1F("eeMassSS",";m^{ee} (GeV/c^{2});Events",30,60,120);

  TH1F *eePt = new TH1F("eePt",";p_{T}^{ee} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  TH1F *eeY = new TH1F("eeY",";y^{ee};Events",nYBins,-2.6,2.6);
  TH1F *eePhi = new TH1F("eePhi",";#phi^{ee};Events",20,-pi,pi);

  int nMuPair = 0;
  int nElPair = 0;

  TFile fin(infilename);

  TFile *fout = new TFile(outfilename,"recreate");

  TTree *inggTree = (TTree*)fin.Get("ggHiNtuplizer/EventTree");
  if(!inggTree){
    cout<<"Could not access gg tree!"<<endl;
    return;
  }
  initggTree(inggTree);

  TTree *injetTree = (TTree*)fin.Get("ak4PFJetAnalyzer/t");
  if(!injetTree) injetTree = (TTree*)fin.Get("akPu4PFJetAnalyzer/t");
  // TTree *injetTree = (TTree*)fin.Get("ak3PFJetAnalyzer/t");
  // if(!injetTree) injetTree = (TTree*)fin.Get("akPu3PFJetAnalyzer/t");
  if(!injetTree){
    cout<<"Could not access jet tree!"<<endl;
    return;
  }
  initjetTree(injetTree);
  
  TTree *evttree = (TTree*)fin.Get("hiEvtAnalyzer/HiTree");
  if(!evttree){
    cout<<"Could not access event tree!"<<endl;
    return;
  }
  evttree->SetBranchAddress("weight", &weight);
  
  TTree * tracktree_                     = (TTree*) fin.Get("anaTrack/trackTree");
  if( tracktree_ == 0 ) tracktree_        = (TTree*) fin.Get("ppTrack/trackTree");
  if(!tracktree_){
    cout<<"Could not access track tree!"<<endl;
    return;
  }
  initTrackTree(tracktree_);
  
  /*
  TTree *inevtTree = (TTree*)fin.Get("hiEvtAnalyzer/HiTree");
  if(!inevtTree){
      cout<<"Could not access tree!"<<endl;
      return;
  }
  inevtTree->SetBranchAddress("hiBin", &hiBin);
  */

  int nEv = inggTree->GetEntries();

  for (int j=0; j<nEv; j++) {

    inggTree->GetEntry(j);
    injetTree->GetEntry(j);
    tracktree_->GetEntry(j);
    evttree->GetEntry(j);
    //inevtTree->GetEntry(j);
    if(j%20000 == 0) cout << "Processing event: " << j << endl;
    bool flagMu = 0; bool flagEle = 0;
    njet = 0;

    TLorentzVector muon1, muon2;
    TLorentzVector ele1, ele2;

    for(int i1 = 0; i1 < nMu; i1++) {

      if(muPt->at(i1)>leptonptcut && fabs(muEta->at(i1))<2.4 && goodMuon(i1)) {

        for(int i2 = i1+1; i2 < nMu; i2++) {

          if(muPt->at(i2)>leptonptcut && fabs(muEta->at(i2))<2.4 && goodMuon(i2)) {

            muon1.SetPtEtaPhiM(muPt->at(i1), muEta->at(i1), muPhi->at(i1), 0.105658);
            muon2.SetPtEtaPhiM(muPt->at(i2), muEta->at(i2), muPhi->at(i2), 0.105658);
            TLorentzVector pair = muon1 + muon2;

            if(pair.M()>60 && pair.M()<120) {

              Ztype = 1;
              Zmass = pair.M();
              Zpt = pair.Pt();
              Zrapidity = pair.Rapidity();
              Zeta = pair.Eta();
              Zphi = pair.Phi();
              Zcharge = muCharge->at(i1) + muCharge->at(i2);
              flagMu = 1;
              Zlepton1Pt = muon1.Pt();
              Zlepton2Pt = muon2.Pt();
              Zlepton1Eta = muon1.Eta();
              Zlepton2Eta = muon2.Eta();
              Zlepton1Phi = muon1.Phi();
              Zlepton2Phi = muon2.Phi();
              if(muCharge->at(i1) != muCharge->at(i2)) {
                mmMass->Fill(pair.M());
                mmPt->Fill(pair.Pt());
                mmY->Fill(pair.Rapidity());
                mmPhi->Fill(pair.Phi());
                nMuPair++;
              }
              else {
                mmMassSS->Fill(pair.M());
              }

            }
          }
        }
      }
    } //end of muon loop

    for(int i1 = 0; i1 < nEle; i1++) {

      if(elePt->at(i1)>leptonptcut && fabs(eleSCEta->at(i1))<2.5 && goodElectron(i1) && (fabs(eleSCEta->at(i1))<1.4442 || fabs(eleSCEta->at(i1))>1.566)) {

        for(int i2 = i1+1; i2 < nEle; i2++) {

          if(elePt->at(i2)>leptonptcut && fabs(eleSCEta->at(i2))<2.5 && goodElectron(i2) && (fabs(eleSCEta->at(i2))<1.4442 || fabs(eleSCEta->at(i2))>1.566)) {

            ele1.SetPtEtaPhiM(elePt->at(i1), eleEta->at(i1), elePhi->at(i1), 0.000511);
            ele2.SetPtEtaPhiM(elePt->at(i2), eleEta->at(i2), elePhi->at(i2), 0.000511);
            TLorentzVector pair = ele1 + ele2;

            if(pair.M()>60 && pair.M()<120) {

              Ztype = 2;
              Zmass = pair.M();
              Zpt = pair.Pt();
              Zrapidity = pair.Rapidity();
              Zeta = pair.Eta();
              Zphi = pair.Phi();
              Zcharge = eleCharge->at(i1) + eleCharge->at(i2);
              flagEle = 1;
              Zlepton1Pt = ele1.Pt();
              Zlepton2Pt = ele2.Pt();
              Zlepton1Eta = ele1.Eta();
              Zlepton2Eta = ele2.Eta();
              Zlepton1Phi = ele1.Phi();
              Zlepton2Phi = ele2.Phi();

              if(eleCharge->at(i1) != eleCharge->at(i2)) {
                eeMass->Fill(pair.M());
                eePt->Fill(pair.Pt());
                eeY->Fill(pair.Rapidity());
                eePhi->Fill(pair.Phi());
                nElPair++;
              }
              else {
                eeMassSS->Fill(pair.M());
              }

            }
          }
        }
      }
    } //end of electron loop
    
    if( flagEle==0 && flagMu==0 ) continue;

    njet = 0;
    for(int ij=0; ij<nref; ij++) {

      if(flagMu) {
        float dR1 = sqrt( (muon1.Eta()-jteta[ij])*(muon1.Eta()-jteta[ij]) + getDphi(muon1.Phi(),jtphi[ij])*getDphi(muon1.Phi(),jtphi[ij]) );
        float dR2 = sqrt( (muon2.Eta()-jteta[ij])*(muon2.Eta()-jteta[ij]) + getDphi(muon2.Phi(),jtphi[ij])*getDphi(muon2.Phi(),jtphi[ij]) );
        if(dR1 < 0.4 || dR2 < 0.4) continue; //jets that are equal to the muon rejected
      }

      if(flagEle) {
        float dR1 = sqrt( (ele1.Eta()-jteta[ij])*(ele1.Eta()-jteta[ij]) + getDphi(ele1.Phi(),jtphi[ij])*getDphi(ele1.Phi(),jtphi[ij]) );
        float dR2 = sqrt( (ele2.Eta()-jteta[ij])*(ele2.Eta()-jteta[ij]) + getDphi(ele2.Phi(),jtphi[ij])*getDphi(ele2.Phi(),jtphi[ij]) );
        if(dR1 < 0.4 || dR2 < 0.4) continue; //jets that are equal to the electron rejected
      }


      if(jtpt[ij]>30) {

        jetpt[njet] = jtpt[ij];
        jeteta[njet] = jteta[ij];
        jetphi[njet] = jtphi[ij];
        jetID[njet] = goodJet(ij);
        njet++;

      }

    } //end of jet loop

    int ntracks = 0;
    for(int i = 0 ; i < nTrk_ ; ++i)
    {
      if((trkMVA_[i]<0.5 && trkMVA_[i]!=-99) || (int)trkNHit_[i]<8 || trkPtError_[i]/trkPt_[i]>0.3 || fabs(trkDz1_[i])/trkDzError1_[i]>3 || fabs(trkDxy1_[i])/trkDxyError1_[i]>3) continue;
      trkPt[ntracks] = trkPt_[i];   //[nTrk]
      trkPtError[ntracks] = trkPtError_[i];   //[nTrk]
      trkNHit[ntracks] = trkNHit_[i];   //[nTrk]
      trkNlayer[ntracks] = trkNlayer_[i];   //[nTrk]
      trkEta[ntracks] = trkEta_[i];   //[nTrk]
      trkPhi[ntracks] = trkPhi_[i];   //[nTrk]
      trkCharge[ntracks] = trkCharge_[i];   //[nTrk]
      trkNVtx[ntracks] = trkNVtx_[i];   //[nTrk]
      highPurity[ntracks] = highPurity_[i];   //[nTrk]
      tight[ntracks] = tight_[i];   //[nTrk]
      loose[ntracks] = loose_[i];   //[nTrk]
      trkChi2[ntracks] = trkChi2_[i];   //[nTrk]
      trkNdof[ntracks] = trkNdof_[i];   //[nTrk]
      trkDxy1[ntracks] = trkDxy1_[i];   //[nTrk]
      trkDxyError1[ntracks] = trkDxyError1_[i];   //[nTrk]
      trkDz1[ntracks] = trkDz1_[i];   //[nTrk]
      trkDzError1[ntracks] = trkDzError1_[i];   //[nTrk]
      trkFake[ntracks] = trkFake_[i];   //[nTrk]
      trkAlgo[ntracks] = trkAlgo_[i];   //[nTrk]
      trkOriginalAlgo[ntracks] = trkOriginalAlgo_[i];   //[nTrk]
      trkMVA[ntracks] = trkMVA_[i];   //[nTrk]
      pfType[ntracks] = pfType_[i];   //[nTrk]
      pfCandPt[ntracks] = pfCandPt_[i];   //[nTrk]
      pfEcal[ntracks] = pfEcal_[i];   //[nTrk]
      pfHcal[ntracks] = pfHcal_[i];   //[nTrk]
      ntracks++;
      //if((trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>15 && (trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>pfHcal[i]+pfEcal[i]) continue;} //Calo Matching
    }
    nTrk=ntracks;
    ztree->Fill();

  } //end of loop over events

  for(int i=0; i<nPtBins; i++) {
    double width = PtBins[i+1]-PtBins[i];
    mmPt->SetBinContent(i+1,mmPt->GetBinContent(i+1)/width);
    mmPt->SetBinError(i+1,mmPt->GetBinError(i+1)/width);
  }

  for(int i=0; i<nPtBins; i++) {
    double width = PtBins[i+1]-PtBins[i];
    eePt->SetBinContent(i+1,eePt->GetBinContent(i+1)/width);
    eePt->SetBinError(i+1,eePt->GetBinError(i+1)/width);
  }

  cout << "Number of muon pairs " << nMuPair << endl;
  cout << "same sign = " << mmMassSS->GetEntries() << endl << endl;

  cout << "Number of electron pairs " << nElPair << endl;
  cout << "same sign = " << eeMassSS->GetEntries() << endl << endl;

  fout->cd();
  ztree->Write();
  // injetTree->Write();
  // trackTree->Write();
  mmMass->Write();
  mmMassSS->Write();
  mmPt->Write();
  mmY->Write();
  mmPhi->Write();
  eeMass->Write();
  eeMassSS->Write();
  eePt->Write();
  eeY->Write();
  eePhi->Write();
  fout->Write();

}

int main(int argc, char *argv[])
{
  if((argc != 3))
  {
    std::cout << "Usage: ./ggTrack.exe <inputfile> <outputfile>" << std::endl;
    return 1;
  }
  ggTrack(argv[1], argv[2]);
  return 0;
}

