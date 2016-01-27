#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TArrow.h"

#include <iostream>
#include <vector>

const double pi = 3.14159265358979323846;

using namespace std;

float getDphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if ( dphi > pi )
    dphi = dphi - 2.*pi;
  if ( dphi <= -pi ) 
    dphi = dphi + 2.*pi;
  return dphi;
}

float findNcoll(int hiBin) {

  float w=1;
  const int nbins = 20;
  const float Ncoll[nbins] = {1819, 1433, 1127, 882, 685.2, 526.5, 399.3, 297.5, 217.1, 155.1, 107.9, 73.51, 48.76, 31.46, 19.69, 12.02, 7.042, 3.974, 2.12, 1.164};
  for(int i=0; i<nbins; i++) if(hiBin>=i*(200/nbins) && hiBin<(i+1)*(200/nbins)) w=Ncoll[i];
  return w;

}

float findVertexWeightPbPb(float vz) {

  float w = 1.355*exp(-0.5*((vz-0.6682)/7.729)*((vz-0.6682)/7.729));
  return w;

}
