/************************************************************************
GAMMA DISTRIBUTION FITTING

version: 4
build:   08.10.31

Matt Kauer
************************************************************************/

#include "_userFuncs.hxx"

void quickfit_gamma(TH1F *histio, Int_t xmax=4000){
  
  beautify();
    
  TCanvas *c99 = new TCanvas("c99","c99",700,500);
  TH1F *dist = (TH1F*) cloneTH1(histio);
  dist->SetName("Distribution Fit");
  dist->SetTitle("Distribution Fit");
  
  TF1 *func = new TF1("myfit",f_gamma,0,xmax,3);
  func->SetParNames("norm","k","theta");
  func->SetLineColor(2);
  func->SetLineWidth(2);
  
  func->SetParLimits(0, 0.00001, 100000);
  func->SetParLimits(1, 0.00001, 1000);
  func->SetParLimits(2, 0.00001, 1000);
  
  dist->Fit("myfit","RM");
   
}

Double_t f_gamma(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  //Int_t k=Int_t(xval[0]/par[2]);
  //Double_t l=par[1]/par[2];
  
  return par[0]*(TMath::Power(x,par[1]-1)*TMath::Exp(-(x/par[2])))/(TMath::Power(par[2],par[1])*TMath::Gamma(par[1]));
}
