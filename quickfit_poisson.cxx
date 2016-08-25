/************************************************************************
POISSON DISTRIBUTION FITTING

version: 4
build:   08.10.31

Matt Kauer
************************************************************************/

#include "_userFuncs.hxx"

void quickfit_poisson(TH1F *histio, Int_t xmax=4000, Int_t xmin=0){
  
  beautify();
  
  TH1F *poiss=(TH1F*)cloneTH1(histio);
  //TH1F *poiss = (TH1F*) histio->Clone();
  poiss->SetName("poisson Fit");
  poiss->SetTitle("poisson Fit");
  poiss->Rebin(2);
  Double_t mean=poiss->GetMean();
  Double_t rms=poiss->GetRMS();
  Int_t pmin=mean-4*rms;
  if(pmin<10)pmin=0;
  Int_t pmax=mean+6*rms;
  if(pmax>3990) pmax=4000;
  
  TF1 *func = new TF1("myfit",f_poisson2,xmin,xmax,5);
  func->SetParNames("pois_norm","pois_mean","pois_factor","gaus_norm","gaus_sig");
  func->SetLineColor(2);
  func->SetLineWidth(2);
  
  func->SetParLimits(0, 1.0, 100000);
  func->SetParLimits(1, 0.1, 1000);
  func->SetParLimits(2, 0.1, 500);
  func->SetParLimits(3, 0.0, 1000);
  func->SetParLimits(4, 0.0001, 10);
  
  // tweak this is the fit is being a whore
  func->SetParameters(1000,10,30,10,2);
  
  TCanvas *c99 = new TCanvas("c99","c99",700,500);
  poiss->Draw();
  poiss->Fit("myfit","WWMR");
  poiss->SetAxisRange(pmin,pmax,"X");
  c99->Print("AA_pois-fit.png");
}
