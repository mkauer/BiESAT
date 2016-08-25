/************************************************************************
POISSON DISTRIBUTION FITTING

version: 5
build:   08.10.31

change log: (most recent at top)
------------------------------------------------------------------
~ changed char[]'s to *char[]'s
+ file output to store all parameters

Matt Kauer
************************************************************************/

#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"
#include "_readinData.hxx"

void poissFit(char *runfile, Int_t ped=0, Int_t xmax=4000, Int_t xmin=0, Bool_t write=0){
  
  beautify();
  
  const Int_t daq=0;    // 0==vme : 1==camac
  const Int_t chan=8;   // 0-7==lowRes : 8-15==hiRes
  const Int_t rebin=2;
  
  TH1F *poiss=new TH1F("poisson Fit",runfile,4000,0,4000);
  readfile(runfile,ped,chan,poiss,daq);
  poiss->Rebin(rebin);
  
  Float_t mean=poiss->GetMean();
  Float_t rms=poiss->GetRMS();
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
  
  // tweak this if the fit is being a whore
  func->SetParameters(1000,30,30,10,2);
  
  TCanvas *c99 = new TCanvas("c99","c99",700,500);
  poiss->Draw();
  poiss->Fit("myfit","WWMR");
  poiss->SetAxisRange(pmin,pmax,"X");
  c99->Print("AA_pois-fit.png");
  
  Float_t fit_chi2=func->GetChisquare();
  Float_t fit_ndf=func->GetNDF();
  Float_t pois_mean=func->GetParameter(1);
  Float_t pois_mean_err=func->GetParError(1);
  
  
  if(write){
    FILE *file;
    file=fopen("AA_results.txt","a");
    fprintf(file,"%.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n",mean,rms,pois_mean,pois_mean_err,fit_chi2,fit_ndf);
    fclose(file);
  }
  
  cout<<"\n\n";
  
  
}
