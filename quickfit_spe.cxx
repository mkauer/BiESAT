/******************************************************************
A program to fit a low Photo Electron spectra using K.Lang method.

version: 4
build:   08.10.31

change log: (most recent at top)
+ added pedestal parameter options
~ changed this to the quickfit version!!
~ changed the fit to k.lang version!
+ file output to store all parameters
+ changed the file input format

Matt Kauer
******************************************************************/

#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"
#include "_readinData.hxx"

void quickfit_spe(TH1F *histio, Int_t xmax=4000, Int_t xmin=0, Double_t ped=300, Double_t ped_sig=20){ 
  
  beautify();  
  
  TH1F *spehist=(TH1F*)cloneTH1(histio);
  spehist->SetName("SPE Fit");
  spehist->SetTitle("SPE Fit");
  spehist->Rebin(2);
  
  Float_t mean=spehist->GetMean();
  Float_t rms=spehist->GetRMS();
  Int_t pmin=mean-4*rms;
  if(pmin<xmin-20)pmin=xmin;
  Int_t pmax=mean+7*rms;
  if(pmax>xmax-20) pmax=xmax;
  
  TF1 *lang = new TF1("lang",f_spe,xmin,xmax,8);
  lang->SetParNames("Norm","ped_mean","ped_sigma","Npe","1pe_mean","1pe_sigma","Df","Ds");
  lang->SetLineColor(2);
  lang->SetLineWidth(2);
  
  lang->SetParLimits(0,10,10000000);
  lang->SetParLimits(1,ped-(2*ped_sig),ped+(2*ped_sig));
  lang->SetParLimits(2,ped_sig*0.5,ped_sig*2.0);
  lang->SetParLimits(3,0,15);
  lang->SetParLimits(4,1,100);
  lang->SetParLimits(5,1,50);
  lang->SetParLimits(6,0.0,0.2);
  lang->SetParLimits(7,100,100000);
  
  //lang->FixParameter(1,326.92);
  //lang->FixParameter(2,11.02);
  
  TCanvas *c201 = new TCanvas("c201","c201",700,500);
  spehist->Draw();
  spehist->Fit("lang","","",pmin,pmax);
  spehist->SetAxisRange(pmin,pmax,"X");
  c201->Print("AA_newspe-fit.png");
  
  Float_t fit_chi2=lang->GetChisquare();
  Float_t fit_ndf=lang->GetNDF();
  Float_t ped_mean=lang->GetParameter(1);
  Float_t ped_mean_err=lang->GetParError(1);
  Float_t ped_sigma=lang->GetParameter(2);
  Float_t ped_sigma_err=lang->GetParError(2);
  Float_t npe=lang->GetParameter(3);
  Float_t npe_err=lang->GetParError(3);
  Float_t pe_mean=lang->GetParameter(4);
  Float_t pe_mean_err=lang->GetParError(4);
  Float_t pe_sigma=lang->GetParameter(5);
  Float_t pe_sigma_err=lang->GetParError(5);
  
  Float_t ampgain=40;
  Float_t Gain=(pe_mean*25*10e-15)/(1.602176487*10e-19);
  Float_t pe_corr=1+((pe_sigma*pe_sigma)/(pe_mean*pe_mean));
  Float_t Npe_corr=((TMath::Power(mean-ped_mean,2))/((rms*rms)-(ped_sigma*ped_sigma)))*pe_corr;
  
  cout<<"\n\n\t========================================================= ";
  cout<<"\n\t PMT+Amp Gain = "<<Gain;
  cout<<"\n\t Assume pre-amp gain = "<<ampgain;
  cout<<"\n\t PMT Gain = "<<Gain/ampgain;
  cout<<"\n\t Correction Factor = "<<pe_corr;
  cout<<"\n\t Npe from fit = "<<npe;
  cout<<"\n\t Corrected Npe from (mean/rms) method = "<<Npe_corr;
  cout<<"\n\t========================================================= \n\n";
  
}

