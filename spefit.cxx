/******************************************************************
A program to fit a low Photo Electron spectra using K.Lang method.

version: 4
build:   08.10.31

change log: (most recent at top)
+ added pedsig input and parameter limits
+ added text output to hist for gain and correction
+ changed the fit to k.lang version!
+ file output to store all parameters
+ changed the file input format

Matt Kauer
******************************************************************/
#include "include/manip_histos.hxx"
#include "include/fittingFuncs.hxx"
#include "include/manip_files.hxx"

void spefit(char *fname, Int_t xmax=3000, Int_t xmin=0){ 
  
  beautify();  
  
  Bool_t calcu=true;   // display the calculated gain and correction
  
  const Int_t ped=0;
  const Int_t daq=0;    // 0==vme : 1==camac
  const Int_t chan=8;   // 0-7==lowRes : 8-15==hiRes
  const Int_t rebin=2;
  
  const Float_t pedmean=300;
  const Float_t pedsig=20;
  const Float_t ampgain=40.0;
  
  TCanvas *c201 = new TCanvas("c201","c201",700,500);
  //TH1F *spehist = new TH1F("spehist",fname,4000,0,4000);
  TH1F *spehist;
  spehist=readfile(fname,ped,chan,daq,"spehist");
  spehist->Rebin(rebin);
  
  Float_t mean=spehist->GetMean();
  Float_t rms=spehist->GetRMS();
  Int_t pmin=mean-4*rms;
  if(pmin<10)pmin=0;
  Int_t pmax=mean+7*rms;
  if(pmax>3990) pmax=4000;
  
  TF1 *lang = new TF1("lang",f_spe,xmin,xmax,8);
  lang->SetParNames("Norm","ped_mean","ped_sigma","Npe","1pe_mean","1pe_sigma","Df","Ds");
  lang->SetLineColor(2);
  lang->SetLineWidth(2);
  
  lang->SetParameter(1,pedmean);  // pedestal mean
  lang->SetParameter(2,pedsig);   // pedestal sigma
  lang->SetParameter(4,10);       // 1st PE mean
  lang->SetParameter(5,5);        // 1st PE sigma
  lang->SetParameter(6,1e-4);
  lang->SetParameter(7,900);
  //lang->SetParameter(6,0.1);
  //lang->SetParameter(7,1e5);
  
  
  lang->SetParLimits(0,10,1e7);
  lang->SetParLimits(1,pedmean*0.9,pedmean*1.1); // pedestal mean
  lang->SetParLimits(2,pedsig*0.8,pedsig*1.3);   // pedestal sigma
  lang->SetParLimits(3,0.1,10); //was 0.1,10                  // Number PE
  lang->SetParLimits(4,50,100); //was 1,300                  // 1st PE mean
  lang->SetParLimits(5,1,150); // was 1,150                  // 1st PE sigma
  lang->SetParLimits(6,0,0.2);
  lang->SetParLimits(7,100,1e6);
  
  //lang->FixParameter(1,pedmean);   // pedestal mean
  //lang->FixParameter(2,pedsig);    // pedestal sigma
  //lang->FixParameter(3,2);         // Number PE
  //lang->FixParameter(4,13);        // 1st PE mean
  //lang->FixParameter(5,24);        // 1st PE sigma
  //lang->FixParameter(6,9e-5);
  //lang->FixParameter(7,900);
  
  
  spehist->Draw();
  spehist->Fit("lang","","",pmin,pmax);
  spehist->SetAxisRange(pmin,pmax,"X");
  
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
  
  Float_t Gain=(pe_mean*25*10e-15)/(1.602176487*10e-19);
  Float_t pe_corr=1+((pe_sigma*pe_sigma)/(pe_mean*pe_mean));
  Float_t Npe_corr=((TMath::Power(mean-ped_mean,2))/((rms*rms)-(ped_sigma*ped_sigma)))*pe_corr;
  
  ostringstream oss_gain,oss_corr;
  oss_gain<<"PMT Gain  =  "<<Gain/ampgain;
  oss_corr<<"Corr. Factor  =  "<<pe_corr;
  
  if(calcu){
    Double_t xpos=ped_mean+pe_mean*2;
    Double_t ypos=spehist->GetBinContent(spehist->GetMaximumBin());
    Float_t fsize=0.040;
    TLatex *tex_v1 = new TLatex(xpos,ypos*0.95,oss_gain.str().c_str());
    tex_v1->SetTextSize(fsize);
    tex_v1->Draw();
    TLatex *tex_v2 = new TLatex(xpos,ypos*0.85,oss_corr.str().c_str());
    tex_v2->SetTextSize(fsize);
    tex_v2->Draw();
  }
  
  ostringstream output;
  output<<fname<<".png";
  
  //c201->Print(output.str().c_str());
  c201->Print("spefit-plot.png");
  
  cout<<"\n\n\t========================================================= ";
  cout<<"\n\t PMT+Amp Gain = "<<Gain;
  cout<<"\n\t Assume pre-amp gain = "<<ampgain;
  cout<<"\n\t PMT Gain = "<<Gain/ampgain;
  cout<<"\n\t Correction Factor = "<<pe_corr;
  cout<<"\n\t Npe from fit = "<<npe;
  cout<<"\n\t Corrected Npe from (mean/rms) method = "<<Npe_corr;
  cout<<"\n\t========================================================= \n\n\n";
  
  
  delete lang;
}
