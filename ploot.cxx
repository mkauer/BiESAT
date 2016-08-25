/*******************************************************************
This program will take two files (the gamma+beta, and gamma) and then
subtract out the pedistal and subtract the gamma from the gamma+beta
to leave you with a plot of just the beta.

version: 4
build:   08.10.31

Matt Kauer
*******************************************************************/

#include "_readinData.hxx"
#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"

void ploot(){
  
  beautify();
  
  //////////////////////////////////////////////////////////////////
  // ----------------  USER INPUTS NEEDED ----------------------- //
  //////////////////////////////////////////////////////////////////
  
  const Int_t daq = 0;           // 0=vme : 1=camac
  const Int_t chan= 1;           // 0-7==lowRes : 8-15==hiRes
  
  char *pdrun   = "run24489";
  char *betarun = "run24490";
  char *gammarun= "run24491";
  const Float_t gamma= 220;      // the 400 kev gamma peak position
  const Float_t beta = 570;      // the 976 kev beta peak position
  const Float_t beta2= 340;
  
  const Int_t rebin  = 2;        // rebinning of the histrograms
  const Float_t Eloss= 14.0;     // energy loss of electrons in air
  
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  
  const Int_t hmin=0;
  if(daq==0) const Int_t hmax=4000;
  if(daq==1) const Int_t hmax=2000;
  
  Float_t calib=(beta-gamma)/(596.-Eloss);  // 976-380 keV
  //Float_t calib=(beta-beta2)/(494.);  // 976-482 keV
  cout<<"\n\t Preliminary Calib ==>  "<<calib<<"\n\n\n";
  
  const Int_t xsize=700,ysize=500;
  Int_t xmin=0,xmax=2000,pmin,pmax;
  if(chan<=7 || daq==1) xmax=150;
  if(chan>=8 && daq==0) xmax=600;
  
  /////////////////////////////////////////////////////////////
  //      FIT THE PEDISTAL WITH BASIC GAUSSIAN
  /////////////////////////////////////////////////////////////
  
  TH1F *hped=(TH1F*)readfile(pdrun,0,chan,daq,"Pedestal");
  TF1 *gauss = new TF1("gauss","gaus",xmin,xmax);
  gauss->SetLineColor(2);
  
  TCanvas *c14 = new TCanvas("c14","ploot ped",xsize,ysize);
  hped->Draw();
  hped->Fit("gauss","Q","",xmin,xmax);
  
  Float_t ped = gauss->GetParameter(1);
  Float_t sped = gauss->GetParameter(2);
  Float_t ped_err = gauss->GetParError(1);
  Float_t sped_err = gauss->GetParError(2);
  
  cout<<"\n\n\t Pedestal Mean = "<<ped<<"  /  Sigma = "<<sped<<"\n\n";
  
  pmin=ped-8*sped;
  pmax=ped+10*sped;
  if(pmin<xmin) pmin=xmin;
  hped->SetAxisRange(pmin,pmax,"X");
  c14->Print("ploot-ped.png");
  
  //////////////////////////////////////////////////////////////////
  //      HISTOGRAM CONVERTING AND REBINNING
  //////////////////////////////////////////////////////////////////
  
  TH1F *hgamma=(TH1F*)readfile(gammarun,ped,chan,daq,"Gammas");
  TH1F *hbeta=(TH1F*)readfile(betarun,ped,chan,daq,"Betas");
  
  hgamma->SetLineColor(8);      
  hbeta->SetLineColor(9);
  
  Int_t imin=640*calib;
  Int_t imax=820*calib;
  //Int_t imin=680;
  //Int_t imax=900;
  
  cout<<"\n\n\t Integrate between  "<<imin<<" - "<<imax<<" \n";
  Float_t beta_area = hbeta->Integral(imin,imax);
  cout<<"\t Area under beta hist = "<<beta_area<<" \n";
  Float_t gamma_area = hgamma->Integral(imin,imax);
  cout<<"\t Area under gamma hist = "<<gamma_area<<" \n";
  Float_t area_ratio = beta_area/gamma_area;
  cout<<"\t Gamma scale factor = "<<area_ratio<<" \n\n\n";
  hgamma->Scale(area_ratio);
  
  hgamma->Rebin(rebin);
  hbeta->Rebin(rebin);
  Int_t bins=hbeta->GetNbinsX();
  
  //////////////////////////////////////////////////////////////////
  //      MIX THE PLOTS AND SUBTRACT AND GAMMAS
  //////////////////////////////////////////////////////////////////
  
  TCanvas *c15 = new TCanvas("c15","ploot mix",xsize,ysize);
  TH1F *mix = (TH1F*)cloneTH1(hgamma,"mix","Betas=blue & Gammas=green");
  mix->SetLineColor(8);
  mix->Draw();
  hbeta->Draw("same");
  mix->SetAxisRange(250*calib,1400*calib,"X");
  c15->Print("ploot-mix.png");
  
  TCanvas *c16 = new TCanvas("c16","ploot subtracted",xsize,ysize);
  TH1F *hsub = new TH1F("hsub","Gamma Subtracted",bins,hmin,hmax);
  hsub->Add(hbeta,hgamma,1,-1);
  hsub->Draw();
  hsub->SetAxisRange(250*calib,1400*calib,"X");
  c16->Print("ploot-sub.png");
  
}

