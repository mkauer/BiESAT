/*****************************************************
This uses the same code from the beginning of BiESAT
to normalize and substract this histogramms to leave
you with hist pointer 'sub' to use with the 976keV
fitting macros.

VERSION: 09.07.19

Matt Kauer (kauer@hep.ucl.ac.uk)
*****************************************************/
gROOT->Reset();

#include "./include/fittingFuncs.hxx"
#include "./include/init_memory.cxx"
#include "./include/manip_files.hxx"
#include "./include/manip_params.hxx"
#include "./include/manip_histos.hxx"

#include "./quickfit_976kev.cxx"
#include "./quickfit_landau_976.cxx"

Int_t verbose=1;
Double_t eg_diff=576; // diff between 500keV compton edge and the 976keV electron.

void subtract(string runnumber,Int_t dbtmp=0){
  string workdir="./DBruns/";
  if(dbtmp) workdir="./DBtemp/";
  string configfile=runnumber+".conf";
  rootfile=runnumber+".root";
  
  cout<<"Looking in ==> "<<workdir<<" for config file ==> "<<configfile<<endl;
  gROOT->LoadMacro((workdir+configfile).c_str());
  config();
  
  file=new TFile((workdir+rootfile).c_str(),"READ");
  ptemp=(TH1F*)file->Get(pedrun);
  gtemp=(TH1F*)file->Get(gammarun);
  btemp=(TH1F*)file->Get(betarun);
  
  outfile.open("AA_RESULTS.txt");
  
  beautify();
  
  cout<<"\n\n ===>  "<<describe<<endl;
  cout<<"\n\t+++++++++++++++++++++++++++++++++++++++++++ \n";
  if(lowres)  cout<<"\t  USING FIT MODE  ==>  low-resolution \n";
  if(!lowres) cout<<"\t  USING FIT MODE  ==>  hi-resolution \n";
  if(thick)   cout<<"\t  USING FIT MODE  ==>  thick scint \n";
  if(!thick)  cout<<"\t  USING FIT MODE  ==>  thin scint \n";
  cout<<"\t+++++++++++++++++++++++++++++++++++++++++++ \n";
  
  calib=(beta-gamma)/(eg_diff-Eloss);
  if(calib<=0){
    cout<<"\n\n\t CALIB IS NEGATIVE! ASSUMING ABSOLUTE VALUE! \n\n\n";
    calib=TMath::Abs(calib);
  }
  cout<<"\t  Preliminary Calib    ==>  "<<calib<<"\n";
  cout<<"\t  Assumed Energy Loss  ==>  "<<Eloss<<"\n";
  cout<<"\t+++++++++++++++++++++++++++++++++++++++++++ \n\n\n";
  
  Int_t xmin=0,xmax=maxbin,pmin=0,pmax=maxbin;
  TCanvas *ctmp = new TCanvas("ctmp","ctmp",1,1,10,10);
  ctmp->cd();
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //   FIT THE PEDISTAL WITH BASIC GAUSSIAN
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  
  Float_t ped=0;
  Float_t sped=0;
  Float_t ped_err=0;
  Float_t sped_err=0;
  
  if(ptemp){
    hped=(TH1F*)resizeHist(ptemp,0,maxbin,"hped",pedrun);
    
    TF1 *pedfit = new TF1("pedfit","gaus",xmin,xmax);
    pedfit->SetLineColor(2);
    pedfit->SetLineWidth(2);
    pedfit->SetParLimits(0,0,1e7);
    pedfit->SetParLimits(1,0,1e3);
    pedfit->SetParLimits(2,0,1e2);
        
    TCanvas *c70 = new TCanvas("c70","- Pedistal Fit -",0,0,xsize,ysize);
    c70->cd();
    hped->Draw();
    hped->SetAxisRange(2,1000,"X");
    ped=hped->GetMaximumBin();
    sped=hped->GetRMS();
    xmin=0;
    xmax=ped+(1*sped);
    hped->Fit("pedfit","Q","",xmin,xmax);
    
    ped = pedfit->GetParameter(1);
    sped = pedfit->GetParameter(2);
    xmin=ped-(3*sped);
    xmax=ped+(2*sped);
    hped->Fit("pedfit","QWW","",xmin,xmax);
    
    ped = pedfit->GetParameter(1);
    sped = pedfit->GetParameter(2);
    ped_err = pedfit->GetParError(1);
    sped_err = pedfit->GetParError(2);
    
    cout<<"\n\n\t Pedestal Mean = "<<ped<<"  /  Sigma = "<<sped<<"\n";
    
    pmin=ped-8*sped;
    pmax=ped+10*sped;
    if(pmin<0) pmin=0;
    hped->SetAxisRange(pmin,pmax,"X");
    c70->Update();
    c70->Print("AA_pedistal.png");
    
    ctmp->cd();
    delete pedfit;
  }
  
  const Int_t hmin = 1;
  const Int_t hmax = maxbin-ped-1;
  
  //break;
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //   CONVERTING, SHIFTING, REBINNING AND NORMALIZING
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  
  if(!ptemp){ 
    cout<<"\n\t NO PEDESTAL FILE FOUND ==> using a vaule = "<<ped<<"\n\n";
  }
  
  // CONVERTING
  shiftHist(btemp,ped);
  TH1F *hbeta=(TH1F*)resizeHist(btemp,0,maxbin,"hbeta",betarun);
  shiftHist(gtemp,ped);
  TH1F *hgamma=(TH1F*)resizeHist(gtemp,0,maxbin,"hgamma",gammarun);
  
  TF1 *cfit = new TF1("cfit","gaus",0,maxbin);
  cfit->SetLineColor(2);
  cfit->SetLineWidth(2);
  
  // SHIFTING
  Int_t newrebin=rebin+4;
  if(verbose) cout<<"\n\n\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  btemp->Rebin(newrebin);
  xmin=340*calib;
  xmax=440*calib;
  btemp->SetAxisRange(xmin,xmax,"X");
  if(verbose) cout<<"\t looking for BETA max bin between "<<xmin<<" - "<<xmax<<" bins \n";
  Double_t comp_beta=btemp->GetMaximumBin()*(newrebin);
  btemp->SetAxisRange(0,maxbin,"X");
  xmin=comp_beta-30*calib;
  xmax=comp_beta+30*calib;
  if(verbose) cout<<"\t looking for BETA 500keV compton between "<<xmin<<" - "<<xmax<<" bins \n";
  btemp->Fit("cfit","Q","",xmin,xmax);
  comp_beta=cfit->GetParameter(1);
  
  gtemp->Rebin(newrebin);
  xmin=340*calib;
  xmax=440*calib;
  gtemp->SetAxisRange(xmin,xmax,"X");
  if(verbose) cout<<"\t looking for GAMMA max bin between "<<xmin<<" - "<<xmax<<" bins \n";
  Double_t comp_gamma=gtemp->GetMaximumBin()*(newrebin);
  gtemp->SetAxisRange(0,maxbin,"X");
  xmin=comp_gamma-30*calib;
  xmax=comp_gamma+30*calib;
  if(verbose) cout<<"\t looking for GAMMA 500keV compton between "<<xmin<<" - "<<xmax<<" bins \n";
  gtemp->Fit("cfit","Q","",xmin,xmax);
  comp_gamma=cfit->GetParameter(1);
  if(verbose) cout<<"\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n";
  
  shiftHist(hgamma,comp_gamma-comp_beta);
  
  if(verbose){
    TCanvas *look = new TCanvas("look","look",0,0,xsize,ysize);
    look->Divide(1,2,0.02,0.02,0);
    look->cd(1);
    btemp->Draw();
    look->Update();
    look->cd(2);
    gtemp->Draw();
    look->Update();
    ctmp->cd();
  }
  delete cfit;
  //break;
  
  // REBINNING
  hgamma->Rebin(rebin);
  hbeta->Rebin(rebin);
  
  // NORMALIZING
  Int_t normv=0;   /// if normalization fails, try setting this to 1
  normalizer(hgamma,hbeta,normv);
  
  drawMix(hgamma,hbeta);
  ctmp->cd();
  drawSub(hgamma,hbeta,2);
  ctmp->cd();
  TH1F *sub=(TH1F*)subtract(hgamma,hbeta);

  delete ctmp;  
}

