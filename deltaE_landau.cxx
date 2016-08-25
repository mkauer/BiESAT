/******************************************************************
Trying to get a numerical description of the Bi207 spectra observed
with the "delta-E" trigger setup. It needs to account for the
Landau energy losses through the scintillator.

version: 2
build:   08.10.30

change log: (most recent at top)
------------------------------------------------------------------
30oct + have a working description for the landau energy losses

Matt Kauer (kauer@hep.ucl.ac.uk)
******************************************************************/

#include "_readinData.hxx"
#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"

void deltaE_landau(){
  
  beautify();
  
  //////////////////////////////////////////////////////////////////
  // ----------------  USER INPUTS NEEDED ----------------------- //
  //////////////////////////////////////////////////////////////////
  const Int_t daq  = 0;            // 0=vme : 1=camac
  const Int_t chan = 1;            // 0-7==lowRes : 8-15==hiRes
  
  char *pdrun      = "run66469";   // pedestal run
  char *betarun    = "run66467";   // beta run
  
  char *fitopt[8]  = "";           // global fit options
  Int_t rebin      = 2;            // rebinning of the histrograms
  Double_t El976   = 25;           // mean energy loss at 976 keV
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  
  Double_t El482   = El976*1.25;   // energy lost based on manu simulations
  
  Int_t range;
  if(daq==0) range=4000;
  if(daq==1) range=2000;
  
  const Int_t canx=700,cany=500;
  Int_t fmin=0,fmax=range,pmin=0,pmax=range;
  
  Double_t mean,rms,max1,max2;
  Double_t Gn482,Ge482,Gs482,Gc482,nGe482,nGs482,
    Ge482_err,Gs482_err,Gc482_err,chi2_482,ndf_482;
  Double_t Gn976,Ge976,Gs976,Gc976,nGe976,nGs976,
    Ge976_err,Gs976_err,Gc976_err,chi2_976,ndf_976;
  Double_t quench,calib_err,tot_chi2,tot_ndf;
  
  /////////////////////////////////////////////////////////////
  //     FIT THE PEDISTAL WITH BASIC GAUSSIAN
  /////////////////////////////////////////////////////////////
  
  TH1F *hped=(TH1F*)readfile(pdrun,0,chan,daq,"Pedestal");
  TF1 *gauss = new TF1("gauss","gaus",fmin,fmax);
  gauss->SetLineColor(2);
  gauss->SetLineWidth(2);
  
  TCanvas *l0 = new TCanvas("l0","- Pedistal Fit -",0,0,canx,cany);
  l0->cd();
  hped->Draw();
  if(chan<=7 || daq==1) fmax=150;
  if(chan>=8 && daq==0) fmax=600;
  hped->Fit("gauss","Q","",fmin,fmax);
  
  Double_t ped = gauss->GetParameter(1);
  Double_t sped = gauss->GetParameter(2);
  Double_t ped_err = gauss->GetParError(1);
  Double_t sped_err = gauss->GetParError(2);
  
  cout<<"\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  cout<<" Pedestal Mean = "<<ped<<"  /  Sigma = "<<sped<<"\n";
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
  
  pmin=ped-8*sped;
  pmax=ped+10*sped;
  if(pmin<fmin) pmin=fmin;
  hped->SetAxisRange(pmin,pmax,"X");
  l0->Update();
  l0->Print("deltaE-landau-ped.png");
  
  const Int_t hmin = 1;
  const Int_t hmax = range-ped-1;
  
  //////////////////////////////////////////////////////////////////
  //     FIND THE 482 AND 976 PEAK POSITIONS
  //////////////////////////////////////////////////////////////////
  
  TH1F *hbeta  = (TH1F*)readfile(betarun,ped,chan,daq,"Delta-E Fit");
  hbeta->Rebin(rebin);
  TH1F *hbeta1 = (TH1F*)cloneHist(hbeta,"Find 482 Peak",betarun);
  TH1F *hbeta2 = (TH1F*)cloneHist(hbeta,"Find 976 Peak",betarun);
  hbeta1->Rebin(rebin);
  hbeta2->Rebin(rebin);
  hbeta->Rebin(rebin);
  
  TCanvas *l1=new TCanvas("l1","- Peak Locations -",0,0,canx,cany*1.5);
  l1->Divide(1,2,0.01,0.01,0);
  mean=hbeta->GetMean();
  rms=hbeta->GetRMS();
  
  l1->cd(1);
  hbeta1->SetAxisRange(hmin,Int_t(mean-rms),"X");
  max1=hbeta1->GetMaximumBin();
  max1=max1*(range/hbeta1->GetNbinsX());
  TF1 *gauss1 = new TF1("gauss1","gaus",max1-30,max1+30); ///make this better!
  gauss1->SetLineColor(4);
  hbeta1->Draw();
  hbeta1->Fit("gauss1","QR");
  Ge482=gauss1->GetParameter(1);
  Gs482=gauss1->GetParameter(2);
  hbeta1->Fit("gauss1","Q","",Ge482-Gs482,Ge482+Gs482);
  Ge482=gauss1->GetParameter(1);
  Gs482=gauss1->GetParameter(2);
  hbeta1->SetAxisRange(Ge482-10*Gs482,Ge482+10*Gs482,"X");
  l1->Update();
  
  l1->cd(2);
  hbeta2->SetAxisRange(Int_t(mean-rms),hmax,"X");
  max2=hbeta2->GetMaximumBin();
  max2=max2*(range/hbeta2->GetNbinsX());
  TF1 *gauss2 = new TF1("gauss2","gaus",max2-40,max2+40); ///make this better!
  gauss2->SetLineColor(4);
  hbeta2->Draw();
  hbeta2->Fit("gauss2","QR");
  Ge976=gauss2->GetParameter(1);
  Gs976=gauss2->GetParameter(2);
  hbeta2->Fit("gauss2","Q","",Ge976-Gs976,Ge976+Gs976);
  Ge976=gauss2->GetParameter(1);
  Gs976=gauss2->GetParameter(2);
  hbeta2->SetAxisRange(Ge976-10*Gs976,Ge976+10*Gs976,"X");
  l1->Update();
  l1->Print("deltaE-landau-peaks.png");
  
  calib=(Ge976-Ge482)/(494.-(El976-El482));  // 976-482 keV
  
  cout<<"\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  cout<<" 482keV = "<<max1<<" and 976keV = "<<max2<<"\n";
  cout<<" 482keV = "<<Int_t(Ge482)<<" and 976keV = "<<Int_t(Ge976)<<"\n";
  cout<<" Calib = "<<calib<<" \n";
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    
  //break;
  /////////////////////////////////////////////////////////////
  //     THE MAIN FITTING PROCEDURE
  /////////////////////////////////////////////////////////////
  
  //TF1 *fland = new TF1("fland",f_landau_976,hmin,hmax,8);
  //TF1 *fland = new TF1("fland",f_landau_full,hmin,hmax,15);
  TF1 *fland = new TF1("fland",f_landau_test,hmin,hmax,15);
  
  fland->SetLineColor(2);
  fland->SetLineWidth(2);
  
  fland->SetParName(0, "Gn976");
  fland->SetParName(1, "Ge976");
  fland->SetParName(2, "Gs976");
  fland->SetParName(3, "Gc976");
  fland->SetParName(4, "Ln976");
  fland->SetParName(5, "Le976");
  fland->SetParName(6, "Ls976");
  fland->SetParName(7, "offset");
  fland->SetParName(8, "Gn482");
  fland->SetParName(9, "Ge482");
  fland->SetParName(10,"Gs482");
  fland->SetParName(11,"Gc482");
  fland->SetParName(12,"Ln482");
  fland->SetParName(13,"Le482");
  fland->SetParName(14,"Ls482");
  
  fland->SetParLimits(0,0,1e7);   // 976 keV normal
  fland->SetParameter(1,Ge976);   // 976 keV peak
  fland->SetParLimits(1,Ge976*.8,Ge976*1.2);
  fland->SetParameter(2,Gs976);   // 976 keV sigma
  fland->SetParLimits(2,Gs976*.8,Gs976*1.2);
  fland->SetParameter(3,calib);   // 976 calib
  fland->SetParLimits(3,.5,3);    // 976 calib limits
  fland->FixParameter(3,calib);
  fland->SetParLimits(4,0,1e7);   // 976 landau normal
  fland->SetParameter(5,El976);   // 976 mean energy losses
  fland->SetParLimits(5,5,50);
  fland->SetParameter(6,Gs976);   // 976 landau sigma
  fland->SetParLimits(6,0,100);
  
  fland->SetParameter(7,0);       // offset
  fland->SetParLimits(7,0,50);    // offset
  //fland->FixParameter(7,0);       // offset
  
  fland->SetParLimits(8,0,1e7);   // 482 keV normal
  fland->SetParameter(9,Ge482);   // 482 keV peak
  fland->SetParLimits(9,Ge482*.8,Ge482*1.2);
  fland->SetParameter(10,Gs482);  // 482 keV sigma
  fland->SetParLimits(10,Gs482*.8,Gs482*1.2);
  fland->SetParameter(11,calib);  // 482 calib
  fland->SetParLimits(11,.5,3);   // 482 calib limits
  fland->FixParameter(11,calib);
  
  fland->SetParLimits(12,0,1e7);  // 482 landau normal
  fland->SetParameter(13,El482);  // 482 mean energy losses
  fland->SetParLimits(13,5,70);   //
  fland->SetParameter(14,Gs482);  // 482 landau sigma
  fland->SetParLimits(14,0,100);
    
  hbeta->SetAxisRange(hmin,hmax,"X");
  TCanvas *l2=new TCanvas("l2","- Gaussian + Landau Fit -",canx,0,canx,cany);
  l2->cd();
  hbeta->Draw();
  fmin=200*calib;
  fmax=1200*calib;
  //--------------------
  //fmin=550;
  //fmax=900;
  //--------------------
  //hbeta->Fit("fland","","",fmin,fmax);
  hbeta->Fit("fland","","",100,810);
  hbeta->SetAxisRange(fmin,fmax,"X");
  l2->Update();
  l2->Print("deltaE-landau-fit.png");
    
  Ge976=fland->GetParameter(1);
  Gs976=fland->GetParameter(2);
  Gc976=fland->GetParameter(3);
  El976=fland->GetParameter(5);
  Ge482=fland->GetParameter(9);
  Gs482=fland->GetParameter(10);
  El482=fland->GetParameter(13);
    
  quench = Gc976*((976-El976)-(Ge976/Gc976));
  nGe976=Ge976+quench;
  nGe482=Ge482+quench;
  
  cout<<"\n\n================================================================"<<endl;
  cout<<"\t 482 keV resolution is == "<<Gs482/Ge482*2.35*100<<" % "<<endl;
  cout<<"\t 482 keV landua E loss == "<<El482<<" keV \n"<<endl;
  
  cout<<"\t 976 keV resolution is == "<<Gs976/Ge976*2.35*100<<" % "<<endl;
  cout<<"\t 976 keV landua E loss == "<<El976<<" keV \n"<<endl;
  
  cout<<"\t Quench = "<<quench<<" bins == "<<quench/Gc976<<" keV \n"<<endl;
  
  cout<<"\t 482 keV FWHM @ 1 MeV  == "
      <<(Gs482/nGe482*2.35*100)*TMath::Sqrt(nGe482/Gc976/1000)<<" % "<<endl;
  cout<<"\t 976 keV FWHM @ 1 MeV  == "
      <<(Gs976/nGe976*2.35*100)*TMath::Sqrt(nGe976/Gc976/1000)<<" % "<<endl;
  cout<<"================================================================\n\n"<<endl;
  
}



//--- Fit to the Bi207 976keV and 482keV with Landau tails
Double_t f_landau_test(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  Double_t Ln976=par[4];
  Double_t El976=par[5];
  Double_t Ls976=par[6];
  
  Double_t offset=par[7];
  
  Double_t Gn482=par[8];
  Double_t Ge482=par[9];
  Double_t Gs482=par[10];
  Double_t Gc482=par[11];
  Double_t Ln482=par[12];
  Double_t El482=par[13];
  Double_t Ls482=par[14];
  
  Double_t gaus976=0;
  Double_t land976=0;
  Double_t mpv976=0;
  
  Double_t gaus482=0;
  Double_t land482=0;
  Double_t mpv482=0;
  
  //El482=1.25*El976;  // energy lost based on manu simulations
  //El976=El482*0.75;
  
  mpv482  = Ge482-(El482*Gc482);
  gaus482 = Gn482*(1.52*TMath::Gaus(x,Ge482,Gs482) 
		   + 0.438*TMath::Gaus(x,Ge482+72.1*Gc482,Gs482*TMath::Sqrt(1+73*Gc482/Ge482)) 
		   + 0.147*TMath::Gaus(x,Ge482+84.1*Gc482,Gs482*TMath::Sqrt(1+85*Gc482/Ge482)));
  land482 = Ln482*(TMath::Landau((2*mpv482)-x,mpv482,Ls482));
  
  mpv976  = Ge976-(El976*Gc976);
  gaus976 = Gn976*(7.03*TMath::Gaus(x,Ge976,Gs976) 
		   + 1.84*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976)) 
		   + 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976)));
  land976 = Ln976*(TMath::Landau((2*mpv976)-x,mpv976,Ls976));
  
  return gaus976+land976+offset+gaus482+land482;
}

