/******************************************************************
This seems to generally work for all bi207 spectrum as long as
there isn't anything funny going on, i.e.
- This can't handle changing the x-axis into units of keV or MeV yet.
     (without some bugs)
- This can't handle weird gain shifts during the run that change
  the calculated "calib" parameter.
- If the offset is too high and the 500keV compton gets cut short,
  this function doesn't know what to do.

version: 4
build:   08.10.31

change log: (most recent at top)
~ changed "h976" to "hbeta" to make compatible with "ntupleFit.cxx"
+ added a 3rd fit loop to correct differences in calib
~ changed the offset calculation because it was wrong
+ added in offset calculation
~ moved if(outfile.isopen) to beginning

Matt Kauer
******************************************************************/

#include "_readinData.hxx"
#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"

void deltaE_gaus(){
  
  beautify();
  
  //////////////////////////////////////////////////////////////////
  // ----------------  USER INPUTS NEEDED ----------------------- //
  //////////////////////////////////////////////////////////////////
  
  const Int_t daq=0;             // 0=vme : 1=camac
  const Int_t chan=1;            // 0-7==lowRes : 8-15==hiRes
  
  char *pdrun="run66469";
  char *betarun="run66467";
  Float_t Ek482=340;             // the 482 kev beta peak position
  Float_t Ek976=740;             // the 976 kev beta peak position
  Float_t s482=20;
  Float_t s976=50;
  
  char *fitopt     = "";       // global fit options
  const Int_t rebin  = 2;        // rebinning of the histrograms
  const Float_t Eloss= 30.0;     // energy loss of electrons in air
  
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  
  if(daq==0) const Int_t bins=4000;
  if(daq==1) const Int_t bins=2000;
  
  Float_t calib=(Ek976-Ek482)/494.;  // 976-482 keV
     
  cout<<"Preliminary calibration = "<<1/calib<<" keV/ch"<<endl;
  cout<<"\t--> Calib = "<<calib<<endl;
  
  const Int_t xsize=700,ysize=500;
  Int_t xmin=0,xmax=2000,pmin,pmax;
  if(chan<=7 || daq==1) xmax=150;
  if(chan>=8 && daq==0) xmax=600;
  
  Float_t Nbi500,cal482,nEk482,ns482,Ek482_err,s482_err,cal482_err,chi2_482,ndf_482;
  Float_t Nbi1000,cal976,nEk976,ns976,Ek976_err,sys_s976,s976_err,cal976_err,chi2_976,ndf_976;
  Float_t offset,calib_err,tot_chi2,tot_ndf;
  
  /////////////////////////////////////////////////////////////
  //   FIT THE PEDISTAL WITH BASIC GAUSSIAN
  /////////////////////////////////////////////////////////////
  
  TH1F *hped=(TH1F*)readfile(pdrun,0,chan,daq,"Pedestal");
  TF1 *gauss = new TF1("gauss","gaus",xmin,xmax);
  gauss->SetLineColor(2);
  
  TCanvas *c70 = new TCanvas("c70","- Pedistal Fit -",xsize,ysize);
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
  c70->Print("A_pedistal.png");
  
  const Float_t hmin = 200.*calib;
  const Float_t hmax = bins-ped;
  
  //////////////////////////////////////////////////////////////////
  //      HISTOGRAM CONVERTING AND SCALING
  //////////////////////////////////////////////////////////////////
  
  TH1F *hbeta=(TH1F*)readfile(betarun,ped,chan,daq,"Rough Fit 1");
  hbeta->Rebin(rebin);
  
  //break;
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  ////   SHOULD BE ABLE TO COPY AND PASTE FROM HERE DOWN   ////
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  
  ofstream outfile;
  outfile.open("A_RESULTS.txt");
  if(outfile.is_open()){
    
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    //      1 ROUGH FIT THE BI207 SPECTRUM
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
  
    TH1F *hbeta1 = (TH1F*)cloneTH1(hbeta,"Rough Fit 2",betarun);
    TH1F *hbeta2 = (TH1F*)cloneTH1(hbeta,"Final Fit",betarun);
    
    TF1 *func1 = new TF1("deltaE",deltaE,hmin,hmax,9);
    func1->SetParNames("N482","Ek482","Sk482","Calib482","N976","Ek976","Sk976","calib976","Offset");
    func1->SetLineWidth(2);
    func1->SetLineColor(2);
    
    // N482 scaling
    func1->SetParameter(0,1);
    func1->SetParLimits(0,1,1e6);
  
    // N976 scaling
    func1->SetParameter(4,1);
    func1->SetParLimits(4,1,1e8);
  
    // Offset
    func1->SetParameter(9,1);
    func1->SetParLimits(9,0,500);
  
    // 482 keV Peak
    func1->SetParameter(1,Ek482);
    func1->SetParLimits(1,Ek482-(3.*s482),Ek482+(3.*s482));
    func1->SetParameter(2,s482);
    func1->SetParLimits(2,s482*0.30,s482*1.1);
    func1->SetParameter(3,calib);
    func1->SetParLimits(3,calib*0.10,calib*1.8);
  
    // 976 keV Peak
    func1->SetParameter(5,Ek976);
    func1->SetParLimits(5,Ek976-(3.*s976),Ek976+(3.*s976));
    func1->SetParameter(6,s976);
    func1->SetParLimits(6,s976*0.30,s976*2.0);
    func1->SetParameter(7,calib);
    func1->SetParLimits(7,calib*0.10,calib*1.8);
    
    ///////////////////////////////////////////////////////////////////
    TCanvas *c72 = new TCanvas("c72","- Rough Fit 1 -",xsize,ysize);
    xmin=280*calib;
    xmax=1100*calib;
    if(xmax>hmax) xmax=hmax;
    hbeta->Draw();
    hbeta->Fit("deltaE",fitopt,"",xmin,xmax);
    ///////////////////////////////////////////////////////////////////
    
    N500         = func1->GetParameter(0);
    N1000        = func1->GetParameter(4);
    
    Ek482        = func1->GetParameter(1);
    s482         = func1->GetParameter(2);
    Ek482_err    = func1->GetParError(1);
    s482_err     = func1->GetParError(2);
    cal482       = func1->GetParameter(3);
    cal482_err   = func1->GetParError(3);
    
    Ek976        = func1->GetParameter(5);
    s976         = func1->GetParameter(6);
    Ek976_err    = func1->GetParError(5);
    s976_err     = func1->GetParError(6);
    cal976       = func1->GetParameter(7);
    cal976_err   = func1->GetParError(7);
    
    calib        = (Ek976-Ek482)/494.;
    offset       = calib*((976.-Eloss)-(Ek976/calib));
    
    nEk976       = Ek976+offset;
    nEk482       = Ek482+offset;
    ns976        = TMath::Sqrt((s976*s976)-(sped*sped));
    ns482        = TMath::Sqrt((s482*s482)-(sped*sped));
    
    tot_chi2=func1->GetChisquare();
    tot_ndf=func1->GetNDF();
    
    pmin=xmin-(4*s482);
    pmax=xmax+(4*s976);
    if(pmin<hmin) pmin=hmin;
    if(pmax>hmax) pmax=hmax;
    hbeta->SetAxisRange(pmin,pmax,"X");
    c72->Print("A_bi207_rough1.png");
  
    cout<<"\n\n======================================================================";
    cout<<"\n ROUGH FIT 1 \n";
    cout<<"\t Pedestal      = "<<pdrun<<endl;
    cout<<"\t Beta File     = "<<betarun<<endl;
    //cout<<"\t Gamma File    = "<<gammarun<<endl;
    
    cout<<"\n Energy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    cout<<"\t Calib 976-482 = "<<calib<<endl;
    cout<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
    cout<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
    
    cout<<"\n 976keV peak fit \n";
    cout<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    cout<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err<<" = "<<s976/cal976<<" keV"<<endl;
    cout<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" % \n";
  
    cout<<"\n 482keV peak fit \n";
    cout<<"\t 482keV Mean   = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    cout<<"\t 482keV Sigma  = "<<s482<<" +/- "<<s482_err<<" = "<<s482/cal482<<" keV"<<endl;
    cout<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" % \n";
    
    outfile<<"======================================================================\n";
    outfile<<"\n ROUGH FIT 1 \n";
    outfile<<"\t Pedestal      = "<<pdrun<<endl;
    outfile<<"\t Beta File     = "<<betarun<<endl;
    //outfile<<"\t Gamma File    = "<<gammarun<<endl;
    
    outfile<<"\n Energy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    outfile<<"\t Calib 976-482 = "<<calib<<endl;
    outfile<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
    outfile<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
    
    outfile<<"\n 976keV peak fit \n";
    outfile<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    outfile<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err<<" = "<<s976/cal976<<" keV"<<endl;
    outfile<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" % \n";
  
    outfile<<"\n 482keV peak fit \n";
    outfile<<"\t 482keV Mean   = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    outfile<<"\t 482keV Sigma  = "<<s482<<" +/- "<<s482_err<<" = "<<s482/cal482<<" keV"<<endl;
    outfile<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" % \n";
    
    // fit results
    cout<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    cout<<"======================================================================\n\n\n";
    outfile<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    outfile<<"======================================================================\n\n\n";
    
    
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    //      2 ROUGH FIT THE BI207 SPECTRUM
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    
    // 482 keV Peak
    //func1->SetParameter(1,Ek482);
    //func1->SetParLimits(1,Ek482-(3.*s482),Ek482+(3.*s482));
    //func1->SetParameter(2,s482);
    //func1->SetParLimits(2,s482*0.10,s482*1.0);
    func1->FixParameter(3,calib);
      
    // 976 keV Peak
    //func1->SetParameter(5,Ek976);
    //func1->SetParLimits(5,Ek976-(3.*s976),Ek976+(3.*s976));
    //func1->SetParameter(6,s976);
    //func1->SetParLimits(6,s976*0.10,s976*3.0);
    func1->FixParameter(7,calib);
        
    ///////////////////////////////////////////////////////////////////
    TCanvas *c73 = new TCanvas("c73","- Rough Fit 2 -",xsize,ysize);
    xmin=280*calib;
    xmax=1100*calib;
    if(xmax>hmax) xmax=hmax;
    hbeta1->Draw();
    hbeta1->Fit("deltaE",fitopt,"",xmin,xmax);
    ///////////////////////////////////////////////////////////////////
  
    N500         = func1->GetParameter(0);
    N1000        = func1->GetParameter(4);
    
    Ek482        = func1->GetParameter(1);
    s482         = func1->GetParameter(2);
    Ek482_err    = func1->GetParError(1);
    s482_err     = func1->GetParError(2);
    cal482       = func1->GetParameter(3);
    cal482_err   = func1->GetParError(3);
    
    Ek976        = func1->GetParameter(5);
    s976         = func1->GetParameter(6);
    Ek976_err    = func1->GetParError(5);
    s976_err     = func1->GetParError(6);
    cal976       = func1->GetParameter(7);
    cal976_err   = func1->GetParError(7);
    
    calib        = (Ek976-Ek482)/494.;
    offset       = calib*((976.-Eloss)-(Ek976/calib));
    
    nEk976       = Ek976+offset;
    nEk482       = Ek482+offset;
    ns976        = TMath::Sqrt((s976*s976)-(sped*sped));
    ns482        = TMath::Sqrt((s482*s482)-(sped*sped));
    
    tot_chi2=func1->GetChisquare();
    tot_ndf=func1->GetNDF();
    
    pmin=xmin-(4*s482);
    pmax=xmax+(4*s976);
    if(pmin<hmin) pmin=hmin;
    if(pmax>hmax) pmax=hmax;
    hbeta1->SetAxisRange(pmin,pmax,"X");
    c73->Print("A_bi207_rough2.png");
    
    
    cout<<"\n\n======================================================================";
    cout<<"\n ROUGH FIT 2 \n";
    cout<<"\t Pedestal      = "<<pdrun<<endl;
    cout<<"\t Beta File     = "<<betarun<<endl;
    //cout<<"\t Gamma File    = "<<gammarun<<endl;
    
    cout<<"\n Energy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    cout<<"\t Calib 976-482 = "<<calib<<endl;
    cout<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
    cout<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
    
    cout<<"\n 976keV peak fit \n";
    cout<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    cout<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err<<" = "<<s976/cal976<<" keV"<<endl;
    cout<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" % \n";
  
    cout<<"\n 482keV peak fit \n";
    cout<<"\t 482keV Mean   = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    cout<<"\t 482keV Sigma  = "<<s482<<" +/- "<<s482_err<<" = "<<s482/cal482<<" keV"<<endl;
    cout<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" % \n";
    
    cout<<"\n 976keV ADC Offset and Pedestal Sigma Contribution \n";
    cout<<"\t Offset        = "<<offset<<" bins = "<<offset/cal976<<" keV"<<endl;
    cout<<"\t New 976 Mean  = "<<nEk976<<" bins = "<<nEk976/cal976<<" keV"<<endl;
    cout<<"\t New 976 Sigma = "<<ns976<<" bins = "<<ns976/cal976<<" keV"<<endl;
    cout<<"\t FWHM @ 976keV = "<<ns976/nEk976*2.354*100<<" +/- "<<s976_err/nEk976*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<ns976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" +/- "<<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" % \n";
    
    cout<<"\n 482keV ADC Offset and Pedestal Sigma Contribution \n";
    cout<<"\t Offset        = "<<offset<<" bins = "<<offset/cal482<<" keV"<<endl;
    cout<<"\t New 482 Mean  = "<<nEk482<<" bins = "<<nEk482/cal482<<" keV"<<endl;
    cout<<"\t New 482 Sigma = "<<ns482<<" bins = "<<ns482/cal482<<" keV"<<endl;
    cout<<"\t FWHM @ 482keV = "<<ns482/nEk482*2.354*100<<" +/- "<<s482_err/nEk482*2.354*100<<" % "<<endl;
    cout<<"\t FWHM @ 1MeV   = "<<ns482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" +/- "<<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" % \n";
    
    outfile<<"======================================================================\n";
    outfile<<"\n ROUGH FIT 2 \n";
    outfile<<"\t Pedestal      = "<<pdrun<<endl;
    outfile<<"\t Beta File     = "<<betarun<<endl;
    //outfile<<"\t Gamma File    = "<<gammarun<<endl;
    
    outfile<<"\n Energy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    outfile<<"\t Calib 976-482 = "<<calib<<endl;
    outfile<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
    outfile<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
    
    outfile<<"\n 976keV peak fit \n";
    outfile<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    outfile<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err<<" = "<<s976/cal976<<" keV"<<endl;
    outfile<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" % \n";
  
    outfile<<"\n 482keV peak fit \n";
    outfile<<"\t 482keV Mean   = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    outfile<<"\t 482keV Sigma  = "<<s482<<" +/- "<<s482_err<<" = "<<s482/cal482<<" keV"<<endl;
    outfile<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" % \n";
    
    outfile<<"\n 976keV ADC Offset and Pedestal Sigma Contribution \n";
    outfile<<"\t Offset        = "<<offset<<" bins = "<<offset/cal976<<" keV"<<endl;
    outfile<<"\t New 976 Mean  = "<<nEk976<<" bins = "<<nEk976/cal976<<" keV"<<endl;
    outfile<<"\t New 976 Sigma = "<<ns976<<" bins = "<<ns976/cal976<<" keV"<<endl;
    outfile<<"\t FWHM @ 976keV = "<<ns976/nEk976*2.354*100<<" +/- "<<s976_err/nEk976*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<ns976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" +/- "<<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" % \n";
    
    outfile<<"\n 482keV ADC Offset and Pedestal Sigma Contribution \n";
    outfile<<"\t Offset        = "<<offset<<" bins = "<<offset/cal482<<" keV"<<endl;
    outfile<<"\t New 482 Mean  = "<<nEk482<<" bins = "<<nEk482/cal482<<" keV"<<endl;
    outfile<<"\t New 482 Sigma = "<<ns482<<" bins = "<<ns482/cal482<<" keV"<<endl;
    outfile<<"\t FWHM @ 482keV = "<<ns482/nEk482*2.354*100<<" +/- "<<s482_err/nEk482*2.354*100<<" % "<<endl;
    outfile<<"\t FWHM @ 1MeV   = "<<ns482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" +/- "<<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" % \n";
    
    // fit results
    cout<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    cout<<"======================================================================\n\n\n";
    outfile<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    outfile<<"======================================================================\n\n\n";
    
    
    
    /**
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    //      FINAL FIT THE BI207 SPECTRUM
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    
    func1->ReleaseParameter(4);
    func1->SetParLimits(4,1,100);
    func1->ReleaseParameter(5);
    func1->SetParLimits(5,0,0.02);
  
    func1->FixParameter(10,Nratio);
  
    func1->FixParameter(16,calib);
    func1->FixParameter(20,calib);
  
    func1->SetParameter(14,Ek482);
    func1->SetParLimits(14,Ek482*0.9,Ek482*1.1);
    func1->SetParameter(15,s482);
    func1->SetParLimits(15,s482*0.9,s482*1.1);
    
    func1->SetParameter(18,Ek976);
    func1->SetParLimits(18,Ek976*0.9,Ek976*1.1);
    func1->SetParameter(19,s976);
    func1->SetParLimits(19,s976*0.9,s976*1.1);
  
    ///////////////////////////////////////////////////////////////////
    TCanvas *c74 = new TCanvas("c74","- Final Fit -",xsize,ysize);
    xmin=280*calib;
    xmax=1200*calib;
    //xmax=Ek976+(10*s976);
    if(xmax>hmax) xmax=hmax;
    hbeta2->Draw();
    hbeta2->Fit("fall","Q","",xmin,xmax);
    hbeta2->Fit("fall","M","",xmin,xmax);
    ///////////////////////////////////////////////////////////////////
  
    cedge500     = func1->GetParameter(3);
    cedge1000    = func1->GetParameter(7);
    N500         = func1->GetParameter(2);
    N1000        = func1->GetParameter(6);
    csigma500    = func1->GetParameter(4);
    csigma1000   = func1->GetParameter(8);
    c2sigma500   = func1->GetParameter(5);
    c2sigma1000  = func1->GetParameter(9);
    Ngaus        = func1->GetParameter(10);
    gausmean     = func1->GetParameter(11);
    gaussig      = func1->GetParameter(12);
    Nratio       = func1->GetParameter(10);
  
    Ek976        = func1->GetParameter(18);
    s976         = func1->GetParameter(19);
    Ek976_err    = func1->GetParError(18);
    s976_err     = func1->GetParError(19);
    cal976       = func1->GetParameter(20);
    cal976_err   = func1->GetParError(20);
  
    Ek482        = func1->GetParameter(14);
    s482         = func1->GetParameter(15);
    Ek482_err    = func1->GetParError(14);
    s482_err     = func1->GetParError(15);
    cal482       = func1->GetParameter(16);
    cal482_err   = func1->GetParError(16);
  
    calib        = (Ek976-Ek482)/494.;
    offset       = TMath::Abs(((Ek976+Ek482)/(2.*cal976))-729.);
    nEk976       = Ek976+offset;
    nEk482       = Ek482+offset;
    
    tot_chi2=func1->GetChisquare();
    tot_ndf=func1->GetNDF();
    
    pmin=xmin-(3*csigma500);
    pmax=xmax+(4*s976);
    if(pmin<hmin) pmin=hmin;
    if(pmax>hmax) pmax=hmax;
    hbeta2->SetAxisRange(pmin,pmax,"X");
    c74->Print("A_bi207_final.png");
  
    cout<<"\n\n======================================================================\n";
    cout<<"FINAL FIT \n";
    cout<<"    Pedistal File    = "<<pdrun<<endl;
    cout<<"    Beta+Gamma File  = "<<allrun<<endl;
    cout<<"    Gamma File       = "<<gammarun<<endl;
    
    cout<<"\nEnergy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    cout<<"    Calib from Ek976-Ek482 = "<<calib<<endl;
    cout<<"    Calib from par[16] = "<<cal482<<endl;
    cout<<"    Calib from par[20] = "<<cal976<<endl;
  
    cout<<"\n976keV peak fit \n";
    cout<<"    Mean = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    cout<<"    Sigma = "<<s976<<" +/- "<<s976_err<<endl;
    cout<<"    FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    cout<<"    FWHM @ 1MeV = "<< s976/Ek976*2.354*100*TMath::Sqrt(0.976)<<" +/- "<< s976_err/Ek976*2.354*100*TMath::Sqrt(0.976) <<" /Sqrt(E) % (tabled energy) \n";
    cout<<"    FWHM @ 1MeV = "<< s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" /Sqrt(E) % (measured energy) \n";
  
    cout<<"\n482keV peak fit \n";
    cout<<"    Mean = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    cout<<"    Sigma = "<<s482<<" +/- "<<s482_err<<endl;
    cout<<"    FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    cout<<"    FWHM @ 1MeV = "<< s482/Ek482*2.354*100*TMath::Sqrt(0.482)<<" +/- "<< s482_err/Ek482*2.354*100*TMath::Sqrt(0.482) <<" /Sqrt(E) % (tabled energy) \n";
    cout<<"    FWHM @ 1MeV = "<< s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    cout<<"\n976keV Offset Results \n";
    cout<<"    Offset+Eloss = "<<offset<<" bins = "<<offset/cal976<<" keV"<<endl;
    cout<<"    New 976 Mean = "<<nEk976<<" bins = "<<nEk976/cal976<<" keV"<<endl;
    cout<<"    FWHM @ 976keV = "<<s976/nEk976*2.354*100<<" +/- "<<s976_err/nEk976*2.354*100<<" % "<<endl;
    cout<<"    FWHM @ 1MeV = "<< s976/nEk976*2.354*100*TMath::Sqrt(0.976)<<" +/- "<< s976_err/nEk976*2.354*100*TMath::Sqrt(0.976) <<" /Sqrt(E) % (tabled energy) \n";
    cout<<"    FWHM @ 1MeV = "<< s976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" +/- "<<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    cout<<"\n482keV Offset Results \n";
    cout<<"    Offset+Eloss = "<<offset<<" bins = "<<offset/cal482<<" keV"<<endl;
    cout<<"    New 482 Mean = "<<nEk482<<" bins = "<<nEk482/cal482<<" keV"<<endl;
    cout<<"    FWHM @ 482keV = "<<s482/nEk482*2.354*100<<" +/- "<<s482_err/nEk482*2.354*100<<" % "<<endl;
    cout<<"    FWHM @ 1MeV = "<< s482/nEk482*2.354*100*TMath::Sqrt(0.482)<<" +/- "<< s482_err/nEk482*2.354*100*TMath::Sqrt(0.482) <<" /Sqrt(E) % (tabled energy) \n";
    cout<<"    FWHM @ 1MeV = "<< s482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" +/- "<<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    outfile<<"======================================================================\n";
    outfile<<"FINAL FIT \n";
    outfile<<"    Pedistal File    = "<<pdrun<<endl;
    outfile<<"    Beta+Gamma File  = "<<allrun<<endl;
    outfile<<"    Gamma File       = "<<gammarun<<endl;
    
    outfile<<"\nEnergy Calibration from Ek976-Ek482 = "<<1/calib<<" keV/ch"<<endl;
    outfile<<"    Calib from Ek976-Ek482 = "<<calib<<endl;
    outfile<<"    Calib from par[16] = "<<cal482<<endl;
    outfile<<"    Calib from par[20] = "<<cal976<<endl;
  
    outfile<<"\n976keV peak fit \n";
    outfile<<"    Mean = "<<Ek976<<" +/- "<<Ek976_err<<" = "<<Ek976/cal976<<" keV"<<endl;
    outfile<<"    Sigma = "<<s976<<" +/- "<<s976_err<<endl;
    outfile<<"    FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "<<s976_err/Ek976*2.354*100<<" % "<<endl;
    outfile<<"    FWHM @ 1MeV = "<< s976/Ek976*2.354*100*TMath::Sqrt(0.976)<<" +/- "<< s976_err/Ek976*2.354*100*TMath::Sqrt(0.976) <<" /Sqrt(E) % (tabled energy) \n";
    outfile<<"    FWHM @ 1MeV = "<< s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" +/- "<<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)<<" /Sqrt(E) % (measured energy) \n";
  
    outfile<<"\n482keV peak fit \n";
    outfile<<"    Mean = "<<Ek482<<" +/- "<<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
    outfile<<"    Sigma = "<<s482<<" +/- "<<s482_err<<endl;
    outfile<<"    FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
    outfile<<"    FWHM @ 1MeV = "<< s482/Ek482*2.354*100*TMath::Sqrt(0.482)<<" +/- "<< s482_err/Ek482*2.354*100*TMath::Sqrt(0.482) <<" /Sqrt(E) % (tabled energy) \n";
    outfile<<"    FWHM @ 1MeV = "<< s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" +/- "<<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    outfile<<"\n976keV Offset Results \n";
    outfile<<"    Offset+Eloss = "<<offset<<" bins = "<<offset/cal976<<" keV"<<endl;
    outfile<<"    New 976 Mean = "<<nEk976<<" bins = "<<nEk976/cal976<<" keV"<<endl;
    outfile<<"    FWHM @ 976keV = "<<s976/nEk976*2.354*100<<" +/- "<<s976_err/nEk976*2.354*100<<" % "<<endl;
    outfile<<"    FWHM @ 1MeV = "<< s976/nEk976*2.354*100*TMath::Sqrt(0.976)<<" +/- "<< s976_err/nEk976*2.354*100*TMath::Sqrt(0.976) <<" /Sqrt(E) % (tabled energy) \n";
    outfile<<"    FWHM @ 1MeV = "<< s976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" +/- "<<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    outfile<<"\n482keV Offset Results \n";
    outfile<<"    Offset+Eloss = "<<offset<<" bins = "<<offset/cal482<<" keV"<<endl;
    outfile<<"    New 482 Mean = "<<nEk482<<" bins = "<<nEk482/cal482<<" keV"<<endl;
    outfile<<"    FWHM @ 482keV = "<<s482/nEk482*2.354*100<<" +/- "<<s482_err/nEk482*2.354*100<<" % "<<endl;
    outfile<<"    FWHM @ 1MeV = "<< s482/nEk482*2.354*100*TMath::Sqrt(0.482)<<" +/- "<< s482_err/nEk482*2.354*100*TMath::Sqrt(0.482) <<" /Sqrt(E) % (tabled energy) \n";
    outfile<<"    FWHM @ 1MeV = "<< s482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" +/- "<<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)<<" /Sqrt(E) % (measured energy) \n";
    
    // fit results
    cout<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    cout<<"======================================================================\n\n\n";
    outfile<<"\nFIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "<<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
    outfile<<"======================================================================\n\n\n";
    
    **/
    
    outfile.close();
  }
  else cout<<"\n\n----- COULD NOT OPEN OUTPUT FILE -----\n\n";
}

