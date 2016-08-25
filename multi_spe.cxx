/******************************************************************
A little program from Zori to fit a low Photo Electron spectra.

version: 4
build:   08.10.31

change log: (most recent at top)
+ lots of minor fixes and tweaks
+ added while loop to go through all runs
+ file output to store all parameters
+ changed the file input format

Matt Kauer
******************************************************************/

#include "_fittingFuncs.hxx"
#include "_userFuncs.hxx"
#include "_readinData.hxx"

void multi_spe(char *runfile="\n", Int_t xmax=4000, Int_t xmin=0){ 
  
  beautify();  
  
  const Int_t daq=0;    // 0=vme : 1=camac
  const Int_t chan=8;   // 0-7==lowRes : 8-15==hiRes
  const Int_t rebin=2;
  
  const Int_t xcan=400;
  const Int_t ycan=300;
  
  char *runlist="runlist-spe.dat";
  char *pedfile="run24141";
  
  ////////////////////////////////////////////////////
  //  FIT THE PEDESTAL
  ////////////////////////////////////////////////////
  TH1F *pedhist = new TH1F("pedhist",pedfile,4000,0,4000);
  readfile(pedfile,0,chan,pedhist,daq);
  pedhist->Rebin(rebin);
  
  Float_t ped_mean=pedhist->GetMean();
  Float_t ped_rms=pedhist->GetRMS();
  Int_t pmin=ped_mean-6*ped_rms;
  if(pmin<10)pmin=0;
  Int_t pmax=ped_mean+8*ped_rms;
  if(pmax>3990) pmax=4000;
  
  TF1 *pgauss = new TF1("pgauss","gaus",xmin,xmax);
  pgauss->SetLineColor(2);
  pgauss->SetLineWidth(2);
  
  TCanvas *c200 = new TCanvas("c200","c200",xcan,ycan);
  pedhist->Draw();
  pedhist->Fit("pgauss","Q","",pmin,pmax);
  pedhist->SetAxisRange(pmin,pmax,"X");
  
  Float_t ped = pgauss->GetParameter(1);
  Float_t sped = pgauss->GetParameter(2);
  Float_t ped_err = pgauss->GetParError(1);
  Float_t sped_err = pgauss->GetParError(2);
  cout<<"\n\t Pedestal = "<<ped<<"  +/-  "<<sped<<"\n\n";
  
  
  ////////////////////////////////////////////////////
  //  OPEN THE RUNLIST AND GET IT READY
  ////////////////////////////////////////////////////
  Float_t prev_spe=1.1;
  Float_t prev_pois=1.1;
  
  FILE *outfile=fopen("multi_spe_results.txt","w");
  
  ifstream infile;
  infile.open(runlist);
  Float_t NDfilter=0;
  while(infile>>runfile){
    infile>>NDfilter;
    cout<<"\n\t Analyzing --> "<<runfile<<" with ND filter = "<<NDfilter<<"\n\n";
    
    fprintf(outfile,"%.1f  ",NDfilter);
    
    TH1F *hist=new TH1F("hist",runfile,4000,0,4000);
    readfile(runfile,0,chan,hist,daq);
    hist->Rebin(rebin);
    Float_t hist_mean=hist->GetMean();
    Float_t hist_rms=hist->GetRMS();
    fprintf(outfile,"%.3f  %.3f  ",hist_mean,hist_rms);
    
    
    ////////////////////////////////////////////////////
    //  FIT WITH POISSON DISTRIBUTION
    ////////////////////////////////////////////////////
    TH1F *poiss=new TH1F("Poisson Fit",runfile,4000,0,4000);
    readfile(runfile,TMath::CeilNint(ped),chan,poiss,daq);
    poiss->Rebin(rebin);
  
    Float_t dist_mean=poiss->GetMean();
    Float_t dist_rms=poiss->GetRMS();
    pmin=dist_mean-4*dist_rms;
    if(pmin<10)pmin=0;
    pmax=dist_mean+7*dist_rms;
    if(pmax>3830-Int_t(ped)) pmax=3830-Int_t(ped);
    
    TF1 *func = new TF1("myfit",f_poisson,xmin,xmax,5);
    func->SetParNames("pois_norm","pois_mean","pois_factor","gaus_norm","gaus_sig");
    func->SetLineColor(2);
    func->SetLineWidth(2);
    
    func->SetParLimits(0, 1.0, 100000);
    func->SetParLimits(1, 0.1, 1000);
    func->SetParLimits(2, 0.1, 500);
    func->SetParLimits(3, 0.0, 1000);
    func->SetParLimits(4, 0.0001, 20);
    
    func->SetParameter(0,1000);
    func->SetParameter(1,prev_pois);
    func->SetParameter(2,30);
    func->SetParameter(3,10);
    func->FixParameter(4,sped);
    
    TCanvas *c202 = new TCanvas("c202","c202",xcan,ycan);
    poiss->Draw();
    poiss->Fit("myfit","WW","",pmin,pmax);
    poiss->SetAxisRange(pmin,pmax,"X");
    c202->Print("AA_pois-fit.png");
    
    prev_pois=func->GetParameter(1);
    
    Float_t pois_chi2=func->GetChisquare();
    Float_t pois_ndf=func->GetNDF();
    Float_t pois_mean=func->GetParameter(1);
    Float_t pois_mean_err=func->GetParError(1);
    
    fprintf(outfile,"%.3f  %.3f  ",pois_mean,pois_mean_err);
    
    
    fprintf(outfile,"\n");
  } // end of while loop
  
  
  fclose(outfile);
  infile.close();
  cout<<"\n\n DONE DONE DONE \n\n";
  
}
