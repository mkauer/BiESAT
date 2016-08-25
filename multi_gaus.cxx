/******************************************************************
A little program from Zori to fit a low Photo Electron spectra.

version: 4
build:   08.10.31

change log: (most recent at top)
+ added errors for mean and Npe
+ added correction factor to Npe calculation
~ flipped axis around
~ use the gauss mean instead of hist mean for adc axis
+ added full plot canvas to show all fits
~ changed this to run on gaussian histos
+ lots of minor fixes and tweaks
+ added while loop to go through all runs
+ file output to store all parameters
+ changed the file input format

Matt Kauer
******************************************************************/

#include "_userFuncs.hxx"
#include "_readinData.hxx"
#include "_fittingFuncs.hxx"

void multi_gaus(char *runlist="runlist-gaus5.dat", Int_t xmax=3800, Int_t xmin=0){ 
  
  beautify();  
  
  const Int_t daq=0;    // 0=vme : 1=camac
  const Int_t chan=0;   // 0-7==lowRes : 8-15==hiRes
  const Int_t rebin=1;
  
  const Int_t xcan=300;
  const Int_t ycan=200;
  
  //const Float_t corr_factor=1.0;  // runlist-gaus4.dat
  const Float_t corr_factor=1.2;  // runlist-gaus5.dat
  
  Bool_t correct=true;
  const Int_t npoints=46;
  const Int_t fit_points=11;
  Bool_t fixpar0=false;
  
  //char *pedfile="run24192"; // runlist-gaus1.dat
  //char *pedfile="run24205"; // runlist-gaus2.dat
  //char *pedfile="run24221"; // runlist-gaus3.dat
  //char *pedfile="run24314"; // runlist-gaus4.dat
  char *pedfile="run24420"; // runlist-gaus5.dat
  
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
  if(pmax>3800) pmax=3800;
  
  TF1 *pgauss = new TF1("pgauss","gaus",xmin,xmax);
  pgauss->SetLineColor(2);
  pgauss->SetLineWidth(2);
  
  TCanvas *c200 = new TCanvas("c200","c200",xcan,ycan);
  pedhist->Draw();
  pedhist->Fit("pgauss","Q","",pmin,pmax);
  pedhist->SetAxisRange(pmin,pmax,"X");
  c200->Update();
  
  Float_t ped = pgauss->GetParameter(1);
  Float_t sped = pgauss->GetParameter(2);
  Float_t ped_err = pgauss->GetParError(1);
  Float_t sped_err = pgauss->GetParError(2);
  cout<<"\n\t Pedestal = "<<ped<<"  +/-  "<<sped<<"\n\n";
  
  
  ////////////////////////////////////////////////////
  //  OPEN THE RUNLIST AND WRITE OUT RESULTS
  ////////////////////////////////////////////////////
  
  FILE *outfile=fopen("multi-gaus-results.txt","w");
  //fprintf(outfile,"NDF\t hmean\t hrms\t gmean\t gsig\t gmerr\t gserr \n");
  fprintf(outfile,"Filename    ADC    Npe  errADC errNpe \n");
  fprintf(outfile,"========= ====== ====== ====== ====== \n");
  
  const Int_t nlines=getLines(runlist);
  ifstream infile;
  infile.open(runlist);
  
  //const Int_t npoints=nlines;
  Double_t hist_mean[npoints];
  Double_t npe_mean_rms[npoints];
  Double_t gmean_ped[npoints];
  Double_t npe_mean_rms_err[npoints];
  Double_t gmean_ped_err[npoints];
  
  Int_t k=0;
  
  Float_t maxvalue=0;
  Float_t NDfilter=0;
  char *runfile="\n";
  
  Int_t col=5;
  Int_t rows=TMath::Ceil(nlines/Double_t(col));
  TCanvas *c555 = new TCanvas("c555","c555",xcan*col,ycan*rows);
  c555->Divide(col,rows,0.002,0.002,0);
  while(infile>>runfile){
    //cout<<"\n\n\t"<<runfile[0]<<" "<<runfile[1]<<"\n\n";
    if(runfile[0]!='#'){
      //infile>>NDfilter;
      //cout<<"\n\t Analyzing --> "<<runfile<<" with ND filter = "<<NDfilter<<"\n\n";
      cout<<"\n\t Analyzing --> "<<runfile<<"\n\n";
      //fprintf(outfile,"%.1f\t ",NDfilter);
      fprintf(outfile,"%-10s",runfile);
      
      TH1F *histio=new TH1F("histio",runfile,4000,0,4000);
      readfile(runfile,0,chan,histio,daq);
      histio->Rebin(rebin);
      hist_mean[k]=histio->GetMean()-ped;
      Float_t hist_rms=histio->GetRMS();
      //fprintf(outfile,"%.3f\t %.3f\t ",hist_mean[k],hist_rms);
      
      TF1 *hgaus = new TF1("hgaus","gaus",xmin,xmax);
      hgaus->SetLineColor(2);
      hgaus->SetLineWidth(2);
      hgaus->SetParameter(1,hist_mean[k]);
      hgaus->SetParameter(2,hist_rms);
      
      pmin=ped+hist_mean[k]-(3*TMath::Sqrt(hist_mean[k]));
      if(pmin<=ped+sped)pmin=ped+sped;
      pmax=ped+hist_mean[k]+(4*TMath::Sqrt(hist_mean[k]));
      if(pmax>3800) pmax=3800;
      
      //TCanvas *c201 = new TCanvas("c201","c201",xcan,ycan);
      c555->cd(Int_t(k+1));
      histio->Draw();
      histio->Fit("hgaus","","",pmin,pmax);
      
      Float_t gmean = hgaus->GetParameter(1);
      Float_t gsig = hgaus->GetParameter(2);
      Float_t gmean_err = hgaus->GetParError(1);
      Float_t gsig_err = hgaus->GetParError(2);
      
      pmin=gmean-(5*gsig);
      if(pmin<=ped+sped)pmin=ped+sped;
      pmax=gmean+(7*gsig);
      if(pmax>3800) pmax=3800;
      histio->Fit("hgaus","","",pmin,pmax);
      
      Float_t gmean = hgaus->GetParameter(1);
      Float_t gsig = hgaus->GetParameter(2);
      Float_t gmean_err = hgaus->GetParError(1);
      Float_t gsig_err = hgaus->GetParError(2);
      
      histio->SetAxisRange((gmean-(5*gsig)),(gmean+(7*gsig)),"X");
      
      //fprintf(outfile,"%.3f\t %.3f\t %.3f\t %.3f ",gmean,gsig,gmean_err,gsig_err);
      
      gmean_ped[k]=gmean-ped;
      gmean_ped_err[k]=TMath::Sqrt(gmean_err*gmean_err+ped_err*ped_err);
      npe_mean_rms[k]=TMath::Power(gmean-ped,2)/((gsig*gsig)-(sped*sped));
      //npe_mean_rms[k]=TMath::Power((gmean-ped)/gsig,2);
      if(correct) npe_mean_rms[k]=npe_mean_rms[k]*corr_factor;
      if(npe_mean_rms[k]>maxvalue) maxvalue=npe_mean_rms[k];
      npe_mean_rms_err[k]=npe_mean_rms[k]*0.01;  //--- Error on correction factor
      
      fprintf(outfile,"%6i %6.3f %6i %6.2f",gmean_ped[k],gmean_ped_err[k],npe_mean_rms[k],npe_mean_rms_err[k]);
      fprintf(outfile,"\n");
      
      k++;
    } // end of if statement
    else cout<<"\n\t "<<runfile<<" was skipped!\n";
  } // end of while loop
  
  c555->Print("multi-gaus-plots.png");
  fclose(outfile);
  infile.close();
  
  Int_t newmax=(TMath::CeilNint(maxvalue/1000)+1)*1000;
  ///////////////////////////////////////////////////
  //  PLOT MEAN/SIGMA VS ADC-MEAN
  ///////////////////////////////////////////////////
  TCanvas *c_v10=new TCanvas("c_v10","Results",700,500);
  c_v10->SetFillColor(0);
  c_v10->SetGrid();
  c_v10->GetFrame()->SetFillColor(0);
  c_v10->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v10=c_v10->DrawFrame(0,0,newmax,4000);
  mg_v10->SetXTitle("Npe");
  mg_v10->SetYTitle("ADC Mean");
  //TH1F *mg_v10=c_v10->DrawFrame(0,0,4000,newmax);
  //mg_v10->SetXTitle("ADC Mean");
  //mg_v10->SetYTitle("Npe");
  mg_v10->SetTitle("(Mean/Sigma)^2");
  if(correct) mg_v10->SetTitle("(Mean/Sigma)^2 * Correction");
  
  TGraphErrors *gain_v10=new TGraphErrors(npoints,npe_mean_rms,gmean_ped,npe_mean_rms_err,gmean_ped_err);
  //TGraph *gain_v10=new TGraph(npoints,npe_mean_rms,gmean_ped);
  //TGraph *gain_v10=new TGraph(npoints,gmean_ped,npe_mean_rms);
  gain_v10->SetMarkerColor(kBlue);
  gain_v10->SetMarkerStyle(20);
  gain_v10->SetMarkerSize(1);
  gain_v10->Draw("P");
  
  TGraphErrors *fitpts=new TGraphErrors(fit_points,npe_mean_rms,gmean_ped,npe_mean_rms_err,gmean_ped_err);
  //TGraph *fitpts=new TGraph(fit_points,npe_mean_rms,gmean_ped);
  //TGraph *fitpts=new TGraph(fit_points,gmean_ped,npe_mean_rms);
  fitpts->SetMarkerColor(kBlack);
  fitpts->SetMarkerStyle(20);
  fitpts->SetMarkerSize(1);
  fitpts->Draw("P");
  
  TF1 *func_v10=new TF1("func_v10",f_line,0,newmax,2);
  //TF1 *func_v10=new TF1("func_v10",f_line,0,4000,2);
  func_v10->SetLineColor(2);
  func_v10->SetLineWidth(2);
  if(fixpar0) func_v10->FixParameter(0,0);
  
  fitpts->Fit("func_v10","Q","",0,newmax);
  c_v10->Print("multi-gaus-linearity.png");
  
}


//--- basic pol2
Double_t f_pol2(Double_t *xval, Double_t *par)
{
  return par[0]+(par[1]*xval[0])+(par[2]*xval[0]*xval[0]);
}

