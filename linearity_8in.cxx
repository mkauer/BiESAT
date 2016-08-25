#include "_userFuncs.hxx"
#include "_fittingFuncs.hxx"

void linearity_8in(char *fname="multi_spe_results.dat"){
  
  beautify();
  
  Double_t pe_mean=25.74;
  Double_t correction=1.099;
  Double_t ped=323.6;
  Double_t ped_rms=10.55;
  
  const Int_t ncols=getCol(fname);
  const Int_t nlines=getLines(fname);
  cout<<"\n Columns --> "<<ncols<<"";
  cout<<"\n Lines   --> "<<nlines<<"\n\n";
  
  Double_t arr2d[nlines][ncols];
  ifstream infile;
  infile.open(fname);
  Int_t i=0,j=0;
  Double_t n;
  while(infile>>n){
    arr2d[i][j]=n;
    j++;
    if(j==ncols){
      i++;
      j=0;
    }
    //cout<<n<<"  "<<i<<" "<<j<<endl;
  }
  
  const Int_t npoints=nlines;
  Double_t ndf[npoints];
  Double_t ndferr[npoints];
  Double_t trans[npoints];
  Double_t transerr[npoints];
  
  Double_t hist_mean[npoints];
  Double_t hist_meanped[npoints];
  
  Double_t npe_mean_rms[npoints];
  Double_t nperr_mean_rms[npoints];
  
  Double_t npe_pois[npoints];
  Double_t nperr_pois[npoints];
  
  Double_t npe_mean_1pe[npoints];
  Double_t nperr_mean_1pe[npoints];
  
  
  for(Int_t k=0;k<npoints;k++){
    
    ndf[k]=arr2d[k][0];
    ndferr[k]=ndf[k]*0.10;
    
    trans[k]=1./TMath::Power(10,ndf[k]);
    transerr[k]=trans[k]*0.10;
    
    hist_mean[k]=arr2d[k][1];
    hist_meanped[k]=arr2d[k][1]-ped;
    
    npe_mean_rms[k]=((TMath::Power(arr2d[k][1]-ped,2))/((arr2d[k][2]*arr2d[k][2])-(ped_rms*ped_rms)))*correction;
    nperr_mean_rms[k]=npe_mean_rms[k]*0.05;
    
    npe_pois[k]=arr2d[k][3]*correction;
    nperr_pois[k]=npe_pois[k]*0.05;
    
    npe_mean_1pe[k]=(arr2d[k][1]-ped)/pe_mean;
    nperr_mean_1pe[k]=npe_mean_1pe[k]*0.05;
    
  }
  
  Int_t ymax=200;
  Double_t xmax=1.2;
  
  //----- mean-ped/rms-ped_rms * corr vs NDfilter
  ///////////////////////////////////////////////////
  TCanvas *c_v0=new TCanvas("c_v0","Results",700,500);
  c_v0->SetFillColor(0);
  c_v0->SetGrid();
  c_v0->GetFrame()->SetFillColor(0);
  c_v0->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v0=c_v0->DrawFrame(0,0,2.4,100);
  mg_v0->SetXTitle("ND Filter");
  mg_v0->SetYTitle("Npe");
  mg_v0->SetTitle("Mean-ped/RMS-pedRMS * corr");
    
  //TGraph *gain_v0=new TGraph(npoints,ndf,npe_mean_rms);
  TGraphErrors *gain_v0=new TGraphErrors(npoints,ndf,npe_mean_rms,ndferr,nperr_mean_rms);
  gain_v0->SetMarkerColor(kBlue);
  gain_v0->SetMarkerStyle(20);
  gain_v0->SetMarkerSize(1);
  gain_v0->Draw("P");
  
  TF1 *func_v0=new TF1("func_v0",f_expo,0,2.4,3);
  func_v0->SetLineColor(2);
  func_v0->SetLineWidth(2);
  
  gain_v0->Fit("func_v0","WR");
  c_v0->Print("8in_expo.png");
  
  
  //----- mean-ped/rms-ped_rms * corr vs transmission
  ///////////////////////////////////////////////////
  TCanvas *c_v1=new TCanvas("c_v1","Results",700,500);
  c_v1->SetFillColor(0);
  c_v1->SetGrid();
  c_v1->GetFrame()->SetFillColor(0);
  c_v1->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v1=c_v1->DrawFrame(0,0,xmax,ymax);
  mg_v1->SetXTitle("Transmission");
  mg_v1->SetYTitle("Npe");
  mg_v1->SetTitle("Mean-ped/RMS-pedRMS * corr");
    
  //TGraph *gain_v1=new TGraph(npoints,trans,npe_mean_rms);
  TGraphErrors *gain_v1=new TGraphErrors(npoints,trans,npe_mean_rms,transerr,nperr_mean_rms);
  gain_v1->SetMarkerColor(kBlue);
  gain_v1->SetMarkerStyle(20);
  gain_v1->SetMarkerSize(1);
  gain_v1->Draw("P");
  
  TF1 *func_v1=new TF1("func_v1",f_line,0,xmax,2);
  func_v1->SetLineColor(2);
  func_v1->SetLineWidth(2);
  func_v1->FixParameter(0,0);
  
  gain_v1->Fit("func_v1","R");
  c_v1->Print("8in_mean_rms.png");
  
  
  //----- pois_mean * corr vs transmission
  ///////////////////////////////////////////////////
  TCanvas *c_v2=new TCanvas("c_v2","Results",700,500);
  c_v2->SetFillColor(0);
  c_v2->SetGrid();
  c_v2->GetFrame()->SetFillColor(0);
  c_v2->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v2=c_v2->DrawFrame(0,0,xmax,ymax);
  mg_v2->SetXTitle("Transmission");
  mg_v2->SetYTitle("Npe");
  mg_v2->SetTitle("Poisson_Mean * corr");
  
  //TGraph *gain_v2=new TGraph(npoints,trans,npe_pois);
  TGraphErrors *gain_v2=new TGraphErrors(npoints,trans,npe_pois,transerr,nperr_pois);
  gain_v2->SetMarkerColor(kBlue);
  gain_v2->SetMarkerStyle(20);
  gain_v2->SetMarkerSize(1);
  gain_v2->Draw("P");
  
  TF1 *func_v2=new TF1("func_v2",f_line,0,xmax,2);
  func_v2->SetLineColor(2);
  func_v2->SetLineWidth(2);
  func_v2->FixParameter(0,0);
  
  gain_v2->Fit("func_v2","R");
  c_v2->Print("8in_pois.png");
  
  
  //----- mean-ped/1pe vs transmission
  ///////////////////////////////////////////////////
  TCanvas *c_v3=new TCanvas("c_v3","Results",700,500);
  c_v3->SetFillColor(0);
  c_v3->SetGrid();
  c_v3->GetFrame()->SetFillColor(0);
  c_v3->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v3=c_v3->DrawFrame(0,0,xmax,ymax);
  mg_v3->SetXTitle("Transmission");
  mg_v3->SetYTitle("Npe");
  mg_v3->SetTitle("Mean-ped/1pe_mean");
  
  //TGraph *gain_v3=new TGraph(npoints,trans,npe_mean_1pe);
  TGraphErrors *gain_v3=new TGraphErrors(npoints,trans,npe_mean_1pe,transerr,nperr_mean_1pe);
  gain_v3->SetMarkerColor(kBlue);
  gain_v3->SetMarkerStyle(20);
  gain_v3->SetMarkerSize(1);
  gain_v3->Draw("P");
  
  TF1 *func_v3=new TF1("func_v3",f_line,0,xmax,2);
  func_v3->SetLineColor(2);
  func_v3->SetLineWidth(2);
  func_v3->FixParameter(0,0);
  
  gain_v3->Fit("func_v3","R");
  c_v3->Print("8in_mean_1pe.png");
  
  
  //----- pois_mean * corr vs adc_mean
  ///////////////////////////////////////////////////
  TCanvas *c_v5=new TCanvas("c_v5","Results",700,500);
  c_v5->SetFillColor(0);
  c_v5->SetGrid();
  c_v5->GetFrame()->SetFillColor(0);
  c_v5->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v5=c_v5->DrawFrame(0,0,4000,ymax);
  mg_v5->SetXTitle("ADC Mean");
  mg_v5->SetYTitle("Npe");
  mg_v5->SetTitle("Poisson_Mean * corr");
  
  TGraph *gain_v5=new TGraph(npoints,hist_meanped,npe_pois);
  gain_v5->SetMarkerColor(kBlue);
  gain_v5->SetMarkerStyle(20);
  gain_v5->SetMarkerSize(1);
  gain_v5->Draw("P");
  
  TF1 *func_v5=new TF1("func_v5",f_line,0,4000,2);
  func_v5->SetLineColor(2);
  func_v5->SetLineWidth(2);
  func_v5->FixParameter(0,0);
  
  gain_v5->Fit("func_v5","WR");
  c_v5->Print("8in_pois_adc.png");
  
  
  //----- mean-ped/rms-ped_rms * corr vs adc_mean
  ///////////////////////////////////////////////////
  TCanvas *c_v6=new TCanvas("c_v6","Results",700,500);
  c_v6->SetFillColor(0);
  c_v6->SetGrid();
  c_v6->GetFrame()->SetFillColor(0);
  c_v6->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v6=c_v6->DrawFrame(0,0,4000,ymax);
  mg_v6->SetXTitle("ADC Mean");
  mg_v6->SetYTitle("Npe");
  mg_v6->SetTitle("Mean-ped/RMS-pedRMS * corr");
    
  TGraph *gain_v6=new TGraph(npoints,hist_meanped,npe_mean_rms);
  gain_v6->SetMarkerColor(kBlue);
  gain_v6->SetMarkerStyle(20);
  gain_v6->SetMarkerSize(1);
  gain_v6->Draw("P");
  
  TF1 *func_v6=new TF1("func_v6",f_line,0,4000,2);
  func_v6->SetLineColor(2);
  func_v6->SetLineWidth(2);
  func_v6->FixParameter(0,0);
  
  gain_v6->Fit("func_v6","WR");
  c_v6->Print("8in_mean_rms_adc.png");
  
  
  
}

