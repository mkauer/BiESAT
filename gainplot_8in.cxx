#include "./include/manip_histos.hxx"
#include "./include/fittingFuncs.hxx"

void gainplot_8in(){
  
  beautify();
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin  (0.10);
  gStyle->SetPadTopMargin   (0.11);
  gStyle->SetPadRightMargin (0.06);
  
  Int_t select=0;
  
  if(select==0){
    const Int_t npoints=4;
    Double_t volts[npoints]={1700.,1800.,1900.,2000.};
    Double_t voltserr[npoints]={1.,1.,1.,1.};
    
    Double_t gain[npoints]={463969.,589339.,742587.,940410.};
    Double_t err=0.01; // error in %
    Double_t gainerr[npoints]={gain[0]*err,gain[1]*err,gain[2]*err,gain[3]*err};
  }
  
  if(select==1){
    const Int_t npoints=5;
    Double_t volts[npoints]={1500.,1600.,1700.,1800.,1900.};
    Double_t voltserr[npoints]={1.,1.,1.,1.,1.};
    
    Double_t gain[npoints]={258677.,341651.,450391.,577375.,732462.};
    Double_t err=0.01; // error in %
    Double_t gainerr[npoints]={gain[0]*err,gain[1]*err,gain[2]*err,gain[3]*err,gain[4]*err};
  }
  
  TCanvas *c_v0=new TCanvas("c_v0","Results",700,500);
  c_v0->SetFillColor(0);
  c_v0->SetGrid();
  c_v0->GetFrame()->SetFillColor(0);
  c_v0->GetFrame()->SetBorderSize(12);
  
  TH1F *mg_v0;
  if(select==0) mg_v0=c_v0->DrawFrame(1600.,400000.,2100.,1100000.);
  if(select==1) mg_v0=c_v0->DrawFrame(1400.,200000.,2000.,900000.);
  mg_v0->SetXTitle("Voltage");
  mg_v0->SetYTitle("Gain");
  mg_v0->SetTitle("Gain of 8in Ham-SBA PMT");
  
  TGraphErrors *gain_v0=new TGraphErrors(npoints,volts,gain,voltserr,gainerr);
  gain_v0->SetMarkerColor(kBlue);
  gain_v0->SetMarkerStyle(20);
  gain_v0->SetMarkerSize(1);
  gain_v0->Draw("P");
  
  TF1 *func_v0;
  if(select==0) func_v0=new TF1("func_v0",f_pmgain,1650,2050,2);
  if(select==1) func_v0=new TF1("func_v0",f_pmgain,1550,1950,2);
  func_v0->SetLineColor(2);
  func_v0->SetLineWidth(2);
  
  //func_v0->SetParameter(0,1e-9);
  //func_v0->SetParLimits(0,1.e-20,1.);
  func_v0->SetParameter(1,5.123);
  func_v0->SetParLimits(1,3.,12.);
  
  gain_v0->Fit("func_v0","ME");
  c_v0->Print("gainplot_8in.png");
  
  delete func_v0;
  
}

