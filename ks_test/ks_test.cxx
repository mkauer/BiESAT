
gROOT->Reset();

#include "../include/fittingFuncs.hxx"
#include "../include/manip_files.hxx"
#include "../include/manip_params.hxx"
#include "../include/manip_histos.hxx"

Int_t pmin,pmax,maxbin=4000,xsize=600,ysize=400,xmin=0,xmax=4000,rebin=1;
Float_t calib=1.632;

void ks_test(){
  
  beautify();
  
  string rldir="./real/";
  rldir="./";
  string realroot="run195";
  string wkdir="./mc/";
  string mcroot="1M_results";
  
  const Int_t party=21; // max=21
  
  //string mcroot[party]={"6.0pc_5x5","6.5pc_5x5","7.0pc_5x5"};
  
  char *pedrun="run194";
  char *gammarun="run196";
  char *betarun="run195";
  betarun="beta2";
  
  TFile *realfile=new TFile((rldir+realroot+".root").c_str(),"READ");
  TH1F *ptemp=(TH1F*)realfile->Get(pedrun);
  TH1F *gtemp=(TH1F*)realfile->Get(gammarun);
  TH1F *btemp=(TH1F*)realfile->Get(betarun);
  
  TH1F *hbeta=(TH1F*)realfile->Get(betarun);
  
  //hbeta->Sumw2();
    
  //TFile *mcfile=new TFile(mcroot.c_str(),"READ");
  //TTree *tree=new TTree("theRunTree","tree");
  //TBranch *branch=tree->GetBranch("RunStatistics");
  //TLeaf *leaf=branch->GetLeaf("Deposit_prim_scint");
  
  TH1F *test[party];
  
  TFile *mcfile=new TFile((wkdir+mcroot+".root").c_str(),"READ");
  
  ostringstream hname[party];
  ostringstream j[party];
  
  Double_t results[party];
  TCanvas *canvas[party];
  
  Float_t imin=hbeta->FindBin(0.95);
  Float_t imax=hbeta->FindBin(1.01);
  Float_t realarea,mcarea,scale;
  
  //hbeta->Sumw2();
  string title="";
  
  for(Int_t i=0;i<party;i++){
    hname[i]<<"hSmearedEnergy_"<<i;
    cout<<hname[i].str()<<endl;
    test[i]=(TH1F*)mcfile->Get(hname[i].str().c_str());
    
    j[i]<<"canvas "<<i;
    //canvas[i]=new TCanvas(j[i].str().c_str(),j[i].str().c_str(),100*i,100*i,600,400);
    canvas[i]=new TCanvas(j[i].str().c_str(),j[i].str().c_str(),10,10,600,400);
    canvas[i]->cd();
    
    mcarea=test[i]->Integral(imin,imax);
    realarea=hbeta->Integral(imin,imax);
    scale=mcarea/realarea;
    hbeta->Scale(scale);
    
    hbeta->Rebin(rebin);
    hbeta->SetLineWidth(1);
    hbeta->SetLineColor(kBlue);
    
    test[i]->Rebin(rebin);
    test[i]->SetLineWidth(1);
    test[i]->SetLineColor(kRed);
    
    //hbeta->Sumw2();
    //test[i]->Sumw2();
    
    results[i]=hbeta->KolmogorovTest(test[i],"ND");
    
    cout<<endl;
        
    test[i]->Draw();
    hbeta->Draw("SAME");
    
    title=test[i]->GetTitle();
    canvas[i]->Update();
    canvas[i]->Print((title+".png").c_str());
    
    
    
  }
  
  
  
  
  
  
}

