//***********************************************************
// SOME USEFUL FUNCTION TO MAKE LIFE EASY
//
// version: 5
// build:   08.10.31
//
// Matt Kauer
//***********************************************************

//--- SET ROOT DEFAULTS TO LOOK GOOD
void beautify()
{
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(2211);
  gStyle->SetOptFit(112);
  gStyle->SetStatFormat("5.4g");
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatW(0.13);
  gStyle->SetStatH(0.08);
  gStyle->SetPalette(1);
}

//--- RETURN A CLONED HISTOGRAM BETTER THAN ROOT
TH1F* cloneHist(TH1F *inhist,char *name="Name",char *title="Title")
{
  Int_t nbins=inhist->GetXaxis()->GetNbins();
  Int_t hmin=inhist->GetXaxis()->GetBinLowEdge(1);
  Int_t hmax=inhist->GetXaxis()->GetBinUpEdge(nbins);
  Double_t content;
  Double_t nevents=inhist->GetEntries();
  
  //cout<<"\n\n\t nbins = "<<nbins<<"  hmin = "<<hmin<<"  hmax = "<<hmax<<"\n\n";
  
  TH1F *outhist=new TH1F(name,title,nbins,hmin,hmax);
  
  for(Int_t i=0;i<nbins;i++){
    content=inhist->GetBinContent(i);
    outhist->SetBinContent(i,content);
  }  
  outhist->SetEntries(nevents);

  return outhist;
}

//--- SHIFT A HISTOGRAM
void shiftHist(TH1F *hist, Double_t double_shift)
{
  Int_t bins=hist->GetXaxis()->GetNbins();
  Double_t xmin=hist->GetXaxis()->GetBinLowEdge(1);
  Double_t xmax=hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift=double_shift*(bins/(xmax-xmin));
  Int_t shift;
  if(double_shift<0) shift=TMath::FloorNint(double_shift);
  if(double_shift>0) shift=TMath::CeilNint(double_shift);
  if(shift>0){
    for(Int_t i=1; i<=bins; i++){
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
  }
  if(shift<0){
    for(Int_t i=bins; i>0; i--){
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }    
  }
} 

//--- RETURN A CLONED HISTOGRAM BETTER THAN ROOT
TH1F* resizeHist(TH1F *inhist,Int_t nhmin,Int_t nhmax,char *name="Name",char *title="Title")
{
  Int_t bins=inhist->GetXaxis()->GetNbins();
  Float_t hmin=inhist->GetXaxis()->GetBinLowEdge(1);
  Float_t hmax=inhist->GetXaxis()->GetBinUpEdge(bins);
  Int_t content;
  Int_t nevents=0;
  
  cout<<"\n\t bins = "<<bins<<"  hmin = "<<hmin<<"  hmax = "<<hmax<<"\n";
  
  Float_t unit=hmax/bins;
  Int_t lobin,hibin,nbins;
  if(nhmin!=0) lobin=TMath::CeilNint(nhmin/unit);
  hibin=TMath::FloorNint(nhmax/unit);
  nbins=hibin-lobin;
  
  cout<<"\t nbins = "<<nbins<<"  nhmin = "<<nhmin<<"  nhmax = "<<nhmax<<"\n";
  
  TH1F *outhist=new TH1F(name,title,nbins,nhmin,nhmax);
  
  
  for(Int_t i=lobin+1;i<=hibin;i++){
    content=inhist->GetBinContent(i-lobin);
    outhist->SetBinContent(i-lobin,content);
    nevents+=content;
    //cout<<"bin = "<<i-lobin<<"  value = "<<content<<"  total events = "<<nevents<<endl;
  }  
  outhist->SetBinContent(0,0);
  outhist->SetBinContent(hibin+1,0);
  outhist->SetEntries(nevents);
  
  return outhist;
}

//--- NORMALIZE THE HISTOS
void normalizer(TH1F *hgamma,TH1F *hbeta,Bool_t rootread)
{
  Int_t imin=620*calib;
  Int_t imax=800*calib;
  Double_t beta_area,gamma_area,area_ratio;
  cout<<"\n\n\t Integrate between  "<<imin<<" - "<<imax<<" \n";
  
  if(!rootread) gamma_area = hgamma->TH1::Integral(imin,imax);
  if(rootread) gamma_area = hgamma->TH1::Integral(imin,imax,"width");
  cout<<"\t Area under gamma hist = "<<gamma_area<<" \n";
  
  if(!rootread) beta_area = hbeta->TH1::Integral(imin,imax);
  if(rootread) beta_area = hbeta->TH1::Integral(imin,imax,"width");
  cout<<"\t Area under beta hist = "<<beta_area<<" \n";
 
  area_ratio = beta_area/gamma_area;
  cout<<"\t Gamma scale factor = "<<area_ratio<<" \n\n\n";
  if(!rootread) hgamma->TH1::Scale(area_ratio);
  if(rootread) hgamma->TH1::Scale(area_ratio,"width");
  
}

// DRAW A CANVAS DISPLAYING THE OVERLAP OF TWO HISTOS
void drawMix(TH1F *hgamma, TH1F *hbeta)
{
  TCanvas *c69 = new TCanvas("c69","- Scaled Mix -",xsize,ysize);
  TH1F *mix1 = (TH1F*)cloneHist(hgamma,"mix","Betas=blue & Gammas=green");
  mix1->SetLineColor(8);
  TH1F *mix2 = (TH1F*)cloneHist(hbeta,"mix","Betas=blue & Gammas=green");
  mix2->SetLineColor(9);
  mix1->SetAxisRange(200*calib,1400*calib,"X");
  mix2->SetAxisRange(200*calib,1400*calib,"X");
  if(mix2->GetBinContent(mix2->GetMaximumBin()) > mix2->GetBinContent(mix1->GetMaximumBin())){
    mix2->Draw();
    mix1->Draw("same");
  }
   if(mix1->GetBinContent(mix1->GetMaximumBin()) > mix2->GetBinContent(mix2->GetMaximumBin())){
    mix1->Draw();
    mix2->Draw("same");
  }
  
  c69->Print("A_mixed.png");
}


//--- RETURN THE NUMBER OF COLUMNS IN A FILE
int getCol(char *datafile)
{
  int info=0;
  FILE *infile=fopen(datafile,"r");
  int temp=0,prev=0,first=1;
  while((temp=getc(infile))!=EOF){
    
    //std::cout<<dec<<temp<<"  ";
    
    // if the first character(s) is a space or tab ignore it
    if((temp==32 || temp==9) && first==1 && prev==0){
      first=0;
      prev=temp;
      continue;
    }
    first=0;
    // if there's a space or tab and -NOT- and previous space or tab
    if((temp==32 || temp==9) && prev!=32 && prev!=9){
      info++;
      prev=temp;
      continue;
    }
    // if there's a newline and -NOT- a previous space or tab
    if(temp==10 && prev!=32 && prev!=9){
      info++;
      break;
    }
    // if there's a newline and -YES- a previous space or tab
    if(temp==10 && (prev==32 || prev==9)) break;
    prev=temp;
  
  } // end of while loop
  
  fclose(infile);
  return info;
}

//--- RETURN THE NUMBER OF LINES IN A FILE
int getLines(char *datafile)
{
  int info=0;
  ifstream infile;
  infile.open(datafile);
  string line;
  while(getline(infile,line)) info++;
  infile.close();
  return info;
}


