/*************************************************************
A FUNCTIN TO READ IN RAW DATA INTO A ROOTFILE

version: 1
build:   08.11.05

Matt Kauer
*************************************************************/
gROOT->Reset();
gROOT->Reset();

Int_t numchans=16;

Int_t binmax=4000;
Float_t histmin=0;
Float_t histmax=4000;

Float_t datmax=4000;

void data2root(string file2="run195", Int_t chan=0)
{
  beautify();
  
  ifstream infile2,infile3;
  infile2.open((file2).c_str());
  infile3.open((file2).c_str());
  
  TString rootoutfile=file2;
  TFile *outfile=gROOT->FindObject(rootoutfile+".root");
  if(outfile) outfile->Close();
  outfile=new TFile(rootoutfile+".root","RECREATE");
  
  TH1F *beta=new TH1F("beta","beta",binmax,histmin,histmax);
    
  Int_t j;
  Float_t n,junk;	      
  
  cout<<"READING IN FILE ==> "<<file2<<" ==> beta "<<endl;
  j=0;
  n=0;	      
  //while(infile2>>n){
  //  if(j==numchans) j=0;
  //  if(j==chan) beta->Fill(n);
  //  //junk=n;
  //  j++;
  //}
  
  //beta->Draw();
  //infile2.close();
  //infile2.open((file2).c_str());
  //infile2.seekg(ios::beg);
  //infile2.seekg(0);
  infile2.clear();
  infile2.seekg(0,ios::beg);
  
  TCanvas *look=new TCanvas("look","look",0,0,800,600);
  look->Divide(1,2,.02,.02,0);
  look->cd(1);
  beta->Draw();
  look->Update();
  
  TH1F *beta2=new TH1F("beta2","beta2",200,0.050,2.00);
  cout<<"READING IN FILE ==> "<<file2<<" ==> beta2 "<<endl;
  j=0;
  n=0;	      
  while(infile2>>n){
    //cout<<n;
    if(j==numchans) j=0;
    if(j==chan) beta2->Fill((n-30)/1570);
    //if(j==numchans) j=0;
    j++;		
  }
  
  look->cd(2);
  beta2->Draw();
  look->Update();
  
  beta->Write();
  beta2->Write();
 
  infile2.close();
  
  //gSystem->Sleep(2000);
  //gROOT->Delete("look");
  //delete look;
  //outfile->Close("R");
  
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"CREATED A ROOTFILE CALLED ==>  "<<TString(rootoutfile+".root")<<endl;
  cout<<"=====================================================================\n"<<endl;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void beautify(){
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

