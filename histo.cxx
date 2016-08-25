/**************************************************************************
This macro will just make a histogram of the camac data channel that you
specify with the "chan" variable.

version: 4
build:   09.04.23

Matt Kauer
**************************************************************************/

//#include "_readinData.hxx"
//#include "_userFuncs.hxx"
#include "./include/manip_files.hxx"
#include "./include/manip_histos.hxx"


void histo(char *datfile="", Int_t resolution=0, Int_t lowchan=1, Int_t ped=0, Int_t daq=0){
  cout<<"-----------------------------------------------------------------------\n";
  cout<<"USAGE ==> histo.cxx (filename, res_select, channel, pedestal, daq_type)\n";
  cout<<"-----------------------------------------------------------------------\n";
  
  beautify();
  if(string(datfile)==""){
    cout<<"\nERROR: Specify a filename to use!\n\n";
    return;
  }

  
  //---- VME ----//
  if(daq==0){
    if(resolution==0 || resolution==1){
      TH1F *LowRes = (TH1F*)readfile(datfile,ped,lowchan,daq,"LowRes");
      TCanvas *c10 = new TCanvas("c10","c10",700,500);
      LowRes->Draw();
      c10->Update();
    }
    if(resolution==0 || resolution==2){
      TH1F *HiRes = (TH1F*)readfile(datfile,ped,lowchan+8,daq,"HiRes");
      TCanvas *c11 = new TCanvas("c11","c11",700,500);
      HiRes->Draw();
      c11->Update();
    }
  }
  
  //---- CAMAC ----//
  if(daq==1){
    TH1F *camac = readfile(datfile,ped,lowchan,daq,"camac");
    TCanvas *c12 = new TCanvas("c12","c12",700,500);
    camac->Draw();
    c12->Update();   
  }
  
}

