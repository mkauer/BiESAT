/*************************************************************
A FUNCTIN TO READ IN VME OR CAMAC FILES INTO A HISTOGRAM

version: 5
build:   08.10.31

Matt Kauer
*************************************************************/

TH1F* readfile(char *source, Float_t ped=0, Int_t chan=1, Int_t daq=0, char *name="Name")
{
  // VME Part
  if(daq==0){
    const Int_t numchans=16;  
    const Int_t binmax=4000;
    TH1F *hadc=new TH1F(name,source,binmax,0,binmax);
    ifstream infile;
    infile.open(source);
    if(infile.is_open()){
      cout<<"\n\tREADING IN FILE ---> "<<source<<"\n";
      Int_t adc[numchans];
      Int_t j=0, events;  
      Int_t n=0;	      
      while(infile>>n){
	if(j==numchans){
	  j=0; 
	}
	adc[j]=n;
	if(j==chan){
	  if(n<=binmax) hadc->Fill(adc[j]-ped); 
	}  
	j++;		
      }
      cout<<"\n\tDONE READING FILE ---> "<<source<<"\n\n";
      infile.close();  
    }
    else cout<<"\n********** Could Not Read in File --> "<<source<<endl;
  }
  
  // CAMAC Part
  if(daq==1){ 
    const Int_t numchans=12;  
    const Int_t binmax=2000;
    TH1F *hadc=new TH1F(name,source,binmax,0,binmax);
    FILE *infile = fopen(source,"r");
    cout<<"\n\tREADING IN FILE ---> "<<source<<"\n\n";
    Int_t adc[numchans];
    char cevent[256], cmodule[256];
    Int_t evnum, ncols;
    while(1){
      ncols = fscanf(infile,"\n%5c %i", &cevent,&evnum);
      if(ncols < 0) break;
      fscanf(infile,"\n%4c", &cmodule);
      for(Int_t j=0; j<12; j++){
	fscanf(infile," %i", &adc[j]);
      }
      if(adc[chan]<=binmax) hadc->Fill(adc[chan]-ped);
      //hadc->Fill(adc[chan]-ped);
    }
    cout<<"\n\tDONE READING FILE ---> "<<source<<"\n\n";
    fclose(infile);  
  }
  return hadc;
}


