// modified: 08.11.04
void config(){
  
  daq        = 0;         // 0=vme(4000 bins) : 1=camac=(2000 bins)
  chan       = 0;         // 0-7==lowRes : 8-15==hiRes
  maxbin     = 4000;      // maximum range you want in your histos
  
  rootfile   = "run195.root";  // obvious?
  pedrun     = "ped";     // run-name or hist-name in rootfile
  betarun    = "beta";    // run-name or hist-name in rootfile
  gammarun   = "gamma";   // run-name or hist-name in rootfile
  
  gamma      = 650;       // the 400 kev gamma peak position
  beta       = 1600;      // the 976 kev beta peak position
  
  lowres     = 0;         // if resolution is >8% set this true
  thick      = 0;         // if using thick scint, set this true
  fitopt     = "";        // global fit options (chi^2 or likelyhood etc)
  rebin      = 4;         // rebinning of the histrograms
  Eloss      = 14.0;      // energy loss of electrons in air

}

