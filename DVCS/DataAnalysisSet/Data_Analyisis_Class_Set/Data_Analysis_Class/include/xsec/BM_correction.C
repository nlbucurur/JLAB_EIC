TH1F* BDT::BM_correction(TH1* hist, int bin, int select){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  int Nphi = hist->GetNbinsX();

  TH1F *BM_corr = new TH1F("BM_corr","",Nphi,0,360);
  TH1F *BM_corr_Sys = new TH1F("BM_corr_Sys","",Nphi,0,360);
  TH1F *ratio = new TH1F("Phi_norec","F_noRec",Nphi,0,360);
  TH1F *ratio_Sys = new TH1F("Phi_norec_Sys","F_noRec_Sys",Nphi,0,360);

  double x, y, xErr, yErr;
  int hbin;

  std::string line;
  std::ifstream inputFile;
  std::ifstream inputFile_Sys;

  //Read File with Bin-Migrated values
  std::cout<<"BM correction"<<endl;
  if(select==1)
  {
    inputFile.open(Form("Systematics_BM/bin_%i/Entries_Most.txt",bin)); // Replace with the name of the block file you want to read
    inputFile_Sys.open(Form("Systematics_BM/bin_%i/Entries_Most_Sys.txt",bin)); // Replace with the name of the block file you want to read
  }
  else
  {
    inputFile.open(Form("Systematics_BM/bin_%i/Entries_Maxi.txt",bin)); // Replace with the name of the block file you want to read
    inputFile_Sys.open(Form("Systematics_BM/bin_%i/Entries_Maxi_Sys.txt",bin)); // Replace with the name of the block file you want to read
  } 

  hbin=1;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    iss >> x >> y >> xErr >> yErr;
    BM_corr->SetBinContent(hbin, x+y);
    //BM_corr->SetBinError(hbin, 0*sqrt(1.0/(xErr + yErr))); //I remove this error as BM will account only for systematic error.
    BM_corr->SetBinError(hbin, 0);
    hbin++;
  }
  inputFile.close();
    
  hbin=1;
  while (std::getline(inputFile_Sys, line)) {
    std::istringstream iss(line);
    iss >> x >> y >> xErr >> yErr;
    BM_corr_Sys->SetBinContent(hbin, x+y);
    //BM_corr->SetBinError(hbin, 0*sqrt(1.0/(xErr + yErr))); //I remove this error as BM will account only for systematic error.
    BM_corr_Sys->SetBinError(hbin, 0);
    hbin++;
  }
  inputFile_Sys.close();

  ratio->Divide(hist, BM_corr,1.0,1.0);
  ratio_Sys->Divide(hist, BM_corr_Sys,1.0,1.0);

  for(int k=1; k<=Nphi; k++)
  {
    if(ratio->GetBinContent(k)==0)
      ratio->SetBinError(k,0);
    else
      ratio->SetBinError(k,abs(ratio->GetBinContent(k)-ratio_Sys->GetBinContent(k))/ratio->GetBinContent(k));

  }

  for (int k=1; k<=Nphi; k++)
    std::cout<<ratio->GetBinContent(k)<<" "<<ratio->GetBinError(k)<<std::endl;
  

  delete BM_corr;
  delete BM_corr_Sys;
  delete ratio_Sys;


  return ratio;
}
