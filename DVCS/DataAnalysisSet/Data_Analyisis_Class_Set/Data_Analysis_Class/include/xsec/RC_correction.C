std::vector<vector<double>> BDT::ReadMeansFile(TString file) {
  // Variables to store the column values
  double col1, col2, col3, col4, col5, col6, col7, col8;
  std::vector<vector<double>> means;
  std::vector<double> line;
  means.clear();

  // Open the input file
  std::ifstream inputFile(file);
  if (!inputFile.is_open()) {
    std::cerr << "Error: Unable to open " << file<< std::endl;
    return means;
  }

    
  // Read and discard the first line (header)
  std::string header;
  std::getline(inputFile, header);

  // Read and process the remaining lines
  while (inputFile >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8) 
    {
      line.clear();
      line.push_back(col1);
      line.push_back(col2);
      line.push_back(col3);
      line.push_back(col4);
      line.push_back(col5);
      line.push_back(col6);
      line.push_back(col7);
      line.push_back(col8);
      means.push_back(line);
    }

  // Close the input file
  inputFile.close();

  return means;

}


TH1F* BDT::F_4Dto3D(double xmean, double Qmean, double tmean, int Nphi, int select) {

  std::vector<vector<double>> means;
  char command[200];
  int n;
  TH1F *ratio_bin = new TH1F("ratio_bin","F_423",Nphi,0,360);
  
  if(select==1)
    means = ReadMeansFile(Folder + TString("../means_most.txt"));
  else
    means = ReadMeansFile(Folder + TString("../means_maxi.txt"));

  for(int m=0; m<Nphi; m++)
    {
      printProgress(m*1.0/Nphi);

      double pm=means.at(m).at(0);
      double tm=means.at(m).at(2);
      double Qm=means.at(m).at(4);
      double xm=means.at(m).at(6);
      
      double Phimean = (360./Nphi)*(m+0.5);
      std::string base_dir(dir.Data());
      n=sprintf (command, "python3 %s/include/xsec/4Dto3D.py %f %f %f %f %f %f %f %f",base_dir.c_str(), xmean, Qmean, tmean, Phimean, xm, Qm, tm, pm);

      //printf ("[%s] is a string %d chars long\n",command,n);
      // Open a pipe and execute the command
      FILE* pipe = popen(command, "r");
      if (!pipe) {
	fprintf(stderr, "Error executing Python script.\n");
	std::cout<<command<<endl;
	return ratio_bin;
      }

      // Read the output from the pipe
      float F_bin;
      int readCount = fscanf(pipe, "%f", &F_bin);

      if (readCount != 1) {
	fprintf(stderr, "Error reading numbers from Python script.\n");
	std::cout<<command<<endl;
      }
      // Close the pipe
      pclose(pipe);
  
  
      //std::cout<<"corr factor "<<m+1<<" "<<F_bin<<endl;
      ratio_bin->SetBinContent(m+1,F_bin);
      ratio_bin->SetBinError(m+1,0.0);
    }
  
  std::cout<<"4D to 3D correction"<<endl;
  for(int k=1; k<=Nphi; k++)
    std::cout<<ratio_bin->GetBinContent(k)<<std::endl;
  
  return ratio_bin;
}


//TH1F* BDT::F_Bin(int bin, double xmean, double Qmean, double tmean, int select) {
TH1F* BDT::F_Bin(int bin, TH1F* factor_t, TH1F* factor_Q, TH1F* factor_x, int select) {
    int Nphi=Nphibins[bin-1];
  std::vector<vector<double>> means;
  char command1[200];
  char command2[200];
  int n1, n2;
  TH1F *ratio_bin1 = new TH1F("ratio_bin1","F_bin1",Nphi,0,360);
  TH1F *ratio_bin2 = new TH1F("ratio_bin2","F_bin2",Nphi,0,360);
  
  if(select==1)
    means = ReadMeansFile(Folder + TString("../means_most.txt"));
  else
    means = ReadMeansFile(Folder + TString("../means_maxi.txt"));

  double pm, tm, Qm, xm;
  
  for(int m=0; m<Nphi; m++)
    {    
      pm=means.at(m).at(0);
      tm=means.at(m).at(2) * factor_t->GetBinContent(m+1);
      Qm=means.at(m).at(4) * factor_Q->GetBinContent(m+1);
      xm=means.at(m).at(6) * factor_x->GetBinContent(m+1);

      printProgress(m*1.0/Nphi);
      std::string base_dir(dir.Data());

      if(pm==0 || tm==0 || Qm==0 || xm==0)
      {
        ratio_bin1->SetBinContent(m+1,0);
        ratio_bin2->SetBinContent(m+1,0);
        ratio_bin1->SetBinError(m+1,0); 
        ratio_bin2->SetBinError(m+1,0);
   
      }
      else
      {
        //n=sprintf (command, "python3 %s/include/xsec/KM15.py %f %f %f %f %f %f %f %f %f %f %f %f",base_dir.c_str(),bins[bin-1][0],bins[bin-1][1],bins[bin-1][2],bins[bin-1][3],bins[bin-1][4],bins[bin-1][5], (360./Nphi)*m, (360./Nphi)*(m+1), xmean, tmean, Qmean, (360./Nphi)*(m+0.5));
        n1=sprintf (command1, "python3 %sinclude/xsec/KM15.py %f %f %f %f %f %f %f %f %f %f %f %f",base_dir.c_str(),bins[bin-1][0],bins[bin-1][1],bins[bin-1][2],bins[bin-1][3],bins[bin-1][4],bins[bin-1][5], (360./Nphi)*m, (360./Nphi)*(m+1), xm, tm, Qm, pm);
        n2=sprintf (command2, "python3 %sinclude/xsec/KM10.py %f %f %f %f %f %f %f %f %f %f %f %f",base_dir.c_str(),bins[bin-1][0],bins[bin-1][1],bins[bin-1][2],bins[bin-1][3],bins[bin-1][4],bins[bin-1][5], (360./Nphi)*m, (360./Nphi)*(m+1), xm, tm, Qm, pm);

        //printf ("[%s] is a string %d chars long\n",command,n);
        // Open a pipe and execute the command
        FILE* pipe1 = popen(command1, "r");
        FILE* pipe2 = popen(command2, "r");
        if (!pipe1 || !pipe2) 
        {
	        fprintf(stderr, "Error executing Python script.\n");
	        std::cout<<command1<<" "<<command2<<endl;
	        return ratio_bin1;
        }

        // Read the output from the pipe
        float F_bin1;
        float F_bin2;
        int readCount1 = fscanf(pipe1, "%f", &F_bin1);
        int readCount2 = fscanf(pipe2, "%f", &F_bin2);

        if (readCount1 != 1 || readCount2 != 1) 
        {
	        fprintf(stderr, "Error reading numbers from Python script.\n");
	        std::cout<<command1<<" "<<command2<<endl;
        }
        // Close the pipe
        pclose(pipe1);
        pclose(pipe2);
  
        //std::cout<<"corr factor "<<m+1<<" "<<F_bin<<endl;
        ratio_bin1->SetBinContent(m+1,F_bin1);
        ratio_bin2->SetBinContent(m+1,F_bin2);
        ratio_bin1->SetBinError(m+1,abs(F_bin1-F_bin2)/F_bin1); //delta_F_bin/F_bin //This is a systematic error
        ratio_bin2->SetBinError(m+1,0.0);
      }
    }
  
    std::cout<<"\nBin correction"<<endl;
    for(int k=1; k<=Nphi; k++)
      std::cout<<ratio_bin1->GetBinContent(k)<<" "<<ratio_bin2->GetBinContent(k)<<endl;

  delete ratio_bin2;
      
  return ratio_bin1;
}


TH1F* BDT::F_RC(int Nphi) {
  TChain *Data_Tree= new TChain("tree");
  Data_Tree->Add(MC_DVCS_RC);
  
  TH1F *xborn = new TH1F("xborn","",Nphi,0,360);
  TH1F *xrad  = new TH1F("xrad","",Nphi,0,360);

  TH1F *nevs = new TH1F("nevs","",Nphi,0,360);
  TH1F *xsc2_born  = new TH1F("xsc2_born ","",Nphi,0,360);
  TH1F *xsc2_rad  = new TH1F("xsc2_rad ","",Nphi,0,360);

  TH1F *var_born1 = new TH1F("var_born1","",Nphi,0,360);
  TH1F *var_rad1 = new TH1F("var_rad1","",Nphi,0,360);
  TH1F *var_born2 = new TH1F("var_born2","",Nphi,0,360);
  TH1F *var_rad2 = new TH1F("var_rad2","",Nphi,0,360);
  TH1F *var_born = new TH1F("var_born","",Nphi,0,360);
  TH1F *var_rad = new TH1F("var_rad","",Nphi,0,360);

  xborn->Sumw2();
  xrad->Sumw2();
  nevs->Sumw2();
  xsc2_born->Sumw2();
  xsc2_rad->Sumw2();
  var_born1->Sumw2();
  var_born2->Sumw2();
  var_born->Sumw2();
  var_rad1->Sumw2();
  var_rad2->Sumw2();
  var_rad->Sumw2();

  TH1F *ratio_rad = new TH1F("ratio_rad","F_rad",Nphi,0,360);
  ratio_rad->Sumw2();

  TString String_cut = cut_bin.GetTitle();
  String_cut.ReplaceAll("bestCandidateFlag==1 && ", "");
  String_cut.ReplaceAll("strip_Q2", "Q2_meas");
  String_cut.ReplaceAll("strip_Xbj", "xB_meas");
  String_cut.ReplaceAll("t_Ph", "t_meas");
  TCut cut = TCut(String_cut);

  var_born1->Divide(xsc2_born,xborn,1,1);
  var_born2->Divide(nevs,xborn,1,1);
  var_born2->Multiply(var_born2, var_born2,1,1);
  var_born->Add(var_born1,var_born2,1,-1);

  var_rad1->Divide(xsc2_rad,xborn,1,1);
  var_rad2->Divide(nevs,xborn,1,1);
  var_rad2->Multiply(var_rad2, var_rad2,1,1);
  var_rad->Add(var_rad1,var_rad2,1,-1);
  
  Data_Tree->Project("xborn", "phi_meas", cut * TCut("1.0/xsec_born"));
  Data_Tree->Project("xrad" , "phi_meas", cut * TCut("1.0/xsec_rad") );

  //std::cout<<cut * TCut("1.0/xsec_born").GetTitle()<<" "<<Nphi<<endl;

  Data_Tree->Project("nevs" , "phi", cut );

  Data_Tree->Project("xsc2_born" , "phi_meas", cut * TCut("xsec_born"));
  Data_Tree->Project("xsc2_rad" , "phi_meas", cut * TCut("xsec_rad"));


  //Here I do not divide by the number of events on each bin as it cancels in the ratio.
  
  ratio_rad->Divide(xrad, xborn,1,1);

  for(int k=1; k<=Nphi; k++)
  {
    if(ratio_rad->GetBinContent(k)==0)
      ratio_rad->SetBinError(k,0);
    else
      ratio_rad->SetBinError( k,ratio_rad->GetBinContent(k)*sqrt(pow(xborn->GetBinError(k)/xborn->GetBinContent(k),2) + pow(xrad->GetBinError(k)/xrad->GetBinContent(k),2)) );
  }
    //ratio_rad->SetBinError(k,sqrt(1.0/xborn->GetBinContent(k) + 1.0/xrad->GetBinContent(k)));
    //ratio_rad->SetBinError(k,sqrt(var_born->GetBinContent(k)/xborn->GetBinContent(k) + var_rad->GetBinContent(k)/xrad->GetBinContent(k)));

  std::cout<<"RC correction"<<endl;
  for(int k=1; k<=Nphi; k++)
    std::cout<<ratio_rad->GetBinContent(k)<<" "<<ratio_rad->GetBinError(k)<<" "<<ratio_rad->GetBinContent(k)*sqrt(xborn->GetBinError(k)/xborn->GetBinContent(k))<<" "<<ratio_rad->GetBinContent(k)*sqrt(xrad->GetBinError(k)/xrad->GetBinContent(k))<<std::endl;

  delete xsc2_born;
  delete xsc2_rad;
  delete nevs;
      
  delete var_born1;
  delete var_born2;
  delete var_born;
  delete var_rad1;
  delete var_rad2;
  delete var_rad;

  delete xborn;
  delete xrad;
  delete Data_Tree;

  return ratio_rad;
}

/*
TH1F* BDT::compute_kin_shift_factor(int Nphi, TString var)
{
  TString var_meas = var + TString("_meas");

  TString String_cut = cut_bin.GetTitle();
  String_cut.ReplaceAll("bestCandidateFlag==1 && ", "");
  String_cut.ReplaceAll("strip_Q2", "Q2_meas");
  String_cut.ReplaceAll("strip_Xbj", "xB_meas");
  String_cut.ReplaceAll("t_Ph", "t_meas");
  TCut cut = TCut(String_cut);

  TChain *Data_Tree= new TChain("tree");
  Data_Tree->Add(MC_DVCS_RC);

  TH1F *ratio = new TH1F("ratio"+var,"ratio"+var,Nphi,0,360);

  for(int i=0; i<Nphi;i++)
  {
    TH1F *hvar = new TH1F("hvar","",100,0,0);
    TH1F *hvar_meas  = new TH1F("hvar_meas","",100,0,0);

    Data_Tree->Project("hvar" , var, cut + TCut(Form("phi_meas>%f && phi_meas<%f",360.*i/Nphi,360.*(i+1)/Nphi)));
    Data_Tree->Project("hvar_meas" , var_meas, cut + TCut(Form("phi_meas>%f && phi_meas<%f",360.*i/Nphi,360.*(i+1)/Nphi)));  
    
    ratio->SetBinContent(i+1,hvar->GetMean()/hvar_meas->GetMean());
    ratio->SetBinError(i+1,ratio->GetBinContent(i+1)*sqrt(pow(hvar->GetMeanError()/hvar->GetMean(),2) + pow(hvar_meas->GetMeanError()/hvar_meas->GetMean(),2)));

    delete hvar;
    delete hvar_meas;
    }

  delete Data_Tree;

  return ratio;
}
*/



TH1F* BDT::compute_kin_shift_factor(int Nphi, TString var)
{
  std::cout<<"\nComputing kinematic correction of "<<var<<endl;
  
  TH1F *ratio = new TH1F("ratio_"+var,"ratio"+var,Nphi,0,360);
  for(int i=1; i<=Nphi; i++)
  {
    ratio->SetBinContent(i,1);
    ratio->SetBinError(i,0);    
  }
  return ratio; 

  TChain *Data_Tree= new TChain("pDVCS");
  Data_Tree->Add(Folder + TString("../../") + TRC_Sim); 
  
  TChain *MCData_Tree= new TChain("pDVCS");
  MCData_Tree->Add(MC_RC_Sim);

  for(int i=0; i<Nphi;i++)
  {
    printProgress(i*1.0/Nphi);

    TH1F *hvar = new TH1F("hvar","",100,0,0);
    TH1F *hvar_meas  = new TH1F("hvar_meas","",100,0,0);

    MCData_Tree->Project("hvar" , var, cut_bin  + TCut(Form("Phi_Ph>%f && Phi_Ph<%f",360.*i/Nphi,360.*(i+1)/Nphi)));
    Data_Tree->Project("hvar_meas" , var, cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)) + TCut(Form("Phi_Ph>%f && Phi_Ph<%f",360.*i/Nphi,360.*(i+1)/Nphi)));  
    
    ratio->SetBinContent(i+1,hvar->GetMean()/hvar_meas->GetMean());
    ratio->SetBinError(i+1,ratio->GetBinContent(i+1)*sqrt(pow(hvar->GetMeanError()/hvar->GetMean(),2) + pow(hvar_meas->GetMeanError()/hvar_meas->GetMean(),2)));

    delete hvar;
    delete hvar_meas;
    }

  delete Data_Tree;
  delete MCData_Tree;

  return ratio;
}


TH1F* BDT::compute_kin_shift_factor(TString var)
{
  std::cout<<"\nComputing kinematic correction (3D) of "<<var<<endl;

  TH1F *ratio = new TH1F("ratio3D_"+var,"ratio"+var,1,0,360);
  ratio->SetBinContent(1,1);
  ratio->SetBinError(1,0);
  return ratio;

  TChain *Data_Tree= new TChain("pDVCS");
  Data_Tree->Add(Folder + TString("../../") + TRC_Sim); 
  
  TChain *MCData_Tree= new TChain("pDVCS");
  MCData_Tree->Add(MC_RC_Sim);

  TH1F *hvar = new TH1F("hvar","",100,0,0);
  TH1F *hvar_meas  = new TH1F("hvar_meas","",100,0,0);

  MCData_Tree->Project("hvar" , var, cut_bin );
  Data_Tree->Project("hvar_meas" , var, cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)) );  
    
  ratio->SetBinContent(1,hvar->GetMean()/hvar_meas->GetMean());
  ratio->SetBinError(1,ratio->GetBinContent(1)*sqrt(pow(hvar->GetMeanError()/hvar->GetMean(),2) + pow(hvar_meas->GetMeanError()/hvar_meas->GetMean(),2)));

  delete hvar;
  delete hvar_meas;

  delete Data_Tree;
  delete MCData_Tree;

  return ratio;
}
