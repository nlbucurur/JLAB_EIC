int BDT::get_index(double t, double Q, double x, double p)
{
  int index=0;
  for(int k=0; k<Nbins; k++)
    {
      if(k>0)
	index+=Nphibins[k-1];
	  
      if(t > bins[k][0] && t <= bins[k][1] && Q > bins[k][2] && Q <= bins[k][3] && x > bins[k][4] && x <= bins[k][5])
	{
	  for(int l=0; l<Nphibins[k]; l++)
	    {
	      if(p > (360./Nphibins[k])*l && p <= (360./Nphibins[k])*(l+1) )
		{
		  index=index + l;
		  return index;
		}
	    }
	}
    }
  std::cout<<"Error "<<t<<" "<<Q<<" "<<x<<" "<<p<<endl;
  return -1;
}

void BDT::Compute_BM_Matrix() 
{
  int N=0;
  for(int i=0; i<Nbins; i++)
  {
  N+=Nphibins.at(i);
  }
  //Line to merge the bin files
  //Merge_Bins(TBM_Sim);

  //Load merged file. It is faster than adding all files to a TChain
  TChain *Data= new TChain("pDVCS");
  Data->Add(Folder + TBM_Sim);

  int Hel;
  static vector<double> *t;
  static vector<double> *Q;
  static vector<double> *x;
  static vector<double> *p;

  static vector<double> *t_MC;
  static vector<double> *Q_MC;
  static vector<double> *x_MC;
  static vector<double> *p_MC;

  static vector<int> *flag;
  static double BDTvar;

  Data->SetBranchAddress("t_Ph",&t);
  Data->SetBranchAddress("strip_Q2",&Q);
  Data->SetBranchAddress("strip_Xbj",&x);
  Data->SetBranchAddress("Phi_Ph",&p);

  Data->SetBranchAddress("t_Ph_MC",&t_MC);
  Data->SetBranchAddress("strip_Q2_MC",&Q_MC);
  Data->SetBranchAddress("strip_Xbj_MC",&x_MC);
  Data->SetBranchAddress("Phi_Ph_MC",&p_MC);

  Data->SetBranchAddress("Helicity",&Hel);
  Data->SetBranchAddress("bestCandidateFlag",&flag);
  Data->SetBranchAddress("_strip_Nuc_BDT",&BDTvar);
  
  int entries=Data->GetEntries();
  TH1D* input_p = new TH1D("input_p", "input_p",N,0,N);
  TH1D* input_m = new TH1D("input_m", "input_m",N,0,N);
  TH1D* input_p_2 = new TH1D("input_p_2", "input_p_2",N,0,N);
  TH1D* input_m_2 = new TH1D("input_m_2", "input_m_2",N,0,N);
  input_p->Sumw2();
  input_m->Sumw2();

  RooUnfoldResponse response (N, 0, N);  

  int input_index1=1;
  int input_index2=1;
  int input_index3=1;
  int input_index4=1;
  int input_index=1;
  double nentriesbin;
  
  //Fill input BSAs
  cout << "==================================== Fill input ===================================" << endl;
  for(int i=1; i<=Nbins; i++)
    {
      //std::ifstream inputFile(Folder + Form("bin_%i/BSA_Most_Values.txt",i));
      std::ifstream inputFilep(Folder + Form("bin_%i/entries_most_p.txt",i));
      std::ifstream inputFilem(Folder + Form("bin_%i/entries_most_m.txt",i));
      std::ifstream inputFilep_2(Folder + Form("bin_%i/entries_maxi_p.txt",i));
      std::ifstream inputFilem_2(Folder + Form("bin_%i/entries_maxi_m.txt",i));
      if (!inputFilep.is_open() || !inputFilem.is_open() || !inputFilep_2.is_open() || !inputFilem_2.is_open() ) {
	std::cerr << "Error: Unable to open "<<Folder + Form("bin_%i/entries_most_X.txt",i)<< std::endl;
      }
	 
      //while (inputFile >> nentriesbin >> comma >> BSA >> comma >> Err) 
      while (inputFilep >> nentriesbin) 
	{
	  input_p->SetBinContent(input_index1, nentriesbin);
	  input_p->SetBinError(input_index1, sqrt(nentriesbin));
	  input_index1++;
	}
      while (inputFilem >> nentriesbin) 
	{
	  input_m->SetBinContent(input_index2, nentriesbin);
	  input_m->SetBinError(input_index2, sqrt(nentriesbin));
	  input_index2++;
	}
      while (inputFilep_2 >> nentriesbin) 
	{
	  input_p_2->SetBinContent(input_index3, nentriesbin);
	  input_p_2->SetBinError(input_index3, sqrt(nentriesbin));
	  input_index3++;
	}
      while (inputFilem_2 >> nentriesbin) 
	{
	  input_m_2->SetBinContent(input_index4, nentriesbin);
	  input_m_2->SetBinError(input_index4, sqrt(nentriesbin));
	  input_index4++;
	}

      inputFilep.close();	  
      inputFilem.close();	  
      inputFilep_2.close();	  
      inputFilem_2.close();	  
    }
    
  TH1 *BA_in= input_p->GetAsymmetry(input_m);
  TH1 *BA_in_2= input_p_2->GetAsymmetry(input_m_2);


  //Fill Matrix
  cout << "==================================== Fill Matrix ===================================" << endl;
  int binx, biny;
  int col, row;
  for(int i=0; i<entries; i++)
    {
      printProgress(i*1.0/entries);
      Data->GetEntry(i);
      for(int j=0; j<flag->size(); j++)
	{
	  if(flag->at(j)==1 && BDTvar > BDT_value)
	    {
	      //Get index of the event, rows are generated, columns are reconstructed
	      //Generated has to be in columns (left to right), reconstructed in rows (top to bottom)
	      row=get_index(t->at(j), Q->at(j), x->at(j), p->at(j));
	      col=get_index(t_MC->at(j), Q_MC->at(j), x_MC->at(j), p_MC->at(j));  

	      for(int k=0; k<Nbins; k++)
		{
		  if(t->at(j) > bins[k][0] && t->at(j) < bins[k][1] && Q->at(j) > bins[k][2] && Q->at(j) < bins[k][3] && x->at(j) > bins[k][4] && x->at(j) < bins[k][5])
			{
		    binx=k;
			}
		  if(t_MC->at(j) > bins[k][0] && t_MC->at(j) < bins[k][1] && Q_MC->at(j) > bins[k][2] && Q_MC->at(j) < bins[k][3] && x_MC->at(j) > bins[k][4] && x_MC->at(j) < bins[k][5])
			{
		    biny=k;
			}
		}
    
	      if(abs(binx-biny)<9999)//Only 2 bins away
			{
			response.Fill(row*1.0, col*1.0);
			}
	    }
	}
    }


  cout << "\n==================================== UNFOLD ===================================" << endl;
  RooUnfoldBayes   unfold_p (&response, input_p, 50);
  TH1D* output_p= (TH1D*) unfold_p.Hunfold();

  RooUnfoldBayes   unfold_m (&response, input_m, 50);
  TH1D* output_m= (TH1D*) unfold_m.Hunfold();

  RooUnfoldBayes   unfold_p_2 (&response, input_p_2, 50);
  TH1D* output_p_2= (TH1D*) unfold_p_2.Hunfold();

  RooUnfoldBayes   unfold_m_2 (&response, input_m_2, 50);
  TH1D* output_m_2= (TH1D*) unfold_m_2.Hunfold();
  
  RooUnfoldSvd   unfold_p_sys (&response, input_p);
  TH1D* output_p_sys= (TH1D*) unfold_p_sys.Hunfold();

  RooUnfoldSvd   unfold_m_sys (&response, input_m);
  TH1D* output_m_sys= (TH1D*) unfold_m_sys.Hunfold();

  RooUnfoldSvd   unfold_p_2_sys (&response, input_p_2);
  TH1D* output_p_2_sys= (TH1D*) unfold_p_2_sys.Hunfold();

  RooUnfoldSvd   unfold_m_2_sys (&response, input_m_2);
  TH1D* output_m_2_sys= (TH1D*) unfold_m_2_sys.Hunfold();

  //RG-A data has the helicity flipped
  TH1 *BA_out= output_m->GetAsymmetry(output_p);
  TH1 *BA_out_2= output_m_2->GetAsymmetry(output_p_2);

  TH1 *BA_out_sys= output_m_sys->GetAsymmetry(output_p_sys);
  TH1 *BA_out_2_sys= output_m_2_sys->GetAsymmetry(output_p_2_sys);

  auto* R = response.HresponseNoOverflow();
  auto* c1 = new TCanvas();
  R->SetStats(0);
  //gPad->SetLogz();
  R->Draw("colz");
  c1->Draw();
  c1->SaveAs("response.png");
  
  //Write output
  cout << "==================================== Write output ===================================" << endl;
  std::ofstream outFile1;
  std::ofstream outFile2;
  std::ofstream outFile3;
  std::ofstream outFile4;

  std::ofstream outFile5;
  std::ofstream outFile6;
  std::ofstream outFile7;
  std::ofstream outFile8;
  input_index=1;
  gSystem->Exec("mkdir -p Systematics_BM");
  for(int i=1; i<=Nbins;i++)
    {
      gSystem->Exec(Form("mkdir -p Systematics_BM/bin_%i",i));
      outFile1.open(Form("Systematics_BM/bin_%i/BSA_Most_Values.txt",i));
      outFile2.open(Form("Systematics_BM/bin_%i/BSA_Maxi_Values.txt",i));
      outFile3.open(Form("Systematics_BM/bin_%i/Entries_Most.txt",i));
      outFile4.open(Form("Systematics_BM/bin_%i/Entries_Maxi.txt",i));

      outFile5.open(Form("Systematics_BM/bin_%i/BSA_Most_Values_Sys.txt",i));
      outFile6.open(Form("Systematics_BM/bin_%i/BSA_Maxi_Values_Sys.txt",i));
      outFile7.open(Form("Systematics_BM/bin_%i/Entries_Most_Sys.txt",i));
      outFile8.open(Form("Systematics_BM/bin_%i/Entries_Maxi_Sys.txt",i));

      for(int j=1; j<=Nphibins[i-1];j++)
	{
	  outFile1<<(output_p  ->GetBinContent(input_index) + output_m  ->GetBinContent(input_index))<<", "<<BA_out  ->GetBinContent(input_index)<<", "<<BA_out  ->GetBinError(input_index)<<endl;	  
	  outFile2<<(output_p_2->GetBinContent(input_index) + output_m_2->GetBinContent(input_index))<<", "<<BA_out_2->GetBinContent(input_index)<<", "<<BA_out_2->GetBinError(input_index)<<endl;	  
	  outFile3<< output_p  ->GetBinContent(input_index) <<" "<< output_m  ->GetBinContent(input_index)<<" "<< output_p  ->GetBinError(input_index) <<" "<< output_m  ->GetBinError(input_index)<<endl;
	  outFile4<< output_p_2->GetBinContent(input_index) <<" "<< output_m_2->GetBinContent(input_index)<<" "<< output_p_2->GetBinError(input_index) <<" "<< output_m_2->GetBinError(input_index)<<endl;

	  outFile5<<(output_p_sys  ->GetBinContent(input_index) + output_m_sys  ->GetBinContent(input_index))<<", "<<BA_out_sys  ->GetBinContent(input_index)<<", "<<BA_out_sys  ->GetBinError(input_index)<<endl;	  
	  outFile6<<(output_p_2_sys->GetBinContent(input_index) + output_m_2_sys->GetBinContent(input_index))<<", "<<BA_out_2_sys->GetBinContent(input_index)<<", "<<BA_out_2_sys->GetBinError(input_index)<<endl;	  
	  outFile7<< output_p_sys  ->GetBinContent(input_index) <<" "<< output_m_sys  ->GetBinContent(input_index)<<" "<< output_p_sys  ->GetBinError(input_index) <<" "<< output_m_sys  ->GetBinError(input_index)<<endl;
	  outFile8<< output_p_2_sys->GetBinContent(input_index) <<" "<< output_m_2_sys->GetBinContent(input_index)<<" "<< output_p_2_sys->GetBinError(input_index) <<" "<< output_m_2_sys->GetBinError(input_index)<<endl;
    //std::cout<<BA_out->GetBinContent(input_index)<<" "<<BA_out_2->GetBinContent(input_index)<<" "<<BA_in->GetBinContent(input_index)<<" "<<BA_in_2->GetBinContent(input_index)<<endl; 
	  input_index++;
	}
      outFile1.close();
      outFile2.close();
      outFile3.close();
      outFile4.close();

      outFile5.close();
      outFile6.close();
      outFile7.close();
      outFile8.close();
    }
  std::cout<<input_index<<endl;
  
  delete input_p;
  delete input_m;
  delete input_p_2;
  delete input_m_2;
  delete output_p;
  delete output_m;
  delete output_p_2;
  delete output_m_2;
  delete BA_out;
  delete BA_out_2;

  delete output_p_sys;
  delete output_m_sys;
  delete output_p_2_sys;
  delete output_m_2_sys;
  delete BA_out_sys;
  delete BA_out_2_sys;

  delete Data;

  
}
