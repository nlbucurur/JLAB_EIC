
void BDT::Training_on_bins(TString Data, int NumEv, int bin=0)//, TCut cut, TString DVCS, TString Pi0, vector<TString> vars)
{
  int NBINS_t=3;
  int NBINS_x=3;
  //int NumEv=29000;

  /*
  TFile *File = new TFile(Data);
  TTree *ch1 = (TTree *)(File->Get("pDVCS"));
  
  TH1F *temp2 = new TH1F("temp2","Histogram",10,-200,0);
  ch1->Project("temp2", "t_Ph",cut);
  int DATASZ=temp2->GetEntries();
  delete temp2;
  
  int NBINS=DATASZ/(NumEv*NBINS_t);
  int NBINS_Q = NBINS%NBINS_x >0 ? (NBINS/NBINS_x) + 1 : (NBINS/NBINS_x);
  
  int NumEv_Corr=DATASZ/(NBINS*NBINS_t);  
  std::cout<<"Total number of events is "<<DATASZ<<endl;
  std::cout<<"Number of events per bin is "<<NumEv_Corr<<endl;

  delete ch1;
  File->Close();
  delete File;

  */
  
  std::vector<vector<double>> bins;
  std::vector<double> binst;

  double *bins_t;
  double *bins_Q;
  double *bins_x;
  int bin_number=1;
  TString Folder_old=Folder;
    
  binst=bin_C(2,NBINS_t,ch1);
  bins_t=binst.data();


  //Bin_View(Data, bins_t, ch1, NumEv);
  int k0=0, kN=NBINS_t, m0=0, mN=NBINS_Q, n0=0, nN;
  int in_range;
  if(bin!=0)
    {
      bin_number=bin;
      k0=(bin_number -1)/9 +1 - 1;
      kN=(bin_number -1)/9 +1;
      in_range = (bin_number -1)%9 +1; //Convert the bin number into a bin_number between 1 and 9
      m0=(in_range -1)/NBINS_x - 1;
      mN=(in_range -1)/NBINS_x;
      n0=(in_range-1)%NBINS_x - 1;

    }
    
  for(int k=k0;k<kN;k++)
    {
      bins=Grid_Bins(k, bins_t, NBINS_t, NumEv,NBINS_x, ch1);
      bins_Q=bins.at(0).data();
      for(int m=m0;m<mN;m++)
	{
	  bins_x=bins.at(m+1).data();

	  nN = (bin==0) ? min(NBINS_x,NBINS - m*NBINS_x) : (in_range-1)%NBINS_x;
	    
	  for(int n=n0;n<nN;n++)
	    {
	        TH1 *Orig;
		TH1 *Most;
		TH1 *Maxi;

	      TCut cut_bin = TCut(Form("bestCandidateFlag==1  && t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1]));
	      std::cout<<Form("\n\n t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1])<<endl;
	      boundaries.clear();
	      boundaries.push_back(bins_t[k]);
	      boundaries.push_back(bins_t[k+1]);
	      boundaries.push_back(bins_Q[m]);
	      boundaries.push_back(bins_Q[m+1]);
	      boundaries.push_back(bins_x[n]);
	      boundaries.push_back(bins_x[n+1]);
	      Folder = Folder_old + TString("bin_")+Form("%i/",bin_number);
	      gSystem->Exec(TString("mkdir -p ") + Folder);

	      //Training(cut + cut_bin, DVCS, Pi0, Vars);
	      //Training_vars(Data, DVCS, Pi0, cut + cut_bin);
	      //Add_BDT_var(cut + cut_bin, Data, TData, Vars);
	      //Add_BDT_var(cut + cut_bin, DVCS, TDVCS, Vars);
	      //Add_BDT_var(cut + cut_bin, Pi0, TPi0, Vars);
	      ////Add_BDT_var_float(cut + cut_bin, Pi0, TPi0);
	      //Explore(TData, TDVCS, cut + cut_bin);
	      //Filter(TData, cut + cut_bin + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)), TString("Data_NP_Theta_g_5.root"));

	      std::cout<<"\n Get Contamination BDT way"<<endl;
	      Get_Contamination(cut + cut_bin, BDT_value);

	      std::cout<<"\n Get Contamination Mostafa way"<<endl;
	      //Add_BDT_var(cut + cut_bin, "/home/munoz/Datasets/DVMP/P/Quality_Pi_as_DVCS_P.root", "Tested_1gamma.root", Vars);
	      //Add_BDT_var(cut + cut_bin, "/home/munoz/Datasets/DVMP/P/Quality_Pi_as_DVCS_P_1.root", "Tested_1gamma_1.root", Vars);
	      
	      Most=Get_Contamination_Mostafa(cut + cut_bin, BDT_value);
	      //boundaries.push_back(1);
	      //boundaries.push_back(1);

	      std::cout<<"\n Get Contamination Maxime way"<<endl;
	      Maxi=Maxime(cut + cut_bin, BDT_value, bin_number);

	      
	      Orig=Single_BSA("Data_NP_Theta_g_5.root",N_Phi);
	      Orig->SetTitle("Before");
	      //gStyle->SetOptFit(0);
	      gStyle->SetOptTitle(0);
	      TCanvas* c3 = new TCanvas("c3","Histograms");
	      Orig->Draw();
	      Most->Draw("SAME");
	      Maxi->Draw("SAME");
	      c3->BuildLegend();

	      c3->Print(Folder + TString("befo_vs_Most.pdf"));
	      gStyle->SetOptTitle(1);

	      Compare_three(cut + cut_bin, "TMaxime_pi0.root", "TMostafa_pi0.root");
	      delete c3;
	      delete Orig;
	      delete Most;
	      delete Maxi;
	      
	      bin_number+=1;

	      std::ofstream outFile(Folder + TString("Amplitudes.txt"));
	      outFile<<"type value fit"<<endl;	      
	      outFile<<"Raw "<<BSA_Amplitude<<" "<<BSA_Amplitude_fit<<endl;	      
	      outFile<<"Mostafa "<<BSA_Amplitude_most<<" "<<BSA_Amplitude_most_fit<<endl;	      
	      outFile<<"Maxime "<<BSA_Amplitude_maxi<<" "<<BSA_Amplitude_maxi_fit<<endl;	      
	      outFile.close();

	      std::cout<<"type value fit"<<endl;	      
	      std::cout<<"Raw "<<BSA_Amplitude<<" "<<BSA_Amplitude_fit<<endl;	      
	      std::cout<<"Mostafa "<<BSA_Amplitude_most<<" "<<BSA_Amplitude_most_fit<<endl;	      
	      std::cout<<"Maxime "<<BSA_Amplitude_maxi<<" "<<BSA_Amplitude_maxi_fit<<endl;	      
	    }
	}

      bins.clear();
    }

  Folder = Folder_old;

  delete bins_t;
  delete bins_Q;
  delete bins_x;
  return;
}
