void BDT::Systematic_phi_BDT(int Nbins=27, TString path="/home/munoz/Documents/RG-A:2/P/Pi0_Training/"){

  std::vector<vector<double>> bins;
  std::vector<double> binst;
  std::vector<double> aux;

  double bins_t[4] = {-13.908200, -0.781163, -0.339554, -0.000041};
  double *bins_Q;
  double *bins_x;
  int bin_number=1;
  TString Folder_old=Folder;
    
  int NBINS_t=3;
  int NBINS_x=3;
  int NBINS_Q=3;
  int NBINS = 10;

  
  //NEED TO GIVE BINS AS INPUT!
  for(int k=0;k<NBINS_t;k++)
    {
      aux.clear();
      if (k==0)
	{
	  aux.push_back(1.003650);
	  aux.push_back(2.250580);
	  aux.push_back(3.232600);
	  aux.push_back(11.189300);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.059511);
	  aux.push_back(0.158438);
	  aux.push_back(0.239953);
	  aux.push_back(0.418176);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.125293);
	  aux.push_back(0.246012);
	  aux.push_back(0.345544);
	  aux.push_back(0.508824);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.181055);
	  aux.push_back(0.346616);
	  aux.push_back(0.461646);
	  aux.push_back(0.775363);
	  bins.push_back(aux);
	  aux.clear();
	}
      if (k==1)
	{
	  aux.push_back(1.001100);
	  aux.push_back(2.000500);
	  aux.push_back(2.748080);
	  aux.push_back(9.674350);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.056738);
	  aux.push_back(0.119315);
	  aux.push_back(0.167430);
	  aux.push_back(0.390013);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.111259);
	  aux.push_back(0.165230);
	  aux.push_back(0.218732);
	  aux.push_back(0.467046);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.153305);
	  aux.push_back(0.228296);
	  aux.push_back(0.295716);
	  aux.push_back(0.689304);
	  bins.push_back(aux);
	  aux.clear();
	}
      if (k==2)
	{
	  aux.push_back(1.000510);
	  aux.push_back(1.823140);
	  aux.push_back(2.389770);
	  aux.push_back(8.928210);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.056515);
	  aux.push_back(0.101721);
	  aux.push_back(0.118623);
	  aux.push_back(0.366043);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.101679);
	  aux.push_back(0.133084);
	  aux.push_back(0.151032);
	  aux.push_back(0.430764);
	  bins.push_back(aux);
	  aux.clear();
	  
	  aux.push_back(0.133748);
	  aux.push_back(0.178059);
	  aux.push_back(0.215117);
	  aux.push_back(0.657785);
	  bins.push_back(aux);
	  aux.clear();
	}
      bins_Q=bins.at(0).data();
      for(int m=0;m<NBINS_Q;m++)
	{
	  bins_x=bins.at(m+1).data();
	  for(int n=0;n<min(NBINS_x,NBINS - m*NBINS_x);n++)
	    {
	      TH1 *Orig;
	      TH1 *Most;
	      TH1 *Maxi;

	      Folder = Folder_old + TString("bin_")+Form("%i/",bin_number);
	      gSystem->Exec(TString("mkdir -p ") + Folder);
		
	      gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("Tested_Quality_Data.root ") + Folder);
	      gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("Tested_1gamma.root ") + Folder);
	      gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("Tested_1gamma_1.root ") + Folder);
	      gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("Data_NP_Theta_g_5.root ") + Folder);
	      gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("TMostafa_pi0.root ") + Folder);
	      //gSystem->Exec(TString("cp ") + path + TString("bin_")+Form("%i/",bin_number) + TString("TMaxime_pi0.root ") + Folder);
	      gSystem->Exec(TString("cp -r ") + path + TString("bin_")+Form("%i/",bin_number) + TString("dataset ") + Folder);

	      TCut cut_bin = TCut(Form("bestCandidateFlag==1  && t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1]));
	      std::cout<<Form("\n\n t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1])<<endl;
	      boundaries.clear();
	      boundaries.push_back(bins_t[k]);
	      boundaries.push_back(bins_t[k+1]);
	      boundaries.push_back(bins_Q[m]);
	      boundaries.push_back(bins_Q[m+1]);
	      boundaries.push_back(bins_x[n]);
	      boundaries.push_back(bins_x[n+1]);

	      boundaries.push_back(1);
	      boundaries.push_back(1);
	      std::cout<<"\n Get Contamination Mostafa way"<<endl;
	      
	      Most=Get_Contamination_Mostafa(cut + cut_bin, BDT_value);
	      
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

	      delete c3;
	      delete Orig;
	      delete Most;
	      delete Maxi;
	      
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


	      gSystem->Exec(TString("rm ") + Folder + TString("Tested_Quality_Data.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("Tested_1gamma.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("Tested_1gamma_1.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("Data_NP_Theta_g_5.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("TMostafa_pi0.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("TMaxime_pi0.root"));
	      gSystem->Exec(TString("rm ") + Folder + TString("TMaxime_pi0.root"));
	      gSystem->Exec(TString("rm -r ") + Folder + TString("dataset"));

	      bin_number++;


	    }
	}
    }
  

}

