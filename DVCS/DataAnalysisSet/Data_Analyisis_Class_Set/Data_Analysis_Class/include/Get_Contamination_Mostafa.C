TH1* BDT::Get_Contamination_Mostafa(TCut cut, double BDT_cut, int Nphi, bool eta=false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain *Pi01g= new TChain("pDVCS");
  TChain *Pi02g;
  TChain *Pi02gData;
  TString output;
  TCut cut2g;      
  TString branch;

  //Convert cut to string
  TString String_cut = cut.GetTitle();

  if(eta)
  {
  std::cout<<"...:::COMPUTING ETA CONTAMINATION:::...\n"<<endl;
  Pi02g= new TChain("epeta");
  Pi02gData= new TChain("epeta");

  Pi01g->Add(Folder + TString("Tested_1gamma_eta.root"));
  Pi02g->Add(sim_epeta);
  Pi02gData->Add(epeta);

  output = "Mostafa_eta.root";
  branch="Phi_eta";
  
  //replace Ph by Pi0
  String_cut.ReplaceAll("Ph ", "eta ");
  String_cut.ReplaceAll("Ph>", "eta>");
  String_cut.ReplaceAll("Ph<", "eta<");
  String_cut.ReplaceAll("strip_Ph_P", "strip_W");
  String_cut.ReplaceAll("gamma", "eta");

  //Create cut for 2gamma case
  cut2g = TCut(String_cut) + TCut("mm2_egg<1.5 && strip_eta_2DChi2 < 0.04 ");
  }
  else
  {

  Pi02g= new TChain("eppi0");
  Pi02gData= new TChain("eppi0");

  Pi01g->Add(Folder + TString("Tested_1gamma.root"));
  //Pi01g->Add(Folder + TString("Tested_1gamma_1.root"));

  Pi02g->Add(sim_eppi0);
  Pi02g->Add(sim_eppi0_1);

  //Use quality 3 as it has the cuts for RGA eppi0 selection
  Pi02gData->Add(eppi0);
  branch="Phi_Pi0";

  output = "Mostafa_pi0.root";

  //replace Ph by Pi0
  String_cut.ReplaceAll("Ph ", "Pi0 ");
  String_cut.ReplaceAll("Ph>", "Pi0>");
  String_cut.ReplaceAll("Ph<", "Pi0<");
  String_cut.ReplaceAll("strip_Ph_P", "strip_W");
  String_cut.ReplaceAll("gamma", "Pi0");

  //Create cut for 2gamma case
  cut2g = TCut(String_cut) + TCut("mm2_egg<1.5 && strip_Pi0_2DChi2 < 0.04 ");
  }

  cut=cut + cut_ref;

  TChain *Data= new TChain("pDVCS");
  Data->Add(Folder + TData);

  TCut cut_phi;
  TCut cut_phi2;
  double est1;
  double est1_FT;
  double est1_FD;
  double est2;
  double est2_FT;
  double est2_FD;
  std::ofstream outFile;
  std::ofstream outFile2;
  
  if(generate_most)
    {
      TH1F *Pi01g_p = new TH1F("Pi01g_p","",Nphi,0,360);
      TH1F *Pi01g_m = new TH1F("Pi01g_m","",Nphi,0,360);
      TH1F *Pi01g_p_BDT = new TH1F("Pi01g_p_BDT","",Nphi,0,360);
      TH1F *Pi01g_m_BDT = new TH1F("Pi01g_m_BDT","",Nphi,0,360);
  
      TH1F *Pi02g_p = new TH1F("Pi02g_p","",Nphi,0,360);
      TH1F *Pi02g_m = new TH1F("Pi02g_m","",Nphi,0,360);
  
      TH1F *Pi02g_Data_p = new TH1F("Pi02g_Data_p","",Nphi,0,360);
      TH1F *Pi02g_Data_m = new TH1F("Pi02g_Data_m","",Nphi,0,360);

      TH1F *Data_p = new TH1F("Data_p","",Nphi,0,360);
      TH1F *Data_m = new TH1F("Data_m","",Nphi,0,360);
      TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",Nphi,0,360);
      TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",Nphi,0,360);

      TH1F *ratio_p = new TH1F("ratio_p","ratio_p",Nphi,0,360);
      TH1F *ratio_m = new TH1F("ratio_m","ratio_m",Nphi,0,360);
      TH1F *ratio_p_BDT = new TH1F("ratio_p_BDT","ratio_p_BDT",Nphi,0,360);
      TH1F *ratio_m_BDT = new TH1F("ratio_m_BDT","ratio_m_BDT",Nphi,0,360);

      TH1F *bkg_p = new TH1F("bkg_p","bkg_p",Nphi,0,360);
      TH1F *bkg_m = new TH1F("bkg_m","bkg_m",Nphi,0,360);
      TH1F *bkg_p_BDT = new TH1F("bkg_p_BDT","bkg_p_BDT",Nphi,0,360);
      TH1F *bkg_m_BDT = new TH1F("bkg_m_BDT","bkg_m_BDT",Nphi,0,360);

      TH1F *data_p_bkg_free = new TH1F("data_p_bkg_free","data_p_bkg_free",Nphi,0,360);
      TH1F *data_m_bkg_free = new TH1F("data_m_bkg_free","data_m_bkg_free",Nphi,0,360);

      TH1F *data_p_BDT_bkg_free = new TH1F("data_p_BDT_bkg_free","data_p_BDT_bkg_free",Nphi,0,360);
      TH1F *data_m_BDT_bkg_free= new TH1F("data_m_BDT_bkg_free","data_m_BDT_bkg_free",Nphi,0,360);


      Pi01g_p->Sumw2();
      Pi01g_m->Sumw2();
      Pi01g_p_BDT->Sumw2();
      Pi01g_m_BDT->Sumw2();
      Pi02g_p->Sumw2();
      Pi02g_m->Sumw2();
      Pi02g_Data_p->Sumw2();
      Pi02g_Data_m->Sumw2();
      Data_p->Sumw2();
      Data_m->Sumw2();

      Pi01g->Project("Pi01g_p", "Phi_Ph", cut);
      Pi01g->Project("Pi01g_m", "Phi_Ph", cut);
      Pi01g->Project("Pi01g_p_BDT", "Phi_Ph", cut + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
      Pi01g->Project("Pi01g_m_BDT", "Phi_Ph", cut + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));

      Pi02g->Project("Pi02g_p", branch, cut2g + TCut("abs(mm2_eNgg)<0.01"));
      Pi02g->Project("Pi02g_m", branch, cut2g + TCut("abs(mm2_eNgg)<0.01"));
  
      Pi02gData->Project("Pi02g_Data_p", branch, cut2g + TCut("Helicity>0") + TCut("abs(mm2_eNgg)<0.01"));
      Pi02gData->Project("Pi02g_Data_m", branch, cut2g + TCut("Helicity<0") + TCut("abs(mm2_eNgg)<0.01"));

      Data->Project("Data_p", "Phi_Ph", cut + TCut("Helicity>0"));
      Data->Project("Data_m", "Phi_Ph", cut + TCut("Helicity<0"));
      Data->Project("Data_p_BDT", "Phi_Ph", cut + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
      Data->Project("Data_m_BDT", "Phi_Ph", cut + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));

      std::ofstream outFileStat(Folder + TString("Bkg_Most_stats.txt"));
      for (int i=1; i<=Nphi; i++)
        {
          outFileStat<<Pi01g_p_BDT->GetBinContent(i)<<" "<<Pi01g_m_BDT->GetBinContent(i)<<" "<<Pi02g_p->GetBinContent(i)<<" "<<Pi02g_m->GetBinContent(i)<<" "<<Pi02g_Data_p->GetBinContent(i)<<" "<<Pi02g_Data_m->GetBinContent(i)<<" "<<Data_p_BDT->GetBinContent(i)<<" "<<Data_m_BDT->GetBinContent(i)<<endl;
        }	      
      outFileStat.close();

      ratio_p->Divide(Pi01g_p,Pi02g_p,1,1);
      ratio_m->Divide(Pi01g_m,Pi02g_m,1,1);
      ratio_p_BDT->Divide(Pi01g_p_BDT,Pi02g_p,1,1);
      ratio_m_BDT->Divide(Pi01g_m_BDT,Pi02g_m,1,1);

      bkg_p->Multiply(ratio_p,Pi02g_Data_p,1,1);
      bkg_m->Multiply(ratio_m,Pi02g_Data_m,1,1);
      bkg_p_BDT->Multiply(ratio_p_BDT,Pi02g_Data_p,1,1);
      bkg_m_BDT->Multiply(ratio_m_BDT,Pi02g_Data_m,1,1);

      //********************************************
      //output Mostafa tree
      //*******************************************

      static vector<double>* phi_M;
      static vector<int>* flag_M;
      static int hel;

      TFile* out_sim_p=new TFile(Folder + TString("Mostafa_pi0_p.root"), "RECREATE");
      Pi01g->SetBranchStatus("Helicity", 0);
      TTree *tree_p = Pi01g->CopyTree(cut);
      tree_p->SetMaxTreeSize(100000000000LL);
          
      TBranch* Wb_p = tree_p->Branch("Weight", &weight);
      TBranch* He_p = tree_p->Branch("Helicity", &hel);

      tree_p->SetBranchAddress("Phi_Ph",&phi_M);
      tree_p->SetBranchAddress("bestCandidateFlag",&flag_M);

      for(int t=0; t<tree_p->GetEntries();t++)
	{
	  tree_p->GetEntry(t);
	  for(int t1=0; t1<flag_M->size(); t1++)
	    {
	      if(flag_M->at(t1)==1)
		{
		  hel=1;
		  weight = bkg_p->GetBinContent(bkg_p->FindBin(phi_M->at(t1)))*1.0/Pi01g_p->GetBinContent(Pi01g_p->FindBin(phi_M->at(t1)));
		  Wb_p->Fill();
		  He_p->Fill();
		}
	    }
	}

      tree_p->Write();
      delete tree_p;
      out_sim_p->Close();
      delete out_sim_p;
  

      TFile* out_sim_m=new TFile(Folder + TString("Mostafa_pi0_m.root"), "RECREATE");
      Pi01g->SetBranchStatus("Helicity", 0);
      TTree *tree_m = Pi01g->CopyTree(cut);
      tree_m->SetMaxTreeSize(100000000000LL);
          
      TBranch* Wb_m = tree_m->Branch("Weight", &weight);
      TBranch* He_m = tree_m->Branch("Helicity", &hel);

      tree_m->SetBranchAddress("Phi_Ph",&phi_M);
      tree_m->SetBranchAddress("bestCandidateFlag",&flag_M);

      for(int t=0; t<tree_m->GetEntries();t++)
	{
	  tree_m->GetEntry(t);
	  for(int t1=0; t1<flag_M->size(); t1++)
	    {
	      if(flag_M->at(t1)==1)
		{
		  hel=-1;
		  weight = bkg_m->GetBinContent(bkg_m->FindBin(phi_M->at(t1)))*1.0/Pi01g_m->GetBinContent(Pi01g_m->FindBin(phi_M->at(t1)));
		  Wb_m->Fill();
		  He_m->Fill();
		}
	    }
	}

      tree_m->Write();
      delete tree_m;
      out_sim_m->Close();
      delete out_sim_m;
  
      TFile* outFileMost=new TFile(Folder + output, "RECREATE");
      TChain *outChain= new TChain("pDVCS");
      outChain->Add(Folder + TString("Mostafa_pi0_p.root"));
      outChain->Add(Folder + TString("Mostafa_pi0_m.root"));
      TTree *outTree = outChain->CopyTree(cut);

      outTree->Write();
      delete outTree;
      outFileMost->Close();
      delete outFileMost;
      delete outChain;
      gSystem->Exec(TString("rm ") + Folder + TString("Mostafa_pi0_p.root"));
      gSystem->Exec(TString("rm ") + Folder + TString("Mostafa_pi0_m.root"));



      delete Pi01g_p;
      delete Pi01g_m;
      delete Pi01g_p_BDT;
      delete Pi01g_m_BDT;

      delete Pi02g_p;
      delete Pi02g_m;

      delete Pi02g_Data_p;
      delete Pi02g_Data_m;

      delete Data_p;
      delete Data_m;
      delete Data_p_BDT;
      delete Data_m_BDT;

      delete ratio_p;
      delete ratio_m;
      delete ratio_p_BDT;
      delete ratio_m_BDT;

      delete bkg_p;
      delete bkg_m;
      delete bkg_p_BDT;
      delete bkg_m_BDT;

      delete data_p_bkg_free;
      delete data_m_bkg_free;
      delete data_p_BDT_bkg_free;
      delete data_m_BDT_bkg_free;  
  
      Add_BDT_var(cut, Folder + output, TString("T") + output, Vars);
      gSystem->Exec(TString("rm ") + Folder + output);

    }
  
  
  
  TChain *bkg= new TChain("pDVCS");
  bkg->Add(Folder + TString("T") + output);

  TString mvars[4] = {"Phi_Ph","t_Ph","strip_Q2","strip_Xbj"};
  double Lbound[4] = {0.  , boundaries.at(0), boundaries.at(2), boundaries.at(4)};
  double Ubound[4] = {360., boundaries.at(1), boundaries.at(3), boundaries.at(5)};

//Mean on the 3D bin
for(int j=1; j<4; j++)
  {
    TH1F *Data_BDT = new TH1F("Data_BDT","",1000,Lbound[j], Ubound[j]);
    Data_BDT->Sumw2();

    TH1F *Phi_BDT = new TH1F("Phi_BDT","",1000,Lbound[j], Ubound[j]);
    Phi_BDT->Sumw2();

    Data->Project("Data_BDT", mvars[j], cut + cut_phi + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
    bkg->Project("Phi_BDT", mvars[j], (cut + cut_phi + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));

    if(Data_BDT->Integral() < Phi_BDT->Integral())
      {
	for(int k=1; k<=Data_BDT->GetNbinsX(); k++)
	  Phi_BDT->SetBinContent(k, Data_BDT->GetBinContent(k));
      }

    Data_BDT->Add(Data_BDT, Phi_BDT,1 ,-1);
    std::cout<<mvars[j]<<" (Mostafa) mean is: "<<Data_BDT->GetMean()<<endl;
	  
    switch(j)
      {
      case 1:
	tmean1 = Data_BDT->GetMean();
	break;
      case 2:
	Qmean1 = Data_BDT->GetMean();
	break;
      case 3:
	xmean1 = Data_BDT->GetMean();
	break;
      }
    delete Phi_BDT;
    delete Data_BDT;

  }
//end of mean on the 3D bin

  if(means_most)
    {
      TString mvars[4] = {"Phi_Ph","t_Ph","strip_Q2","strip_Xbj"};
      double Lbound[4] = {0., -10., 0., 0.};
      double Ubound[4] = {360., 0., 10., 1.0};
      if(eta)
      {
      outFile.open(Folder + TString("means_most_eta.txt"));
      outFile2.open(Folder + TString("contamination_most_eta.txt"));
      }
      else
      {
      outFile.open(Folder + TString("means_most.txt"));
      outFile2.open(Folder + TString("contamination_most.txt"));
      }
      
      outFile <<"mean_phi error mean_t error mean_Q2 error mean_xB error"<<endl;
      outFile2 <<"bin_number entries_bef entries_aft before after"<<endl;

      for(int i=0; i<Nphi; i++)
	{
	  cut_phi =TCut(Form("Phi_Ph > %f && Phi_Ph < %f", (360./Nphi)*i, (360./Nphi)*(i+1)));
	  for(int j=0; j<4; j++)
	    {

	      TH1F *Data_p = new TH1F("Data_p","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_m = new TH1F("Data_m","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",1000,Lbound[j], Ubound[j]);

	      Data_p->Sumw2();
	      Data_m->Sumw2();

	      Data->Project("Data_p", mvars[j], cut + cut_phi + TCut("Helicity>0"));
	      Data->Project("Data_m", mvars[j], cut + cut_phi + TCut("Helicity<0"));
	      Data->Project("Data_p_BDT", mvars[j], cut + cut_phi + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
	      Data->Project("Data_m_BDT", mvars[j], cut + cut_phi + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));

	      //****************To get estimation before BDT*****************
	      TH1F *Phi_p = new TH1F("Phi_p","",1000,Lbound[j], Ubound[j]);
	      TH1F *Phi_m = new TH1F("Phi_m","",1000,Lbound[j], Ubound[j]);
	      Phi_p->Sumw2();
	      Phi_m->Sumw2();
	      
	      //No helicity cut because it is simulated data
	      bkg->Project("Phi_p", mvars[j], (cut + cut_phi + TCut("Helicity > 0"))*TCut("Weight"));
	      bkg->Project("Phi_m", mvars[j], (cut + cut_phi + TCut("Helicity < 0"))*TCut("Weight"));
	      
	      if(Data_p->Integral() < Phi_p->Integral())
	      {
	      for(int k=1; k<=Data_p->GetNbinsX(); k++)
		      Phi_p->SetBinContent(k, Data_p->GetBinContent(k));
		}

	      if(Data_m->Integral() < Phi_m->Integral())
	      {
	      for(int k=1; k<=Data_m->GetNbinsX(); k++)
		      Phi_m->SetBinContent(k, Data_m->GetBinContent(k));
		}

	      if(Data_p->Integral() + Data_m->Integral() ==0)
	      	est1=0;
	      else
	      	est1=(Phi_p->Integral() + Phi_m->Integral())*1.0/(Data_p->Integral() + Data_m->Integral());  

	      delete Phi_p;
	      delete Phi_m;
	      //*************************************************************
	      TH1F *Phi_p_BDT = new TH1F("Phi_p_BDT","",1000,Lbound[j], Ubound[j]);
	      TH1F *Phi_m_BDT = new TH1F("Phi_m_BDT","",1000,Lbound[j], Ubound[j]);

	      Phi_p_BDT->Sumw2();
	      Phi_m_BDT->Sumw2();

	      bkg->Project("Phi_p_BDT", mvars[j], (cut + cut_phi + TCut(Form("Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
	      bkg->Project("Phi_m_BDT", mvars[j], (cut + cut_phi + TCut(Form("Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
	      
	      if(Data_p_BDT->Integral() < Phi_p_BDT->Integral())
	      {
	      for(int k=1; k<=Data_p_BDT->GetNbinsX(); k++)
		      Phi_p_BDT->SetBinContent(k, Data_p_BDT->GetBinContent(k));
		}

	      if(Data_m_BDT->Integral() < Phi_m_BDT->Integral())
	      {
	      for(int k=1; k<=Data_m_BDT->GetNbinsX(); k++)
		      Phi_m_BDT->SetBinContent(k, Data_m_BDT->GetBinContent(k));
		}

	      entries_bef_most=Data_p_BDT->Integral() + Data_m_BDT->Integral() ;
	      if(Data_p_BDT->Integral() + Data_m_BDT->Integral() ==0)
	      	est2=0;
	      else
	      	est2=(Phi_p_BDT->Integral() + Phi_m_BDT->Integral())*1.0/(Data_p_BDT->Integral() + Data_m_BDT->Integral());

	      Data_p_BDT->Add(Data_p_BDT, Phi_p_BDT,1 ,-1);
	      Data_m_BDT->Add(Data_m_BDT, Phi_m_BDT,1 ,-1);
	      
	      entries_aft_most=Data_p_BDT->Integral() + Data_m_BDT->Integral();

	      TH1F *mean_final= new TH1F("mean_final","mean_final",1000, Lbound[j], Ubound[j]);
	      mean_final->Add(Data_p_BDT, Data_m_BDT, 1, 1);

  
	  
	      if(j==0)
		{
		  std::cout<<"\n Mostafa estimation: Phi bin "<<i+1<<" : "<<est1*100<<"% "<<est2*100<<"%"<<endl;
		  outFile2<<i+1<<" "<<entries_bef_most<<" "<<entries_aft_most<<" "<<est1*100<<"% "<<est2*100<<"%"<<endl;

		  if(entries_aft_most>1)
		  	est1=mean_final->GetMean();
		  else
		  	est1=0;

		  std::cout<<mvars[j]<<" mean in bin "<<i+1<<" is: "<<est1<<endl;
		  outFile<<" "<<est1<<" "<<mean_final->GetStdDev();
		}
	      else
		{
		  if(entries_aft_most>1)
		  	est1=mean_final->GetMean();
		  else
		  	est1=0;

		  std::cout<<mvars[j]<<" mean in bin "<<i+1<<" is: "<<est1<<endl;
		  outFile<<" "<<est1<<" "<<mean_final->GetStdDev();
		}
	  
	      delete Phi_p_BDT;
	      delete Phi_m_BDT;
	      delete Data_p;
	      delete Data_m;
	      delete Data_p_BDT;
	      delete Data_m_BDT;
	      delete mean_final;

	    }
	  outFile<<" "<<endl;
	}
      outFile.close();
    }






  TH1F *Data_p = new TH1F("Data_p","",Nphi,0,360);
  TH1F *Data_m = new TH1F("Data_m","",Nphi,0,360);
  TH1F *Data_p_FT = new TH1F("Data_p_FT","",Nphi,0,360);
  TH1F *Data_m_FT = new TH1F("Data_m_FT","",Nphi,0,360);
  TH1F *Data_p_FD = new TH1F("Data_p_FD","",Nphi,0,360);
  TH1F *Data_m_FD = new TH1F("Data_m_FD","",Nphi,0,360);
  TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",Nphi,0,360);
  TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",Nphi,0,360);
  TH1F *Data_p_BDT_FT = new TH1F("Data_p_BDT_FT","",Nphi,0,360);
  TH1F *Data_m_BDT_FT = new TH1F("Data_m_BDT_FT","",Nphi,0,360);
  TH1F *Data_p_BDT_FD = new TH1F("Data_p_BDT_FD","",Nphi,0,360);
  TH1F *Data_m_BDT_FD = new TH1F("Data_m_BDT_FD","",Nphi,0,360);

  Data_p->Sumw2();
  Data_m->Sumw2();

  Data->Project("Data_p", "Phi_Ph", cut + TCut("Helicity>0"));
  Data->Project("Data_m", "Phi_Ph", cut + TCut("Helicity<0"));
  Data->Project("Data_p_FT", "Phi_Ph", cut + TCut("strip_Ph_Theta < 5 && Helicity>0"));
  Data->Project("Data_m_FT", "Phi_Ph", cut + TCut("strip_Ph_Theta < 5 && Helicity<0"));
  Data->Project("Data_p_FD", "Phi_Ph", cut + TCut("strip_Ph_Theta > 5 && Helicity>0"));
  Data->Project("Data_m_FD", "Phi_Ph", cut + TCut("strip_Ph_Theta > 5 && Helicity<0"));
  Data->Project("Data_p_BDT", "Phi_Ph", cut + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT", "Phi_Ph", cut + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_p_BDT_FT", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta < 5 && Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT_FT", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta < 5 && Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_p_BDT_FD", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta > 5 && Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT_FD", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta > 5 && Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));

  //****************To get estimation before BDT*****************

  TH1F *Phi_p = new TH1F("Phi_p","",Nphi,0,360);
  TH1F *Phi_m = new TH1F("Phi_m","",Nphi,0,360);
  TH1F *Phi_p_FT = new TH1F("Phi_p_FT","",Nphi,0,360);
  TH1F *Phi_m_FT = new TH1F("Phi_m_FT","",Nphi,0,360);
  TH1F *Phi_p_FD = new TH1F("Phi_p_FD","",Nphi,0,360);
  TH1F *Phi_m_FD = new TH1F("Phi_m_FD","",Nphi,0,360);

  Phi_p->Sumw2();
  Phi_m->Sumw2();
  Phi_p_FT->Sumw2();
  Phi_m_FT->Sumw2();
  Phi_p_FD->Sumw2();
  Phi_m_FD->Sumw2();

  bkg->Project("Phi_p", "Phi_Ph", (cut + TCut("Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m", "Phi_Ph", (cut + TCut("Helicity < 0"))*TCut("Weight"));
  bkg->Project("Phi_p_FT", "Phi_Ph", (cut + TCut("strip_Ph_Theta < 5 && Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m_FT", "Phi_Ph", (cut + TCut("strip_Ph_Theta < 5 && Helicity < 0"))*TCut("Weight"));
  bkg->Project("Phi_p_FD", "Phi_Ph", (cut + TCut("strip_Ph_Theta > 5 && Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m_FD", "Phi_Ph", (cut + TCut("strip_Ph_Theta > 5 && Helicity < 0"))*TCut("Weight"));

  for(int k=1; k<=Data_p->GetNbinsX(); k++)
    {
      if(Data_p->GetBinContent(k) < Phi_p->GetBinContent(k))
	{
	  Phi_p->SetBinContent(k, Data_p->GetBinContent(k));
	    }
      if(Data_m->GetBinContent(k) < Phi_m->GetBinContent(k))
	{
	  Phi_m->SetBinContent(k, Data_m->GetBinContent(k));
	    }
      if(Data_p_FT->GetBinContent(k) < Phi_p_FT->GetBinContent(k))
	{
	  Phi_p_FT->SetBinContent(k, Data_p_FT->GetBinContent(k));
	    }
      if(Data_m_FT->GetBinContent(k) < Phi_m_FT->GetBinContent(k))
	{
	  Phi_m_FT->SetBinContent(k, Data_m_FT->GetBinContent(k));
	    }
      if(Data_p_FD->GetBinContent(k) < Phi_p_FD->GetBinContent(k))
	{
	  Phi_p_FD->SetBinContent(k, Data_p_FD->GetBinContent(k));
	    }
      if(Data_m_FD->GetBinContent(k) < Phi_m_FD->GetBinContent(k))
	{
	  Phi_m_FD->SetBinContent(k, Data_m_FD->GetBinContent(k));
	    }
    }

  if(Data_p->Integral() + Data_m->Integral() ==0)
 	est1=0;
  else
  	est1=(Phi_p->Integral() + Phi_m->Integral())*1.0/(Data_p->Integral() + Data_m->Integral());  

  if(Data_p_FT->Integral() + Data_m_FT->Integral() ==0)
 	est1_FT=0;
  else
  	est1_FT=(Phi_p_FT->Integral() + Phi_m_FT->Integral())*1.0/(Data_p_FT->Integral() + Data_m_FT->Integral());  

  if(Data_p_FD->Integral() + Data_m_FD->Integral() ==0)
 	est1_FD=0;
  else
  	est1_FD=(Phi_p_FD->Integral() + Phi_m_FD->Integral())*1.0/(Data_p_FD->Integral() + Data_m_FD->Integral());  


  delete Phi_p;
  delete Phi_m;
  delete Phi_p_FT;
  delete Phi_m_FT;
  delete Phi_p_FD;
  delete Phi_m_FD;
  delete Data_p_FT;
  delete Data_m_FT;
  delete Data_p_FD;
  delete Data_m_FD;
  //*************************************************************
  TH1F *Phi_p_BDT = new TH1F("Phi_p_BDT","",Nphi,0,360);
  TH1F *Phi_m_BDT = new TH1F("Phi_m_BDT","",Nphi,0,360);
  TH1F *Phi_p_BDT_FT = new TH1F("Phi_p_BDT_FT","",Nphi,0,360);
  TH1F *Phi_m_BDT_FT = new TH1F("Phi_m_BDT_FT","",Nphi,0,360);
  TH1F *Phi_p_BDT_FD = new TH1F("Phi_p_BDT_FD","",Nphi,0,360);
  TH1F *Phi_m_BDT_FD = new TH1F("Phi_m_BDT_FD","",Nphi,0,360);

  Phi_p_BDT->Sumw2();
  Phi_m_BDT->Sumw2();

  bkg->Project("Phi_p_BDT", "Phi_Ph", (cut + TCut(Form("Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT", "Phi_Ph", (cut + TCut(Form("Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_p_BDT_FT", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta < 5 && Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT_FT", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta < 5 && Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_p_BDT_FD", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta > 5 && Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT_FD", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta > 5 && Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));

  for(int k=1; k<=Data_p_BDT->GetNbinsX(); k++)
    {
      if(Data_p_BDT->GetBinContent(k) < Phi_p_BDT->GetBinContent(k))
	{
	  Phi_p_BDT->SetBinContent(k, Data_p_BDT->GetBinContent(k));
	    }
      if(Data_m_BDT->GetBinContent(k) < Phi_m_BDT->GetBinContent(k))
	{
	  Phi_m_BDT->SetBinContent(k, Data_m_BDT->GetBinContent(k));
	    }
      if(Data_p_BDT_FT->GetBinContent(k) < Phi_p_BDT_FT->GetBinContent(k))
	{
	  Phi_p_BDT_FT->SetBinContent(k, Data_p_BDT_FT->GetBinContent(k));
	    }
      if(Data_m_BDT_FT->GetBinContent(k) < Phi_m_BDT_FT->GetBinContent(k))
	{
	  Phi_m_BDT_FT->SetBinContent(k, Data_m_BDT_FT->GetBinContent(k));
	    }
      if(Data_p_BDT_FD->GetBinContent(k) < Phi_p_BDT_FD->GetBinContent(k))
	{
	  Phi_p_BDT_FD->SetBinContent(k, Data_p_BDT_FD->GetBinContent(k));
	    }
      if(Data_m_BDT_FD->GetBinContent(k) < Phi_m_BDT_FD->GetBinContent(k))
	{
	  Phi_m_BDT_FD->SetBinContent(k, Data_m_BDT_FD->GetBinContent(k));
	    }
    }

  entries_bef_most=Data_p_BDT->Integral() + Data_m_BDT->Integral();
  entries_bef_most_FT=Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral();
  entries_bef_most_FD=Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral();
  
  if(Data_p_BDT->Integral() + Data_m_BDT->Integral() ==0)
  	est2=0;
  else
  	est2=(Phi_p_BDT->Integral() + Phi_m_BDT->Integral())*1.0/(Data_p_BDT->Integral() + Data_m_BDT->Integral());

  if(Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral() ==0)
  	est2_FT=0;
  else
  	est2_FT=(Phi_p_BDT_FT->Integral() + Phi_m_BDT_FT->Integral())*1.0/(Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral());

  if(Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral() ==0)
  	est2_FD=0;
  else
  	est2_FD=(Phi_p_BDT_FD->Integral() + Phi_m_BDT_FD->Integral())*1.0/(Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral());
  std::cout<<"Mostafa estimation "<<est1<<" "<<est2<<endl;
  std::cout<<"Mostafa estimation FT "<<est1_FT<<" "<<est2_FT<<endl;
  std::cout<<"Mostafa estimation FD "<<est1_FD<<" "<<est2_FD<<endl;
  
  Data_p_BDT->Add(Data_p_BDT, Phi_p_BDT,1 ,-1);
  Data_m_BDT->Add(Data_m_BDT, Phi_m_BDT,1 ,-1);
  
  entries_aft_most=Data_p_BDT->Integral() + Data_m_BDT->Integral();
  entries_aft_most_FT=Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral();
  entries_aft_most_FD=Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral();

  std::cout<<"Overall "<<entries_bef_most<<" "<<entries_aft_most<<" "<<entries_bef_most_FT<<" "<<entries_bef_most_FD<<" (Total "<<entries_bef_most_FT+entries_bef_most_FD<<") "<<est1*100<<"% "<<est2*100<<"% "<<est1_FT*100<<"% "<<est2_FT*100<<"% "<<est1_FD*100<<"% "<<est2_FD*100<<"% "<<endl;
  boundaries.push_back(est1);
  boundaries.push_back(est2);
  boundaries.push_back(est1_FT);
  boundaries.push_back(est2_FT);
  boundaries.push_back(est1_FD);
  boundaries.push_back(est2_FD);

  if(means_most)
    {
      outFile2<<"Overall "<<entries_bef_most<<" "<<entries_aft_most<<" "<<est1*100<<"% "<<est2*100<<"%"<<endl;
      outFile2.close();
    }
  
  TH1 *BA= Data_m_BDT->GetAsymmetry(Data_p_BDT);
  TCanvas* c2 = new TCanvas("c2","Histograms");

  
  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf = new TF1("fitf","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf->SetParameter(0,0.1);
  fitf->SetParameter(1,-0.3);
  fitf->SetParLimits(0,0,1.0);
  fitf->SetParLimits(1,-1.,1.);
	  
  BA->SetAxisRange(-1., 1.,"Y");	  
  BA->Scale(1.0/Bpol);
  //BA->Fit("fitf","Q");
  BA->SetTitle("RG-A note");


  // Create an output file to save the histogram
  if(eta)
  {
   std::ofstream outFile_eta(Folder + TString("eta_contamination.txt"));
   outFile_eta<<"eta estimation "<< " before/after: "<<entries_bef_most<<" "<<entries_aft_most<<" "<<est1*100<<"% "<<est2*100<<"%"<<endl;  
   outFile_eta.close();   
  }
  else
  {
  TFile *outputFile = new TFile(Folder + TString("Mostafa_Clean.root"), "RECREATE");
  fitf->SetLineColor(kBlack);
  BA->SetLineColor(kBlack);
  BA->SetMarkerColor(kBlack);

  // Write the histogram to the output file
  BA->Write();

  // Close the output file
  outputFile->Close();

  std::ofstream outFile5(Folder + TString("BSA_Most_Values.txt"));
  std::ofstream outFile6(Folder + TString("entries_most_p.txt"));
  std::ofstream outFile7(Folder + TString("entries_most_m.txt"));
  for(int k=1; k<=BA->GetNbinsX(); k++)
    {
      outFile5<<Data_p_BDT->GetBinContent(k) + Data_m_BDT->GetBinContent(k)<<", "<<BA->GetBinContent(k)<<", "<<BA->GetBinError(k)<<endl;
      outFile6<<Data_p_BDT->GetBinContent(k)<<endl;
      outFile7<<Data_m_BDT->GetBinContent(k)<<endl;
    }
  outFile5.close();
  outFile6.close();
  outFile7.close();

  delete outputFile;
}



  fitf->SetLineColor(kRed);
  BA->SetLineColor(kRed);
  BA->SetMarkerColor(kRed);
  BA->SetAxisRange(-1.0, 1.0,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi(deg)");
  BA->GetYaxis()->SetTitle("BSA");
  BA->SetTitle("Method 1 estimation");
  BA->GetXaxis()->SetTitleSize(0.06);
  BA->GetYaxis()->SetTitleSize(0.06);
  BA->GetXaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetTitleOffset(0.5);
  BA->GetXaxis()->SetTitleOffset(0.5);
  BA->GetYaxis()->SetNdivisions(6);
  BA->GetXaxis()->SetNdivisions(4);

  BA->Draw();

  if(!eta)
  {
  c2->Print(Folder + TString("Mostafa_BSA.pdf"));
  }
  
  BSA_Amplitude_most=BA->GetBinContent(BA->GetMaximumBin());
  BSA_Amplitude_most_fit=fitf->GetParameter(0);
  BSA_Error_most_fit=fitf->GetParError(0);

  std::cout<<"Final amplitudes: value - fit "<<BSA_Amplitude_most<<" "<<BSA_Amplitude_most_fit<<endl;



  // Clean up memory
      
  delete c2;
  delete Phi_p_BDT;
  delete Phi_m_BDT;
  delete Phi_p_BDT_FT;
  delete Phi_m_BDT_FT;
  delete Phi_p_BDT_FD;
  delete Phi_m_BDT_FD;
  delete Data_p;
  delete Data_m;
  delete Data_p_BDT;
  delete Data_m_BDT;
  delete Data_p_BDT_FT;
  delete Data_m_BDT_FT;
  delete Data_p_BDT_FD;
  delete Data_m_BDT_FD;

  delete Data;
  delete bkg;
  
  delete Pi01g;
  delete Pi02g;
  delete Pi02gData;
  return BA;
}
