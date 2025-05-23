void BDT::Stack(TCut cut, TCut cut_add="")
{
  //TFile *input = new TFile("BDT/Tested_Quality_RGH.root","READ");
  TFile *input = new TFile(Folder + TString("Data_NP_Theta_g_5.root"),"READ");
  TTree *Data = (TTree*)input->Get("pDVCS");

  //TFile *input3 = new TFile("~/Documents/RG-A:DVCS/files/Sample_Train_Sim_DVCS.root","READ");
  //TTree *Sim = (TTree*)input3->Get("pDVCS");

  //ask for histograms to draw
  printf("Do you wanna generate the histograms for:\n");
  printf("All (0), selection of histogras (1) or give an input selection (2)\n");
  //&& mm2_eg<4.0 && cos2theta > 0.1 && p_perp < 0.3 && strip_Nuc_Theta > 5 && strip_Nuc_Theta < 65 
  int choose=0;
  int pass[33];
  int i1=1, i2=1, i3=1, i4=1, i5=1, i6=1, i7=1, i8=1;

  printf("Kinematic variables \n");
  pass[0]=i1;
  pass[1]=i1;
  pass[2]=i1;
  pass[3]=i1;
  printf("Electron kinematics \n");
  pass[4]=i2;
  pass[5]=i2;
  pass[6]=i2;
  pass[7]=i2;
  printf("Photon Kinematics \n");
  pass[8]=i3;
  pass[9]=i3;
  pass[10]=i3;
  pass[11]=i3;
  printf("Proton kinematics \n");
  pass[12]=i4;
  pass[13]=i4;
  pass[14]=i4;
  pass[15]=i4;
  printf("Missing masses \n");
  pass[16]=i5;
  pass[17]=i5;
  pass[18]=i5;
  pass[19]=i5;
  printf("Angles\n");
  pass[20]=i6;
  pass[21]=i6;
  pass[22]=i6;
  pass[23]=i6;
  printf("dt, dphi, dcos2theta, phi\n");
  pass[24]=i7;
  pass[25]=i7;
  pass[26]=i7;
  pass[27]=i7;
  printf("N_Ph, N_El, N_Nuc and missing momentum\n");
  pass[28]=i8;
  pass[29]=i8;
  pass[30]=i8;
  pass[31]=i8;
  pass[32]=i8;
 

  //ask for cutted or MC comparison
  //ask for 2x2 grid or single ones
  int grid=1;
  //Ask for printing
  int print=1;

  TCanvas *c[33];
  //THStack Collection[33] = new THStack("Collection","");
  TH1F *hist[33];
  TH1F *hist_cutted[33];

  TString variable[33]={"t_Ph",
			"strip_Q2",
			"strip_Xbj",
			"strip_W",
			"strip_El_P",
			"strip_El_Theta",
			"strip_El_px",
			"strip_El_E",
			"strip_Ph_P",
			"strip_Ph_Theta",
			"strip_Ph_px",
			"strip_Ph_E",
			"strip_Nuc_P",
			"strip_Nuc_Theta",
			"strip_Nuc_px",
			"strip_Nuc_E",
			"mm2_eNg",
			"mm2_eg",
			"mm2_g",
			"mm2_e",
			"theta_N_e",
			"theta_gamma_e",
			"theta_gamma_X",
			"cos2theta",
			"delta_t",
			"delta_Phi",
			"dcos2theta",
			"Phi_Ph",
			"N_Ph",
			"N_El",
			"N_Nuc",
			"miss_mom_eNg",
			"p_perp"
  };

  TString Title[33]={"t (Gev^{2})",
		     "Q^{2} (Gev^{2})",
		     "x_{B}",
		     "W (Gev)",
		     "Electron P (Gev)",
		     "#theta_{e} (deg)",
		     "Electron Px (Gev)",
		     "Electron E (Gev)",
		     "Photon P (Gev)",
		     "#theta_{#gamma} (deg)",
		     "Photon Px (Gev)",
		     "Photon E (Gev)",
		     "Proton P (Gev)",
		     "#theta_{p} (deg)",
		     "Proton Px (Gev)",
		     "Proton E (Gev)",
		     "M_{eNg}^{2} (Gev^{2})",
		     "M_{eg}^{2} (Gev^{2})",
		     "M_{g}^{2} (Gev^{2})",
		     "M_{e}^{2} (Gev^{2})",
		     "#theta_{Ne} (deg)",
		     "#theta_{#gamma e} (deg)",
		     "#theta_{#gamma X} (deg)",
		     "cos^{2}(#theta_{#gamma^{*}N})",
		     "#Delta t (Gev^{2})",
		     "#Delta #phi (deg)",
		     "#Delta cos^{2}#theta_{#gamma^{*}N}",
		     "#phi (deg)",
		     "Number of photons",
		     "Number of electrons",
		     "Number of protons",
		     "P^{miss}_{eNg} (GeV)",
		     "Missing transverse P eNg"	 
  };

  static double Low[33]={-10,0,0,0, 0,0,0,0,     0,0,0,0,     0,0,0,0,     -0.1,0,0,0,   0,0,0,-1,  -2,-2,-1,0, 0,0,0,    0,0};
  static double High[33]={0,10,1,5, 10,45,10,10, 10,45,10,10, 10,140,10,10, 0.1,5,17,17, 80,80,40,1, 2,2,1,360, 20,20,10, 1,1};

  int j=-1;

  for(int i=0;i<33;i++)
    {
      if(pass[i]==1)
	{
	  if(grid==1)
	    {
	      if(i%4==0)
		{
		  j++;
		  c[j]=new TCanvas(Form("c%i",i),"c");
		  c[j]->Divide(2,2);
		}
	    }
	  else
	    {
	      j=i;
	      c[j]=new TCanvas(Form("c%i",i),"c");
	    }
	  hist[i] = new TH1F (variable[i],Title[i],100,Low[i],High[i]);
	  hist_cutted[i] = new TH1F (variable[i]+"_Cutted",Title[i],100,Low[i],High[i]);
	  Data->Project(variable[i],variable[i],cut);
	  Data->Project(variable[i]+"_Cutted",variable[i],cut+cut_add);
	  hist[i]->SetMinimum(0);
	  hist_cutted[i]->SetMinimum(0);
	  hist_cutted[i]->SetLineColor(kRed);

	  if(grid==1) c[j]->cd(i%4 +1);
	  if(hist[i]->GetMaximum() > hist_cutted[i]->GetMaximum())
	    {
	      hist[i]->Draw();
	      hist_cutted[i]->Draw("SAME");
	    }
	  else
	    {
	      hist_cutted[i]->Draw();
	      hist[i]->Draw("SAME");
	    }
      	}      
    }
  if(print==1)
    {
      gSystem->Exec("mkdir -p " + Folder + "Graphs");
      gSystem->Exec("rm "+ Folder +"Graphs/*.pdf");
      for(int k=0; k<=j;k++)
	c[k]->Print(Folder + Form("Graphs/Canvas_%i.pdf",k));
    }
  std::cout<<"Contamination="<<hist_cutted[0]->GetEntries()*100/hist[0]->GetEntries()<<"%"<<endl;
}

