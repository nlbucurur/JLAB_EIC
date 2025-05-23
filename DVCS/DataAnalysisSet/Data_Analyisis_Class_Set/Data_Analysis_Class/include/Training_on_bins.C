
void BDT::Training_on_bins(TString Data, int NBinsPhi=0, int bin=0, bool build=true, bool eta=false)//, TCut cut, TString DVCS, TString Pi0, vector<TString> vars)
{
  std::ifstream inputFile;
  std::string line;
  std::vector<double> xValues;
  std::vector<double> yValues;
  std::vector<double> yErrors;
  std::vector<double> xErrors;
  std::vector<double> phi_bins;

  double x, y, yErr;
  double binWidth;
  double xErr;
  int numPoints;

  
  int bin_number=1;
  TString Folder_old=Folder;
    
  int k0=0, kN=64;
  if(bin!=0)
    {
      bin_number=bin;
      k0=bin-1;
      kN=bin;      
    }
    

  for(int k=k0;k<kN;k++)
    {

      if(NBinsPhi==0)
	{
	  NBinsPhi=Nphibins[k];
	}

      cut_bin = Mbins[k];
      std::cout<<Form("\n\n t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins[k][0],bins[k][1],bins[k][2],bins[k][3],bins[k][4],bins[k][5])<<endl;
      boundaries.clear();
      boundaries.push_back(bins[k][0]);
      boundaries.push_back(bins[k][1]);
      boundaries.push_back(bins[k][2]);
      boundaries.push_back(bins[k][3]);
      boundaries.push_back(bins[k][4]);
      boundaries.push_back(bins[k][5]);
      Folder = Folder_old + Form("bin_%i/",bin_number);
      gSystem->Exec(TString("mkdir -p ") + Folder);

      //Flag for training category or normal      
      //Check if there is enough background in the FT. At least 100 events
     TChain *background= new TChain("pDVCS");
     background->Add(Pi0);
     
      TH1F *hist1 = new TH1F("hist1","",100,0,1);
      TH1F *hist2 = new TH1F("hist2","",100,0,1);
      background->Project("hist1","strip_Xbj",cut + cut_bin + TCut("strip_Ph_Theta < 5"));
      background->Project("hist2","strip_Xbj",cut + cut_bin);
      nft=hist1->GetEntries()/hist2->GetEntries();
      categories=(nft>0.05 && hist1->GetEntries()>50) ? true: false ;
      int nev_bkg = hist2->GetEntries();
      delete hist1;
      delete hist2;
      delete background;
      //categories=false;          

      TH1 *Orig;
      TH1 *Most;
      TH1 *Maxi;

      std::ifstream nodataInp(Folder + TString("nodata.conf"));
      if (nodataInp) nodata = true;
      nodataInp.close();
      std::ofstream nodataFile;	
  
      if(build)
	{
    std::cout<<"\nBuilding ML training..."<<endl;	

    try 
    {
      if(nev_bkg<50)
        throw std::runtime_error("Manually triggered exception!");

      TrainingCategory(cut + cut_bin, DVCS, Pi0, Vars);
    } 
    catch (std::exception& e) 
    {
      std::cerr << "NO DATA \n Impossible to do a BDT classification. \n generating contamination and exiting..." << "\n";
      nodata=true;
      Write_Null(bin_number);

      gSystem->ChangeDirectory(dir);
      gSystem->Exec("pwd");

      nodataFile.open(Folder + TString("nodata.conf"));
      if (!nodataFile.is_open()) std::cerr << "Error opening nodata.conf :" <<Folder + TString("nodata.conf")<< std::endl;
      nodataFile<<"true"<<endl;
      nodataFile.close();

	/*
      if(generate || recast)
      {
      	std::cout<<"\n Get Contamination Maxime way"<<endl;
      	Maxi=Maxime(cut + cut_bin + cut_ref, BDT_value, bin_number,NBinsPhi);
      }
      return;
      */
    }

	  if(pid_sys)
	  	Add_BDT_var(cut + cut_bin + cut_pid, Data, TData, Vars);
	  else
	  	Add_BDT_var(cut + cut_bin, Data, TData, Vars);
	  	
	  Add_BDT_var(cut + cut_bin, DVCS, TDVCS, Vars);
	  Add_BDT_var(cut + cut_bin, Pi0, TPi0, Vars);
	  Training_vars(Folder + TData, Folder + TDVCS, Folder + TPi0, cut + cut_bin);
	  
	  if(!nodata)
	  	Get_BDT_Score();

    ////Add_BDT_var_float(cut + cut_bin, Pi0, TPi0);
	  Explore(TData, TDVCS, cut + cut_bin);

	  Add_BDT_var(cut + cut_bin, sim_pi_as_dvcs, "Tested_1gamma.root", Vars);
	  Add_BDT_var(cut + cut_bin, sim_pi_as_dvcs_1, "Tested_1gamma_1.root", Vars);
	  Add_BDT_var(cut + cut_bin, sim_eta_as_dvcs, "Tested_1gamma_eta.root", Vars);
    Add_BDT_var(cut + cut_bin, BM_Sim, TBM_Sim, Vars);
	  
	  Filter_Pi0(eppi0, l.cut + cut_bin, eppi0_name);	  
	  Filter_Pi0(sim_eppi0, l.cut + cut_bin, TString("sim_")+eppi0_name);	  
          eppi0_hists();	  
	  ////Resolution_Match(cut2g, eppi0_name, TString("sim_")+eppi0_name);

          Kin_vars(TData, TDVCS, TPi0, cut + cut_bin);	
          gSystem->Exec(TString("mv ") + Folder + TString("Kin_Vars.pdf ") + Folder + TString("Kin_Vars_bef_BDT.pdf"));
          Kin_vars(TData, TDVCS, TPi0, cut + cut_bin + TCut("_strip_Nuc_BDT > 0"));	
          gSystem->Exec(TString("mv ") + Folder + TString("Kin_Vars.pdf ") + Folder + TString("Kin_Vars_aft_BDT.pdf"));
	}
      //Refinement cuts
      //To be used at the end if needed
      //cut_ref=Rcuts[k];
/*
      Filter(TData, cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)), TString("Data_NP_Theta_g_5.root"));
      Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TData , TData, TString("pass1_comparison.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder);	  
      //Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TDVCS , TDVCS, TString("pass1_comparison_sim_DVCS.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder,1, true);	
      Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TPi0 , TPi0, TString("pass1_comparison_sim_Pi0.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder,1, true);	
      
      Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TData , TData, TString("pass1_comparison_1.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder,2);	  
      //Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TDVCS , TDVCS, TString("pass1_comparison_sim_DVCS_1.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder,2, true);	
      Compare_two(cut + cut_bin + cut_ref, cut + cut_bin + cut_ref, TPi0 , TPi0, TString("pass1_comparison_sim_Pi0_1.pdf"), TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/"),Folder,2, true);	
*/  	
      /*******************************************************************/	      
    //Add_BDT_var(cut + cut_bin, BM_Sim, TBM_Sim, Vars);
          Add_BDT_var(cut + cut_bin, sim_eta_as_dvcs, "Tested_1gamma_eta.root", Vars);

      if(generate || recast)
      {
      	std::cout<<"\n Get Contamination Maxime way"<<endl;
      	Maxi=Maxime(cut + cut_bin + cut_ref, BDT_value, bin_number,NBinsPhi);
      	return;
      }

      std::cout<<"\n Getting Contamination BDT way..."<<endl;
      Get_Contamination(cut + cut_bin + cut_ref, BDT_value);

      std::cout<<"\n Getting Contamination Mostafa way..."<<endl;      
      Most=Get_Contamination_Mostafa(cut + cut_bin, BDT_value,NBinsPhi, eta); //Cut ref is coded inside
      //boundaries.push_back(1);
      //boundaries.push_back(1);
      Most->GetXaxis()->SetLabelSize(0.05);
      Most->GetYaxis()->SetLabelSize(0.05);
      Most->GetXaxis()->SetTitleSize(0.05);
      Most->GetYaxis()->SetTitleSize(0.05);
      Most->GetXaxis()->SetTitleOffset(0.8);

      TFile *inputFile2 = TFile::Open(TString("/work/clas12/jsalvg/RGA-Analysis/inb/P/Analysis/") +  TString("bin_")+Form("%i/",bin_number)+  TString("Mostafa_Clean.root"));
      TH1 *pass1=dynamic_cast<TH1F*>(inputFile2->Get("Asymmetry_Data_m_BDT-Data_p_BDT"));
      pass1->SetTitle("Pass1");
      pass1->SetLineColor(kRed);
      pass1->SetMarkerColor(kRed);

      std::cout<<"\n Get Contamination Maxime way"<<endl;
      Maxi=Maxime(cut + cut_bin + cut_ref, BDT_value, bin_number,NBinsPhi);
      Maxi->GetXaxis()->SetLabelSize(0.05);
      Maxi->GetYaxis()->SetLabelSize(0.05);
      Maxi->GetXaxis()->SetTitleSize(0.05);
      Maxi->GetYaxis()->SetTitleSize(0.05);
      Maxi->GetXaxis()->SetTitleOffset(0.8);

      /*******************************************************************/	      
      /*
      std::cout<<"\nFiltering signal events..."<<endl;	      
      Filter(TData, cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)), TString("Data_NP_Theta_g_5.root"));
      Filter("TMostafa_pi0.root", cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)), TString("Filtered_Most.root"));
      Filter("TMaxime_pi0.root", cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_value)), TString("Filtered_Maxi.root"));
      */

      /*******************************************************************/	      
      std::cout<<"\nGetting comparison plots..."<<endl;	      
      Orig=Single_BSA("Data_NP_Theta_g_5.root",NBinsPhi);

      gStyle->SetOptFit(0);
      gStyle->SetOptTitle(0);
      Orig->SetAxisRange(-1., 1.,"Y");
      Orig->SetMarkerColor(kBlack);
      Orig->SetLineColor(kBlack);
      Most->SetAxisRange(-1., 1.,"Y");
      Most->SetMarkerColor(kBlue);
      Most->SetLineColor(kBlue);
      if(Most->GetFunction("fitf"))
      	Most->GetFunction("fitf")->SetLineColor(kWhite);
      Maxi->SetAxisRange(-1., 1.,"Y");
      Maxi->SetMarkerColor(kRed);
      Maxi->SetLineColor(kRed);
      if(Maxi->GetFunction("fitf"))
      	Maxi->GetFunction("fitf")->SetLineColor(kWhite);

      Most->SetTitle("Method 1");
      Maxi->SetTitle("Method 2");
	      
      TCanvas* c3 = new TCanvas("c3","Histograms");
      Orig->SetTitle("Before");
      Orig->Draw();
      Most->Draw("SAME");
      Maxi->Draw("SAME");
      c3->BuildLegend();
      c3->Print(Folder + TString("Background_subtraction.pdf"));
      delete c3;

      
      TCanvas* c4 = new TCanvas("c4","Histograms");
      gStyle->SetOptFit(0);
      gStyle->SetOptTitle(0);
      //Plot Maxime BSA
      inputFile.open(extBSA + Form("bin_%i.txt",bin_number)); // Replace with the name of the block file you want to read
      xValues.clear();
      yValues.clear();
      xErrors.clear();
      yErrors.clear();
      while (std::getline(inputFile, line)) {
	    std::istringstream iss(line);
	    iss >> x >> y >> yErr;
	    x=x*180/TMath::Pi();
	    xErr = (360./NBinsPhi) / 2.0;
        
	    xValues.push_back(x);
	    yValues.push_back(y);
	    yErrors.push_back(yErr);
	    xErrors.push_back(xErr);
      }
      numPoints = xValues.size();

      inputFile.close();
      TGraphErrors* graph = new TGraphErrors(numPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
      graph->SetMarkerColor(kBlack);
      graph->SetLineColor(kBlack);
      graph->SetMinimum(-1.);
      graph->SetMaximum(1.);
      graph->GetXaxis()->SetLimits(0,360);
      graph->SetTitle("RG-A Published");

	//To not display null measurements i.e. without events
      for(int t=1; t<=NBinsPhi+1;t++) //Including the overflow bin
      {
      if(Most->GetBinContent(t)==0 || Most->GetBinError(t)==0) Most->SetBinContent(t,10);
      if(Maxi->GetBinContent(t)==0 || Maxi->GetBinError(t)==0) Maxi->SetBinContent(t,10);
      }
      Most->Draw("E0");
      Maxi->Draw("SAME, E0");
      graph->Draw("SAME P");

      c4->BuildLegend();
      c4->Print(Folder + TString("RGA_vs_This.pdf"));
      delete c4;


      TCanvas* c5 = new TCanvas("c5","Histograms");
      Most->SetTitle("Pass2");
      Most->SetLineWidth(2);
      Most->Draw("SAME");
      pass1->Draw("SAME");
      //graph->Draw("SAME P");
      c5->BuildLegend();
      c5->Print(Folder + TString("pass1_BSA_comparison.pdf"));
      delete c5;
      Most->SetTitle("Method 1");
      Most->SetLineWidth(1);

      if(!nodata)
      {
      Single_BSA_2("../Data_NP_Theta_g_5.root","Data_NP_Theta_g_5.root", boundaries, NBinsPhi);
      
      gStyle->SetOptTitle(1);
      Get_excl_vars(cut + cut_bin + cut_ref);
      Compare_three(cut + cut_bin + cut_ref + TCut(Form("_strip_Nuc_BDT>%f",0.0)), "TMaxime_pi0.root", "TMostafa_pi0.root");
      }  
      /*******************************************************************/	      
      std::cout<<"\nSUMMARY:"<<endl;	      

      std::ofstream outFile(Folder + TString("Amplitudes.txt"));
      outFile<<"type value fit error on fit"<<endl;	      
      outFile<<"Raw "<<BSA_Amplitude<<" "<<BSA_Amplitude_fit<<" "<<BSA_Error_fit<<endl;	      
      outFile<<"Mostafa "<<BSA_Amplitude_most<<" "<<BSA_Amplitude_most_fit<<" "<<BSA_Error_most_fit<<endl;	      
      outFile<<"Maxime "<<BSA_Amplitude_maxi<<" "<<BSA_Amplitude_maxi_fit<<" "<<BSA_Error_maxi_fit<<endl;	      

      std::cout<<"type value fit"<<endl;	      
      cout<<"Raw "<<BSA_Amplitude<<" "<<BSA_Amplitude_fit<<" "<<BSA_Error_fit<<endl;	      
      cout<<"Mostafa "<<BSA_Amplitude_most<<" "<<BSA_Amplitude_most_fit<<" "<<BSA_Error_most_fit<<endl;	      
      cout<<"Maxime "<<BSA_Amplitude_maxi<<" "<<BSA_Amplitude_maxi_fit<<" "<<BSA_Error_maxi_fit<<endl;
      
      outFile<<"Entries before BDT (All/FT/FD): "<<entries_bef_BDT<<" "<<entries_bef_BDT_FT<<" "<<entries_bef_BDT_FD<<endl;
      outFile<<"Mostafa entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<entries_bef_most<<" "<<entries_aft_most<<" "<<entries_bef_most_FT<<" "<<entries_bef_most_FD<<" "<<boundaries.at(8)*100<<"% "<<boundaries.at(9)*100<<"%"<<" "<<boundaries.at(10)*100<<"% "<<boundaries.at(11)*100<<"%"<<" "<<boundaries.at(12)*100<<"% "<<boundaries.at(13)*100<<"%"<<endl;  
      outFile<<"Maxime  entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<entries_bef_maxi<<" "<<entries_aft_maxi<<" "<<entries_bef_most_FT<<" "<<entries_bef_most_FD<<" "<<boundaries.at(14)*100<<"% "<<boundaries.at(15)*100<<"% "<<boundaries.at(16)*100<<"% "<<boundaries.at(17)*100<<"%"<<" "<<boundaries.at(18)*100<<"% "<<boundaries.at(19)*100<<"%"<<endl;  
      
      std::cout<<"Entries before BDT (All/FT/FD): "<<entries_bef_BDT<<" "<<entries_bef_BDT_FT<<" "<<entries_bef_BDT_FD<<endl;
      std::cout<<"Mostafa entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<entries_bef_most<<" "<<entries_aft_most<<" "<<entries_bef_most_FT<<" "<<entries_bef_most_FD<<" "<<boundaries.at(8)*100<<"% "<<boundaries.at(9)*100<<"%"<<" "<<boundaries.at(10)*100<<"% "<<boundaries.at(11)*100<<"%"<<" "<<boundaries.at(12)*100<<"% "<<boundaries.at(13)*100<<"%"<<endl;  
      std::cout<<"Maxime  entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<entries_bef_maxi<<" "<<entries_aft_maxi<<" "<<entries_bef_most_FT<<" "<<entries_bef_most_FD<<" "<<boundaries.at(14)*100<<"% "<<boundaries.at(15)*100<<"% "<<boundaries.at(16)*100<<"% "<<boundaries.at(17)*100<<"%"<<" "<<boundaries.at(18)*100<<"% "<<boundaries.at(19)*100<<"%"<<endl;  

      bin_number+=1;

      outFile.close();
      
      if(xsection)
        xsec_on_bins(bin, Orig, Most, Maxi, NBinsPhi);

      delete Orig;
      delete Most;
      delete Maxi;

    }


  Folder = Folder_old;

  return;
}
