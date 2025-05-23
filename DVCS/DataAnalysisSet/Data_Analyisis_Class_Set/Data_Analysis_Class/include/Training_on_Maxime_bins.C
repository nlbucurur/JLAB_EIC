
void BDT::Training_on_Maxime_bins(TString Data)//, TCut cut, TString DVCS, TString Pi0, vector<TString> vars)
{
  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Data);

  TString Folder_old=Folder;
  TCut bins[64];

  bins[0]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[1]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[2]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[3]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[4]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[5]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[6]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[7]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[8]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[9]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[10]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[11]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");

  bins[12]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[13]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[14]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[15]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[16]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[17]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[18]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[19]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[20]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  bins[21]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  bins[22]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  bins[23]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");

  bins[24]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  bins[25]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  bins[26]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  bins[27]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  bins[28]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  bins[29]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  bins[30]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  bins[31]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  bins[32]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  bins[33]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  bins[34]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  bins[35]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");

  bins[36]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  bins[37]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  bins[38]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  bins[39]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  bins[40]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  bins[41]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  bins[42]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  bins[43]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  bins[44]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  bins[45]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  bins[46]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  bins[47]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");

  bins[48]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33");
  bins[49]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  bins[50]= TCut("bestCandidateFlag==1 && t_Ph<0    && t_Ph>-0.2 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33 ");
  bins[51]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0  ");
  bins[52]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.   && strip_Xbj < 0.33");
  bins[53]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  bins[54]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33");
  bins[55]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");

  bins[56]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  bins[57]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  bins[58]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  bins[59]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  bins[60]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  bins[61]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  bins[62]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  bins[63]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");

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
  double xold;
for(int bin_number=1; bin_number<=64;bin_number++)
    {
      TH1 *Orig;
      TH1 *Most;
      //TH1 *Maxi;

      TCut cut_bin = bins[bin_number - 1];
      boundaries.clear();
      boundaries.push_back(1);
      boundaries.push_back(1);
      boundaries.push_back(1);
      boundaries.push_back(1);
      boundaries.push_back(1);
      boundaries.push_back(1);
      Folder = Folder_old + TString("bin_")+Form("%i/",bin_number);
      gSystem->Exec(TString("mkdir -p ") + Folder);

      //Training(cut + cut_bin, DVCS, Pi0, Vars);
      //Training_vars(Data, DVCS, Pi0, cut + cut_bin);
      //Add_BDT_var(cut + cut_bin, Data, TData, Vars);
      //Add_BDT_var(cut + cut_bin, DVCS, TDVCS, Vars);
      //Add_BDT_var(cut + cut_bin, Pi0, TPi0, Vars);
      //Explore(TData, TDVCS, cut + cut_bin);
      //Filter(TData, cut + cut_bin + TCut("_strip_Nuc_BDT > 0.0"), TString("Data_NP_Theta_g_5.root"));

      std::cout<<"\n Get Contamination BDT way"<<endl;
      Get_Contamination(cut + cut_bin, 0.0);

      //Plot Maxime BSA
      inputFile.open(Form("/home/munoz/Documents/RG-A:2/P/Maxime_BSA/bin_%i.txt",bin_number)); // Replace with the name of the block file you want to read
      xold=0;
      phi_bins.clear();
      phi_bins.push_back(xold);
      xValues.clear();
      yValues.clear();
      xErrors.clear();
      yErrors.clear();
      while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        iss >> x >> y >> yErr;
	x=x*180/TMath::Pi();
        // Compute the error on X based on the bin width
        binWidth = (x - xold)*2.0;
        xErr = binWidth / 2.0;
        
        xValues.push_back(x);
        yValues.push_back(y);
        yErrors.push_back(yErr);
        //xErrors.push_back(xErr);
	xErrors.push_back(1.0);

	xold+=binWidth;
	phi_bins.push_back(xold);
	//std::cout<<x<<" "<<xold<<" "<<binWidth<<endl;
      }
      inputFile.close();
      numPoints = xValues.size();
      std::cout<<"Number of BSA points "<<phi_bins.size() -1<<endl;

      std::cout<<"\n Get Contamination Mostafa way"<<endl;
      //Add_BDT_var(cut + cut_bin, "/home/munoz/Datasets/DVMP/Quality_Pi_as_DVCS_P.root", "Tested_1gamma.root", Vars);
      //Add_BDT_var(cut + cut_bin, "/home/munoz/Datasets/DVMP/Quality_Pi_as_DVCS_P_1.root", "Tested_1gamma_1.root", Vars);
	      
      Most=Get_Contamination_Mostafa_Maxime(cut + cut_bin, phi_bins);
	      
      //std::cout<<"\n Get Contamination Maxime way"<<endl;
      //Maxi=Maxime(cut + cut_bin);
      boundaries.push_back(1);
      boundaries.push_back(1);

      Orig=Single_BSA_Maxime("Data_NP_Theta_g_5.root", phi_bins);
      Orig->SetTitle("Before");
      gStyle->SetOptTitle(0);
      TCanvas* c3 = new TCanvas("c3","Histograms");
      //THStack* hs = new THStack("hs","");
      //hs->Add(Orig);
      //hs->Add(Most);
      Most->SetAxisRange(-1., 1.,"Y");
      Orig->SetAxisRange(-1., 1.,"Y");
      
      Most->Draw("E0");
      Orig->Draw("SAME, E0");
      
      //Maxi->Draw("SAME");
      
      TGraphErrors* graph = new TGraphErrors(numPoints, xValues.data(), yValues.data(), xErrors.data(), yErrors.data());
      graph->SetLineColor(kRed);
      graph->SetMinimum(-1.);
      graph->SetMaximum(1.);
      graph->Draw("SAME P");
      graph->SetTitle("RG-A Note");

      //hs->Draw("nostack");
      c3->BuildLegend();

      c3->Print(Folder + TString("Most_vs_Maxi.pdf"));
      gStyle->SetOptTitle(1);

      delete c3;
      delete Orig;
      delete Most;
      //delete Maxi;
      delete graph;
      //delete hs;
    }

  Folder = Folder_old;

  delete ch1;
  return;
}
