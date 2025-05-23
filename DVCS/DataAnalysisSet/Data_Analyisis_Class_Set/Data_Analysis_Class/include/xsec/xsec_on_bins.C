
void BDT::xsec_on_bins(int bin, TH1* Orig, TH1* Most, TH1* Maxi, int NBinsPhi)
{
  //Initial declarations
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);


  std::ifstream inputFile;
  std::string line;
  std::vector<double> xValues;
  std::vector<double> yValues;
  std::vector<double> yErrors;
  std::vector<double> xErrors;
  std::vector<double> phi_bins;

  double x, y, aux, yErr;
  double binWidth;
  double xErr;
  int numPoints;
  
  TString Folder_old=Folder;
  
  cut_bin = Mbins[bin-1];
  gSystem->Exec(TString("mkdir -p ") + Folder + TString("xsec/"));
  Folder = Folder + TString("xsec/");

  std::vector<vector<double>> means1, means2;
  means1 = ReadMeansFile(Folder + TString("../means_most.txt"));
  means2 = ReadMeansFile(Folder + TString("../means_maxi.txt"));

  //Get shift of kinematics
  TGraphErrors* t_bef1 = new TGraphErrors();
  TGraphErrors* Q_bef1 = new TGraphErrors();
  TGraphErrors* x_bef1 = new TGraphErrors();
  TGraphErrors* t_aft1 = new TGraphErrors();
  TGraphErrors* Q_aft1 = new TGraphErrors();
  TGraphErrors* x_aft1 = new TGraphErrors();

  TGraphErrors* t_bef2 = new TGraphErrors();
  TGraphErrors* Q_bef2 = new TGraphErrors();
  TGraphErrors* x_bef2 = new TGraphErrors();
  TGraphErrors* t_aft2 = new TGraphErrors();
  TGraphErrors* Q_aft2 = new TGraphErrors();
  TGraphErrors* x_aft2 = new TGraphErrors();

  TH1F* factor_t = compute_kin_shift_factor(NBinsPhi, "t_Ph");
  TH1F* factor_Q = compute_kin_shift_factor(NBinsPhi, "strip_Q2");
  TH1F* factor_x = compute_kin_shift_factor(NBinsPhi, "strip_Xbj");
  TH1F* factor_p = compute_kin_shift_factor(NBinsPhi, "Phi_Ph");

  t_bef1->SetTitle("Before 1;#phi (deg); t (GeV^{2})");
  Q_bef1->SetTitle("Before 1;#phi (deg); Q (GeV^{2})");
  x_bef1->SetTitle("Before 1;#phi (deg); x_{B}");
  t_aft1->SetTitle("Before 1;#phi (deg); t (GeV^{2})");
  Q_aft1->SetTitle("Before 1;#phi (deg); Q (GeV^{2})");
  x_aft1->SetTitle("Before 1;#phi (deg); x_{B}");
  t_bef2->SetTitle("Before 2;#phi (deg); t (GeV^{2})");
  Q_bef2->SetTitle("Before 2;#phi (deg); Q (GeV^{2})");
  x_bef2->SetTitle("Before 2;#phi (deg); x_{B}");
  t_aft2->SetTitle("Before 2;#phi (deg); t (GeV^{2})");
  Q_aft2->SetTitle("Before 2;#phi (deg); Q (GeV^{2})");
  x_aft2->SetTitle("Before 2;#phi (deg); x_{B}");
  factor_t->SetTitle("Factor t;#phi (deg); Factor");
  factor_Q->SetTitle("Factor Q;#phi (deg); Factor");
  factor_x->SetTitle("Factor x;#phi (deg); Factor");
  factor_p->SetTitle("Factor p;#phi (deg); Factor");

  t_bef1->SetName("t_mean1_bef");
  Q_bef1->SetName("Q_mean1_bef");
  x_bef1->SetName("x_mean1_bef");
  t_aft1->SetName("t_mean1_aft");
  Q_aft1->SetName("Q_mean1_aft");
  x_aft1->SetName("x_mean1_aft");
  t_bef2->SetName("t_mean2_bef");
  Q_bef2->SetName("Q_mean2_bef");
  x_bef2->SetName("x_mean2_bef");
  t_aft2->SetName("t_mean2_aft");
  Q_aft2->SetName("Q_mean2_aft");
  x_aft2->SetName("x_mean2_aft");


  t_bef1->SetLineColor(kBlack);
  Q_bef1->SetLineColor(kBlack);
  x_bef1->SetLineColor(kBlack);
  t_aft1->SetLineColor(kRed);
  Q_aft1->SetLineColor(kRed);
  x_aft1->SetLineColor(kRed);
  t_bef2->SetLineColor(kBlack);
  Q_bef2->SetLineColor(kBlack);
  x_bef2->SetLineColor(kBlack);
  t_aft2->SetLineColor(kRed);
  Q_aft2->SetLineColor(kRed);
  x_aft2->SetLineColor(kRed);

  t_bef1->SetMarkerColor(kBlack);
  Q_bef1->SetMarkerColor(kBlack);
  x_bef1->SetMarkerColor(kBlack);
  t_aft1->SetMarkerColor(kRed);
  Q_aft1->SetMarkerColor(kRed);
  x_aft1->SetMarkerColor(kRed);
  t_bef2->SetMarkerColor(kBlack);
  Q_bef2->SetMarkerColor(kBlack);
  x_bef2->SetMarkerColor(kBlack);
  t_aft2->SetMarkerColor(kRed);
  Q_aft2->SetMarkerColor(kRed);
  x_aft2->SetMarkerColor(kRed);
  
  for(int k=0; k<NBinsPhi; k++)
  {
    t_bef1->SetPoint(k, means1.at(k).at(0),means1.at(k).at(2));
    t_bef1->SetPointError(k, means1.at(k).at(1),means1.at(k).at(3));
    Q_bef1->SetPoint(k, means1.at(k).at(0),means1.at(k).at(4));
    Q_bef1->SetPointError(k, means1.at(k).at(1),means1.at(k).at(5));
    x_bef1->SetPoint(k, means1.at(k).at(0),means1.at(k).at(6));
    x_bef1->SetPointError(k, means1.at(k).at(1),means1.at(k).at(7));

    t_bef2->SetPoint(k, means2.at(k).at(0),means2.at(k).at(2));
    t_bef2->SetPointError(k, means2.at(k).at(1),means2.at(k).at(3));
    Q_bef2->SetPoint(k, means2.at(k).at(0),means2.at(k).at(4));
    Q_bef2->SetPointError(k, means2.at(k).at(1),means2.at(k).at(5));
    x_bef2->SetPoint(k, means2.at(k).at(0),means2.at(k).at(6));
    x_bef2->SetPointError(k, means2.at(k).at(1),means2.at(k).at(7));

    t_aft1->SetPoint(k, means1.at(k).at(0)*factor_p->GetBinContent(k+1),means1.at(k).at(2)*factor_t->GetBinContent(k+1));
    t_aft1->SetPointError(k, means1.at(k).at(1)*factor_p->GetBinContent(k+1),means1.at(k).at(3)*factor_t->GetBinContent(k+1));
    Q_aft1->SetPoint(k, means1.at(k).at(0)*factor_p->GetBinContent(k+1),means1.at(k).at(4)*factor_Q->GetBinContent(k+1));
    Q_aft1->SetPointError(k, means1.at(k).at(1)*factor_p->GetBinContent(k+1),means1.at(k).at(5)*factor_Q->GetBinContent(k+1));
    x_aft1->SetPoint(k, means1.at(k).at(0)*factor_p->GetBinContent(k+1),means1.at(k).at(6)*factor_x->GetBinContent(k+1));
    x_aft1->SetPointError(k, means1.at(k).at(1)*factor_p->GetBinContent(k+1),means1.at(k).at(7)*factor_x->GetBinContent(k+1));

    t_aft2->SetPoint(k, means2.at(k).at(0)*factor_p->GetBinContent(k+1),means2.at(k).at(2)*factor_t->GetBinContent(k+1));
    t_aft2->SetPointError(k, means2.at(k).at(1)*factor_p->GetBinContent(k+1),means2.at(k).at(3)*factor_t->GetBinContent(k+1));
    Q_aft2->SetPoint(k, means2.at(k).at(0)*factor_p->GetBinContent(k+1),means2.at(k).at(4)*factor_Q->GetBinContent(k+1));
    Q_aft2->SetPointError(k, means2.at(k).at(1)*factor_p->GetBinContent(k+1),means2.at(k).at(5)*factor_Q->GetBinContent(k+1));
    x_aft2->SetPoint(k, means2.at(k).at(0)*factor_p->GetBinContent(k+1),means2.at(k).at(6)*factor_x->GetBinContent(k+1));
    x_aft2->SetPointError(k, means2.at(k).at(1)*factor_p->GetBinContent(k+1),means2.at(k).at(7)*factor_x->GetBinContent(k+1));

  }

  TH1F* factor3D_t = compute_kin_shift_factor("t_Ph");
  TH1F* factor3D_Q = compute_kin_shift_factor("strip_Q2");
  TH1F* factor3D_x = compute_kin_shift_factor("strip_Xbj");

  TCanvas* canvas = new TCanvas("canvas","canvas",700,500);
  t_bef1->Draw("ap");
  t_aft1->Draw("p");
  TLegend* legend1 = new TLegend(0.6,0.6,0.88,0.88);
  legend1->AddEntry(t_bef1,"Before","l");
  legend1->AddEntry(t_aft1,"After","l");
  legend1->Draw();
  canvas->Print(Folder + TString("t_corr1.pdf"));

  Q_bef1->Draw("ap");
  Q_aft1->Draw("p");
  canvas->Print(Folder + TString("Q_corr1.pdf"));

  x_bef1->Draw("ap");
  x_aft1->Draw("p");
  canvas->Print(Folder + TString("x_corr1.pdf"));

  t_bef2->Draw("ap");
  t_aft2->Draw("p");
  canvas->Print(Folder + TString("t_corr2.pdf"));

  Q_bef2->Draw("ap");
  Q_aft2->Draw("p");
  canvas->Print(Folder + TString("Q_corr2.pdf"));

  x_bef2->Draw("ap");
  x_aft2->Draw("p");
  canvas->Print(Folder + TString("x_corr2.pdf"));

  //root output
  TFile *outputFileMeans = new TFile(Folder + TString("means.root"), "RECREATE");
  t_bef1->Write();
  Q_bef1->Write();
  x_bef1->Write();
  t_aft1->Write();
  Q_aft1->Write();
  x_aft1->Write();
  t_bef2->Write();
  Q_bef2->Write();
  x_bef2->Write();
  t_aft2->Write();
  Q_aft2->Write();
  x_aft2->Write();
  factor_t->Write();
  factor_Q->Write();
  factor_x->Write();
  factor_p->Write();
  factor3D_t->Write();
  factor3D_Q->Write();
  factor3D_x->Write();
  outputFileMeans->Close();


  delete t_bef1;
  delete Q_bef1;
  delete x_bef1;
  delete t_aft1;
  delete Q_aft1;
  delete x_aft1;
  delete t_bef2;
  delete Q_bef2;
  delete x_bef2;
  delete t_aft2;
  delete Q_aft2;
  delete x_aft2;

  //Transform BSA to event count
  TH1* NOrig = BSA2Nev(Orig);
  TH1* NMost = BSA2Nev(Most);
  TH1* NMaxi = BSA2Nev(Maxi);
	  
  NOrig->SetMarkerColor(kBlack);
  NOrig->SetLineColor(kBlack);

  NMost->SetMarkerColor(kBlue);
  NMost->SetLineColor(kBlue);

  NMaxi->SetMarkerColor(kRed);
  NMaxi->SetLineColor(kRed);

  NMost->SetTitle("Method 1");
  NMaxi->SetTitle("Method 2");
  
  NOrig->Draw();
  NMost->Draw("SAME");
  NMaxi->Draw("SAME");

  NOrig->SetName("Before");
  NMost->SetName("Method 1");
  NMaxi->SetName("Method 2");

  //Compute errors of Nevent histograms
  TH1F* OrigErr1 = new TH1F("OrigErr1","OrigErr1",NBinsPhi,0,360);
  TH1F* OrigErr2 = new TH1F("OrigErr2","OrigErr2",NBinsPhi,0,360);
  std::ifstream StatFile;
  StatFile.open(Folder + TString("../Bkg_Most_stats.txt")); // Replace with the name of the block file you want to read
  double n1g_p, n1g_m, n2g_p, n2g_m, n2gD_p, n2gD_m, nData_p, nData_m;
  double n1g, n2g, n2gD;

  for(int t=1; t<=NBinsPhi; t++)
    {
      std::getline(StatFile, line);
      std::istringstream iss(line);
      iss >> n1g_p >> n1g_m >> n2g_p >> n2g_m >> n2gD_p >> n2gD_m >> nData_p >> nData_m;
      n1g = (n1g_p + n1g_m);
      n2g = (n2g_p + n2g_m);
      n2gD = (n2gD_p + n2gD_m);

      if(NMost->GetBinContent(t)>1 && NOrig->GetBinContent(t)>1 && n1g>1 && n2g>1 && n2gD>1)
      {
        OrigErr1->SetBinContent(t, 1.0);
        OrigErr1->SetBinError(t, NOrig->GetBinContent(t)*sqrt(1.0/(NOrig->GetBinContent(t)*pow(NMost->GetBinContent(t)/NOrig->GetBinContent(t),2)) ));  
        NMost->SetBinError(t, NMost->GetBinContent(t)*sqrt(pow(NOrig->GetBinContent(t)/NMost->GetBinContent(t) - 1,2)*(1./n1g + 1./n2g + 1./n2gD) ));
        }
      else
      {
        OrigErr1->SetBinContent(t, 0.0);
        OrigErr1->SetBinError(t, 0.0);
        NMost->SetBinError(t, 0.0);
      }
      if(NMaxi->GetBinContent(t)>1 && NOrig->GetBinContent(t)>1 && n2gD>1)
      {
        OrigErr2->SetBinContent(t, 1.0);
        OrigErr2->SetBinError(t, NOrig->GetBinContent(t)*sqrt(1.0/(NOrig->GetBinContent(t)*pow(NMaxi->GetBinContent(t)/NOrig->GetBinContent(t),2)) ));
        NMaxi->SetBinError(t, NMaxi->GetBinContent(t)*sqrt(pow(NOrig->GetBinContent(t)/NMaxi->GetBinContent(t) - 1,2)*(1./NMaxi->GetBinContent(t) + 1./n2gD) ));  
      }
      else
      {
        OrigErr2->SetBinContent(t, 0.0);
        OrigErr2->SetBinError(t, 0.0);
        NMaxi->SetBinError(t, 0.0);
      }
      std::cout<<NOrig->GetBinContent(t)<<" "<<OrigErr1->GetBinError(t)<<" "<<NMost->GetBinContent(t)<<" "<<NMost->GetBinError(t)<<" "<<NMaxi->GetBinContent(t)<<" "<<NMaxi->GetBinError(t)<<endl;
    }
  StatFile.close();

  canvas->BuildLegend();
  canvas->Print(Folder + TString("Background_subtraction.pdf"));
  delete canvas;
	  
  //Theory xsection 
  
  Theory(1,bin, xmean1*factor3D_x->GetBinContent(1), Qmean1*factor3D_Q->GetBinContent(1), tmean1*factor3D_t->GetBinContent(1), NBinsPhi); //BH
  Theory(2,bin, xmean1*factor3D_x->GetBinContent(1), Qmean1*factor3D_Q->GetBinContent(1), tmean1*factor3D_t->GetBinContent(1), NBinsPhi); //VGG
  Theory(3,bin, xmean1*factor3D_x->GetBinContent(1), Qmean1*factor3D_Q->GetBinContent(1), tmean1*factor3D_t->GetBinContent(1), NBinsPhi); //KM15
  //Theory(1,bin, 0.5*(bins[bin-1][4]+bins[bin-1][5]), 0.5*(bins[bin-1][2]+bins[bin-1][3]), 0.5*(bins[bin-1][0]+bins[bin-1][1]), NBinsPhi); //BH
  //Theory(2,bin, 0.5*(bins[bin-1][4]+bins[bin-1][5]), 0.5*(bins[bin-1][2]+bins[bin-1][3]), 0.5*(bins[bin-1][0]+bins[bin-1][1]), NBinsPhi); //VGG
  //Theory(3,bin, 0.5*(bins[bin-1][4]+bins[bin-1][5]), 0.5*(bins[bin-1][2]+bins[bin-1][3]), 0.5*(bins[bin-1][0]+bins[bin-1][1]), NBinsPhi); //KM15
  //Theory(4,bin); //GK19
  
  //Plot KM xsec
  inputFile.open(extXSEC + Form("KM/bin_%i.txt",bin)); // Replace with the name of the block file you want to read
  xValues.clear();
  yValues.clear();
  while (std::getline(inputFile, line)) {
	  std::istringstream iss(line);
	  iss >> x >> aux >> aux >> aux >> aux >> y;
    
	  xValues.push_back(x);
	  yValues.push_back(y);
  }
  numPoints = xValues.size();
  inputFile.close();
  TGraph* graph = new TGraph(numPoints, xValues.data(), yValues.data());
  graph->SetTitle("KM15");
  graph->SetMarkerColor(kCyan);
  graph->SetLineColor(kCyan);
  graph->SetLineWidth(2);
  graph->SetTitle("KM15");

  //Plot VGG xsec
  inputFile.open(extXSEC + Form("VGG/bin_%i.txt",bin)); // Replace with the name of the block file you want to read
  xValues.clear();
  yValues.clear();
  while (std::getline(inputFile, line)) {
	  std::istringstream iss(line);
	  iss >> x >> aux >> aux >> aux >> aux >> y;
    
	  xValues.push_back(x);
	  yValues.push_back(y);
  }
  numPoints = xValues.size();
  inputFile.close();
  TGraph* graph2 = new TGraph(numPoints, xValues.data(), yValues.data());
  graph2->SetTitle("VGG");
  graph2->SetMarkerColor(kBlue);
  graph2->SetLineColor(kBlue);
  graph2->SetLineWidth(2);
  graph2->SetTitle("VGG");

  //Plot GK xsec
  inputFile.open(extXSEC + Form("GK/bin_%i.txt",bin)); // Replace with the name of the block file you want to read
  xValues.clear();
  yValues.clear();
  while (std::getline(inputFile, line)) {
	  std::istringstream iss(line);
	  iss >> x >> aux >> aux >> aux >> aux >> y;
    
	  xValues.push_back(x);
	  yValues.push_back(y);
  }
  numPoints = xValues.size();
  inputFile.close();
  TGraph* graph3 = new TGraph(numPoints, xValues.data(), yValues.data());
  graph3->SetTitle("GK19");
  graph3->SetMarkerColor(kRed);
  graph3->SetLineColor(kRed);
  graph3->SetLineWidth(2);
  graph3->SetTitle("GK19");

  //Plot BH xsec
  inputFile.open(extXSEC + Form("BH/bin_%i.txt",bin)); // Replace with the name of the block file you want to read
  xValues.clear();
  yValues.clear();
  while (std::getline(inputFile, line)) {
	  std::istringstream iss(line);
	  iss >> x >> aux >> aux >> aux >> aux >> y;
    
	  xValues.push_back(x);
	  yValues.push_back(y);
  }
  numPoints = xValues.size();
  inputFile.close();
  
  TGraph* graph4 = new TGraph(numPoints, xValues.data(), yValues.data());
  graph4->SetTitle("BH");
  graph4->SetMarkerColor(kRed);
  graph4->SetLineColor(kRed);
  graph4->SetLineWidth(2);
  graph4->SetLineStyle(kDashed);
  graph4->SetTitle("BH");

  gStyle->SetOptTitle(1);
	  
  //Corrections	  
  TH1F* F_noRec = No_Rec_correction(Folder + TString("../../") + TBM_Sim, MC_BM_Sim, BDT_value, NBinsPhi);
  F_noRec->SetName("F_noRec");

  TH1F* F_BM1 = BM_correction(NMost, bin,1);
  F_BM1->SetName("F_BM_1");
  TH1F* F_BM2 = BM_correction(NMaxi, bin,2);
  F_BM2->SetName("F_BM_2");

  
  TH1F* F_bin_KM_model1 =  F_Bin( bin, factor_t, factor_Q, factor_x, 1);
  F_bin_KM_model1->SetName("F_Bin_KM_Model1");
  TH1F* F_bin_KM_model2 =  F_Bin( bin, factor_t, factor_Q, factor_x, 2);
  F_bin_KM_model2->SetName("F_Bin_KM_Model2");
  
  /*
  TH1F* F_bin_KM_model1 = new TH1F("F_bin_KM_model1","F_bin_KM_model1",NBinsPhi,0,360);
  TH1F* F_bin_KM_model2 = new TH1F("F_bin_KM_model2","F_bin_KM_model2",NBinsPhi,0,360);
  for(int t=1; t<=NBinsPhi; t++)
    {
      F_bin_KM_model1->SetBinContent(t, 1.0);
      F_bin_KM_model1->SetBinError(t, 0.1);
      F_bin_KM_model2->SetBinContent(t, 1.0);
      F_bin_KM_model2->SetBinError(t, 0.1);
    }
  */

  TH1F* F_rad = F_RC(NBinsPhi);
  F_rad->SetName("F_RC");
  
  /*
  TH1F* F_423_1 = F_4Dto3D(xmean1, Qmean1, tmean1,NBinsPhi,1);
  F_423_1->SetName("F_423_1");
  TH1F* F_423_2 = F_4Dto3D(xmean2, Qmean2, tmean2,NBinsPhi,2);
  F_423_2->SetName("F_423_2");

  */
  double L=40.09*1e6 ; //40.09 fb-1 in nb-1
  double PS_factor = (boundaries.at(1) - boundaries.at(0))*(boundaries.at(3) - boundaries.at(2))*(boundaries.at(5) - boundaries.at(4))*(2*TMath::Pi()/NBinsPhi);
  std::cout<<"Luminosity times PS_factor "<<L*PS_factor<<endl;

  TH1F *xsec1 = new TH1F("xsec1","xsec1",NBinsPhi,0,360);
  TH1F *xsec2 = new TH1F("xsec2","xsec2",NBinsPhi,0,360);
  
  xsec1->Divide(NMost, F_BM1,1,1);
  xsec2->Divide(NMaxi, F_BM2,1,1);

  xsec1->Divide(xsec1, F_noRec,1,1);
  xsec2->Divide(xsec2, F_noRec,1,1);

  xsec1->Divide(xsec1, F_bin_KM_model1,1,1);
  xsec2->Divide(xsec2, F_bin_KM_model2,1,1);

  xsec1->Divide(xsec1, F_rad,1,1);
  xsec2->Divide(xsec2, F_rad,1,1);

  /*
  xsec1->Divide(xsec1, F_423_1,1,1);
  xsec2->Divide(xsec2, F_423_2,1,1);
  */
    
  xsec1->Scale(1.0/(L*PS_factor));
  xsec2->Scale(1.0/(L*PS_factor));

  double F_eff1;
  double F_eff2;
  int counter=0;
  std::cout<<"Normalization factor"<<endl;
  for(int t=0; t<NBinsPhi; t++)
  {
    if(graph4->GetPointX(t)<30 || graph4->GetPointX(t)>330)
    {
      F_eff1+=xsec1->GetBinContent(t+1)/graph4->GetPointY(t);
      F_eff2+=xsec2->GetBinContent(t+1)/graph4->GetPointY(t);
      std::cout<<xsec1->GetBinContent(t+1)<<" "<<graph4->GetPointY(t)<<endl;
      counter++;
    }
  }
  F_eff1=F_eff1/counter;
  F_eff2=F_eff2/counter;
  std::cout<<"Global normalization factor "<<F_eff1<<" "<<F_eff2<<endl;

  xsec1->Scale(1.0/F_eff1);
  xsec2->Scale(1.0/F_eff2);

  xsec1->SetTitle("Cross-section (nb); #phi (deg); #sigma (nb)");
  xsec2->SetTitle("Cross-section (nb); #phi (deg); #sigma (nb)");
  xsec1->SetLineColor(kBlack);
  xsec2->SetLineColor(kBlue);
  xsec1->SetMarkerColor(kBlack);
  xsec2->SetMarkerColor(kBlue);

  //Final xsection
  TGraphErrors* xsec1gr = new TGraphErrors();
  TGraphErrors* xsec2gr = new TGraphErrors();
  double err1=0, err2=0;

  for(int t=1; t<=NBinsPhi; t++)
    {
      err1=0;
      err2=0;
      if(NOrig->GetBinContent(t)!=0)
        err1+= pow(OrigErr1->GetBinError(t)/NOrig->GetBinContent(t),2);
      if(NMost->GetBinContent(t)!=0)
        err1+= pow(NMost->GetBinError(t)/NMost->GetBinContent(t),2) ;
      if(F_noRec->GetBinContent(t)!=0)
        err1+= pow(F_noRec->GetBinError(t)/F_noRec->GetBinContent(t),2) ;
      if(F_rad->GetBinContent(t)!=0)
        err1+= pow(F_rad->GetBinError(t)/F_rad->GetBinContent(t),2) ;

      if(NOrig->GetBinContent(t)!=0)
        err2+= pow(OrigErr2->GetBinError(t)/NOrig->GetBinContent(t),2) ;
      if(NMaxi->GetBinContent(t)!=0)
        err2+= pow(NMaxi->GetBinError(t)/NMaxi->GetBinContent(t),2) ;
      if(F_noRec->GetBinContent(t)!=0)
        err2+= pow(F_noRec->GetBinError(t)/F_noRec->GetBinContent(t),2) ;
      if(F_rad->GetBinContent(t)!=0)
        err2+= pow(F_rad->GetBinError(t)/F_rad->GetBinContent(t),2) ;

      err1 = sqrt(err1) * xsec1->GetBinContent(t);
      err2 = sqrt(err2) * xsec2->GetBinContent(t);

      xsec1->SetBinError(t,err1);
      xsec2->SetBinError(t,err2);

      std::cout<<" phi "<<means1.at(t-1).at(0)<<" "<<means2.at(t-1).at(0)<<" xsection "<<xsec1->GetBinContent(t)<<" "<<xsec2->GetBinContent(t)<<" Stat errors "<<err1<<" "<<err2<<endl;
      if(means1.at(t-1).at(0)>1)
      {
        xsec1gr->SetPoint(xsec1gr->GetN(), means1.at(t-1).at(0), xsec1->GetBinContent(t));
        xsec1gr->SetPointError(xsec1gr->GetN()-1, means1.at(t-1).at(1), err1);  
      }
      if(means2.at(t-1).at(0)>1)
      {
      xsec2gr->SetPoint(xsec2gr->GetN(), means2.at(t-1).at(0), xsec2->GetBinContent(t));
      xsec2gr->SetPointError(xsec2gr->GetN()-1, means2.at(t-1).at(1), err2);
      }
    }
  xsec1gr->SetTitle("Cross-section (nb); #phi (deg); #sigma (nb)");
  xsec2gr->SetTitle("Cross-section (nb); #phi (deg); #sigma (nb)");
  xsec1gr->SetName("xsec1gr");    
  xsec2gr->SetName("xsec2gr");    

  // Write the histogram to the output file
  //root output
  TFile *outputFile = new TFile(Folder + TString("xsec_dists.root"), "RECREATE");
  NMost->Write();
  NMaxi->Write();
  NOrig->Write();

  F_BM1->Write();
  F_BM2->Write();
  F_noRec->Write();

  F_bin_KM_model1->Write();
  F_bin_KM_model2->Write();

  F_rad->Write();

  /*
  F_423_1->Write();
  F_423_2->Write();
  */

  xsec1->Write();
  xsec2->Write();
  xsec1gr->Write();
  xsec2gr->Write();
  graph->Write();
  graph2->Write();
  graph3->Write();
  graph4->Write();
  // Close the output file
  outputFile->Close();
  std::cout<<"plots saved succesfully !"<<endl;

  //Summary	  
  xsec1gr->SetMarkerStyle(20);
  xsec2gr->SetMarkerStyle(20);
  xsec1gr->SetLineColor(kBlack);
  xsec2gr->SetLineColor(kBlue);
  xsec1gr->SetMarkerColor(kBlack);
  xsec2gr->SetMarkerColor(kBlue);
  xsec1gr->SetMarkerSize(0.5);
  xsec2gr->SetMarkerSize(0.5);
  xsec1gr->SetMarkerSize(0.5);
  xsec2gr->SetMarkerSize(0.5);
  xsec1gr->SetLineWidth(2);
  xsec2gr->SetLineWidth(2);

  xsec1gr->GetXaxis()->SetLimits(0,360);
  xsec2gr->GetXaxis()->SetLimits(0,360);

  TCanvas *canvas2 = new TCanvas("canvas2","canvas2",700,500);
  TMultiGraph* mg = new TMultiGraph();
  //gPad->SetLogy();
  mg->Add(xsec1gr,"P");
  mg->Add(xsec2gr,"P");
  mg->Add(graph,"l");
  mg->Add(graph2,"l");
  //mg->Add(graph3,"l");
  mg->Add(graph4,"l");
  mg->Draw("A");
 
  mg->SetTitle("Cross-section (nb); #phi (deg); #sigma (nb)");
  mg->GetXaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleOffset(0.7);
  mg->GetXaxis()->SetTitleOffset(0.5);
  mg->GetYaxis()->SetNdivisions(6);
  mg->GetXaxis()->SetNdivisions(4);

  TLegend* legend = new TLegend(0.3,0.6,0.6,0.88);
  legend->AddEntry(xsec1,"Method 1","lep");
  legend->AddEntry(xsec2,"Method 2","lep");
  legend->AddEntry(graph,"KM15","l");
  legend->AddEntry(graph2,"VGG","l");
  //legend->AddEntry(graph3,"GK19","l");
  legend->AddEntry(graph4,"BH","l");
  legend->Draw();

  canvas2->Print(Folder + TString("RGA_vs_This_xsec.pdf"));
  
  delete xsec1;
  delete xsec2;

  delete F_BM1;
  delete F_BM2;
  delete F_noRec;
  delete F_bin_KM_model1;
  delete F_bin_KM_model2;
  delete F_rad;
  
  /*
  delete F_423_1;
  delete F_423_2;
  */
  
  delete NOrig;
  delete NMost;
  delete NMaxi;
  delete OrigErr1;
  delete OrigErr2;

  delete factor_t;
  delete factor_Q;
  delete factor_x;
  delete factor_p;

  delete factor3D_t;
  delete factor3D_Q;
  delete factor3D_x;

  delete canvas2;
  
  Folder = Folder_old;

  return;
}
