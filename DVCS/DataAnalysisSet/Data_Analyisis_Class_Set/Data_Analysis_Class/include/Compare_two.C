void BDT::Compare_two(TCut cut1, TCut cut2, TString Data1, TString Data2, TString output, TString path1="", TString path2="", int set=1, bool normalize=false){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","Histograms");
  c1->Divide(3,3);

  TFile *input2 = new TFile(path1 + Data1,"READ");
  TTree *tree1 = (TTree*)input2->Get("pDVCS");
  
  TFile *input3 = new TFile(path2 + Data2,"READ");
  std::cout<<path2 + Data2<<endl;
  TTree *tree2 = (TTree*)input3->Get("pDVCS");

  double factor=1;
  
  THStack* hs1;
  THStack* hs2;
  THStack* hs3;
  THStack* hs4;
  THStack* hs5;
  THStack* hs6;
  THStack* hs7;
  THStack* hs8;
  THStack* hs9;
  
  TH1F* hist1c_1;
  TH1F* hist2c_1;
  TH1F* hist3c_1;
  TH1F* hist4c_1;
  TH1F* hist5c_1;
  TH1F* hist6c_1;
  TH1F* hist7c_1;
  TH1F* hist8c_1;
  TH1F* hist9c_1;

  TH1F* hist1c_2;
  TH1F* hist2c_2;
  TH1F* hist3c_2;
  TH1F* hist4c_2;
  TH1F* hist5c_2;
  TH1F* hist6c_2;
  TH1F* hist7c_2;
  TH1F* hist8c_2;
  TH1F* hist9c_2;

  const char *p1;
  const char *p2;
  const char *p3;
  const char *p4;
  const char *p5;
  const char *p6;
  const char *p7;
  const char *p8;
  const char *p9;

  if(set==1)
  {
  //Exclusivity variables
  hs1 = new THStack("hs1","#bf{q}'(GeV)");
  hs2 = new THStack("hs2","#bf{k}'(GeV)");
  hs3 = new THStack("hs3","M^{2}_{ep#gamma X} (GeV^{2})");
  hs4 = new THStack("hs4","M^{2}_{e#gamma X} (GeV^{2})");
  hs5 = new THStack("hs5","#Delta #phi (deg)");
  hs6 = new THStack("hs6","#Delta t (GeV^{2})");
  hs7 = new THStack("hs7","P_{miss}(ep#gamma) (GeV)");
  hs8 = new THStack("hs8","#theta_{#gamma X} (deg)");
  hs9 = new THStack("hs9","P_{#perp} (GeV)");

  //Variables to plot
  p1="strip_Ph_P";
  p2="strip_El_P";
  p3="mm2_eNg";
  p4="mm2_eg";
  p5="delta_Phi";
  p6="delta_t";
  p7="miss_mom_eNg";
  p8="theta_gamma_X";
  p9="p_perp";

  hist1c_2 = new TH1F("hist1c_2","",100,2,10);
  hist2c_2 = new TH1F("hist2c_2","",100,0,7);
  hist3c_2 = new TH1F("hist3c_2","",100,-0.1,0.1);
  hist4c_2 = new TH1F("hist4c_2","",100,0,3);
  hist5c_2 = new TH1F("hist5c_2","",100,-0.5,0.5);
  hist6c_2 = new TH1F("hist6c_2","",100,-0.2,0.4);  
  hist7c_2 = new TH1F("hist7c_2","",100,0,0.5);
  hist8c_2 = new TH1F("hist8c_2","",100,0,2);
  hist9c_2 = new TH1F("hist9c_2","",100,0,0.06);

  hist1c_1 = new TH1F("hist1c_1","",100,2,10);
  hist2c_1 = new TH1F("hist2c_1","",100,0,7);
  hist3c_1 = new TH1F("hist3c_1","",100,-0.1,0.1);
  hist4c_1 = new TH1F("hist4c_1","",100,0,3);
  hist5c_1 = new TH1F("hist5c_1","",100,-0.5,0.5);
  hist6c_1 = new TH1F("hist6c_1","",100,-0.2,0.4);  
  hist7c_1 = new TH1F("hist7c_1","",100,0,0.5);
  hist8c_1 = new TH1F("hist8c_1","",100,0,2);
  hist9c_1 = new TH1F("hist9c_1","",100,0,0.06);
  }
  else if(set==2)
  {
  //Kinematics
  hs1 = new THStack("hs1","#bf{q}'(GeV)");
  hs2 = new THStack("hs2","#bf{k}'(GeV)");
  hs3 = new THStack("hs3","#bf{p}'(GeV)");
  hs4 = new THStack("hs4","#theta_{#gamma} (deg)");
  hs5 = new THStack("hs5","#theta_{e} (deg)");
  hs6 = new THStack("hs6","#theta_{p} (deg)");
  hs7 = new THStack("hs7","#phi_{#gamma} (deg)");
  hs8 = new THStack("hs8","#phi_{e} (deg)");
  hs9 = new THStack("hs9","#phi_{p} (deg)");

  //Variables to plot
  p1="strip_Ph_P";
  p2="strip_El_P";
  p3="strip_Nuc_P";
  p4="strip_Ph_Theta";
  p5="strip_El_Theta";
  p6="strip_Nuc_Theta";
  p7="strip_Ph_Phi";
  p8="strip_El_Phi";
  p9="strip_Nuc_Phi";

  hist1c_2 = new TH1F("hist1c_2","",100,2,10);
  hist2c_2 = new TH1F("hist2c_2","",100,0,7);
  hist3c_2 = new TH1F("hist3c_2","",100,0,3.0);
  hist4c_2 = new TH1F("hist4c_2","",100,0,35);
  hist5c_2 = new TH1F("hist5c_2","",100,4,35);
  hist6c_2 = new TH1F("hist6c_2","",100,4,70);  
  hist7c_2 = new TH1F("hist7c_2","",100,-180,180);
  hist8c_2 = new TH1F("hist8c_2","",100,-180,180);
  hist9c_2 = new TH1F("hist9c_2","",100,-180,180);

  hist1c_1 = new TH1F("hist1c_1","",100,2,10);
  hist2c_1 = new TH1F("hist2c_1","",100,0,7);
  hist3c_1 = new TH1F("hist3c_1","",100,0,3.0);
  hist4c_1 = new TH1F("hist4c_1","",100,0,35);
  hist5c_1 = new TH1F("hist5c_1","",100,4,35);
  hist6c_1 = new TH1F("hist6c_1","",100,4,70);  
  hist7c_1 = new TH1F("hist7c_1","",100,-180,180);
  hist8c_1 = new TH1F("hist8c_1","",100,-180,180);
  hist9c_1 = new TH1F("hist9c_1","",100,-180,180);
  }


  tree2->Project("hist1c_2", p1, cut2);
  tree2->Project("hist2c_2", p2, cut2);
  tree2->Project("hist3c_2", p3, cut2);
  tree2->Project("hist4c_2", p4, cut2);
  tree2->Project("hist5c_2", p5, cut2);
  tree2->Project("hist6c_2", p6, cut2);
  tree2->Project("hist7c_2", p7, cut2);
  tree2->Project("hist8c_2", p8, cut2);
  tree2->Project("hist9c_2", p9, cut2);

  tree1->Project("hist1c_1", p1, cut1);
  tree1->Project("hist2c_1", p2, cut1);
  tree1->Project("hist3c_1", p3, cut1);
  tree1->Project("hist4c_1", p4, cut1);
  tree1->Project("hist5c_1", p5, cut1);
  tree1->Project("hist6c_1", p6, cut1);
  tree1->Project("hist7c_1", p7, cut1);
  tree1->Project("hist8c_1", p8, cut1);
  tree1->Project("hist9c_1", p9, cut1);

  hs1->Add(hist1c_1);
  hs1->Add(hist1c_2);
  hs2->Add(hist2c_1);
  hs2->Add(hist2c_2);
  hs3->Add(hist3c_1);
  hs3->Add(hist3c_2);
  hs4->Add(hist4c_1);
  hs4->Add(hist4c_2);
  hs5->Add(hist5c_1);
  hs5->Add(hist5c_2);
  hs6->Add(hist6c_1);
  hs6->Add(hist6c_2);
  hs7->Add(hist7c_1);
  hs7->Add(hist7c_2);
  hs8->Add(hist8c_1);
  hs8->Add(hist8c_2);
  hs9->Add(hist9c_1);
  hs9->Add(hist9c_2);

  if(normalize)
  {
  hist1c_1->Scale(factor/hist1c_1->Integral());
  hist2c_1->Scale(factor/hist2c_1->Integral());
  hist3c_1->Scale(factor/hist3c_1->Integral());
  hist4c_1->Scale(factor/hist4c_1->Integral());
  hist5c_1->Scale(factor/hist5c_1->Integral());
  hist6c_1->Scale(factor/hist6c_1->Integral());
  hist7c_1->Scale(factor/hist7c_1->Integral());
  hist8c_1->Scale(factor/hist8c_1->Integral());
  hist9c_1->Scale(factor/hist9c_1->Integral());

  hist1c_2->Scale(factor/hist1c_2->Integral());
  hist2c_2->Scale(factor/hist2c_2->Integral());
  hist3c_2->Scale(factor/hist3c_2->Integral());
  hist4c_2->Scale(factor/hist4c_2->Integral());
  hist5c_2->Scale(factor/hist5c_2->Integral());
  hist6c_2->Scale(factor/hist6c_2->Integral());
  hist7c_2->Scale(factor/hist7c_2->Integral());
  hist8c_2->Scale(factor/hist8c_2->Integral());
  hist9c_2->Scale(factor/hist9c_2->Integral());
  }

  c1->cd(1);
  hist1c_2->SetLineColor(kBlack);
  hist1c_1->SetLineColor(kRed);
  hs1->Draw("hist, nostack");  
  
  c1->cd(2);
  hist2c_2->SetLineColor(kBlack);
  hist2c_1->SetLineColor(kRed);
  hs2->Draw("hist, nostack");  


  c1->cd(3);
  hist3c_2->SetLineColor(kBlack);
  hist3c_1->SetLineColor(kRed);
  hs3->Draw("hist, nostack");  

  c1->cd(4);
  hist4c_2->SetLineColor(kBlack);
  hist4c_1->SetLineColor(kRed);
  hs4->Draw("hist, nostack");  


  c1->cd(5);
  hist5c_1->SetLineColor(kRed);
  hist5c_2->SetLineColor(kBlack);
  hs5->Draw("hist, nostack");  


  c1->cd(6);
  hist6c_2->SetLineColor(kBlack);
  hist6c_1->SetLineColor(kRed);
  hs6->Draw("hist, nostack");  

  
  c1->cd(7);
  hist7c_1->SetLineColor(kRed);
  hist7c_2->SetLineColor(kBlack);
  hs7->Draw("hist, nostack");  

  
  c1->cd(8);
  hist8c_2->SetLineColor(kBlack);
  hist8c_1->SetLineColor(kRed);
  hs8->Draw("hist, nostack");  

  
  c1->cd(9);
  hist9c_2->SetLineColor(kBlack);
  hist9c_1->SetLineColor(kRed);
  hs9->Draw("hist, nostack");  

  c1->Print(Folder + output);


  delete hist1c_1;
  delete hist2c_1;
  delete hist3c_1;
  delete hist4c_1;
  delete hist5c_1;
  delete hist6c_1;
  delete hist7c_1;
  delete hist8c_1;
  delete hist9c_1;

  delete hist1c_2;
  delete hist2c_2;
  delete hist3c_2;
  delete hist4c_2;
  delete hist5c_2;
  delete hist6c_2;
  delete hist7c_2;
  delete hist8c_2;
  delete hist9c_2;
  
  delete hs1;
  delete hs2;
  delete hs3;
  delete hs4;
  delete hs5;
  delete hs6;
  delete hs7;
  delete hs8;
  delete hs9;

  delete c1;
  delete tree1;
  delete tree2;
  input2->Close();
  input3->Close();
  delete input2;
  delete input3;
}

