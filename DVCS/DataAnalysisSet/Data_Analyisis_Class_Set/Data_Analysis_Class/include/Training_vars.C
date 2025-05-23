void BDT::Training_vars(TString Data, TString MC_DVCS, TString MC_Pi0, TCut cut){
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","Histograms");
  c1->Divide(3,2);

  TFile *input2 = new TFile(Data,"READ");
  TTree *pDVCS_Data = (TTree*)input2->Get("pDVCS");
  
  TFile *input3 = new TFile(MC_DVCS,"READ");
  TTree *pDVCS_Sim = (TTree*)input3->Get("pDVCS");

  TFile *input4 = new TFile(MC_Pi0,"READ");
  TTree *pDVCS_Pi0;
  if( input4->GetListOfKeys()->Contains("pDVCS") )
    pDVCS_Pi0     = (TTree*)input4->Get("pDVCS");
  else
    pDVCS_Pi0     = (TTree*)input4->Get("eppi0");

  TString t1="M^{2}_{e#gamma p X} (GeV^{2})";
  TString t2="M^{2}_{e#gamma X} (GeV^{2})";
  TString t3="#Delta t (GeV^{2})";
  TString t4="#Delta #phi (deg)";
  TString t5="#theta_{#gamma X} (deg)";
  TString t6="";
  TString t7="";
  TString t8="";
  TString t9="";

  auto hs1 = new THStack("hs1",t1);
  auto hs2 = new THStack("hs2",t2);
  auto hs3 = new THStack("hs3",t3);
  auto hs4 = new THStack("hs4",t4);
  auto hs5 = new THStack("hs5",t5);
  auto hs6 = new THStack("hs6",t6);
  auto hs7 = new THStack("hs7",t7);
  auto hs8 = new THStack("hs8",t8);
  auto hs9 = new THStack("hs9",t9);

  TH1F *hist1c_Sim = new TH1F("hist1c_Sim","",100,-0.02,0.02);
  TH1F *hist2c_Sim = new TH1F("hist2c_Sim","",100, 0,3);
  TH1F *hist3c_Sim = new TH1F("hist3c_Sim","",100,-0.5,0.5);
  TH1F *hist4c_Sim = new TH1F("hist4c_Sim","",100,-0.2,0.4);
  TH1F *hist5c_Sim = new TH1F("hist5c_Sim","",100, -0.2,2);
  
  TH1F *hist6c_Sim = new TH1F("hist6c_Sim","",100,-0.2,0.4);  
  TH1F *hist7c_Sim = new TH1F("hist7c_Sim","",100,0,100);
  TH1F *hist8c_Sim = new TH1F("hist8c_Sim","",100,0,2);
  TH1F *hist9c_Sim = new TH1F("hist9c_Sim","",100,0,100);
  
  TH1F *hist1c_Data = new TH1F("hist1c_Data","",100,-0.02,0.02);
  TH1F *hist2c_Data = new TH1F("hist2c_Data","",100,0,3);
  TH1F *hist3c_Data = new TH1F("hist3c_Data","",100,-0.5,0.5);
  TH1F *hist4c_Data = new TH1F("hist4c_Data","",100,-0.2,0.4);
  TH1F *hist5c_Data = new TH1F("hist5c_Data","",100, -0.2,2);

  TH1F *hist6c_Data = new TH1F("hist6c_Data","",100,-0.2,0.4);  
  TH1F *hist7c_Data = new TH1F("hist7c_Data","",100,0,100);
  TH1F *hist8c_Data = new TH1F("hist8c_Data","",100,0,2);
  TH1F *hist9c_Data = new TH1F("hist9c_Data","",100,0,100);

  TH1F *hist1c_Pi0 = new TH1F("hist1c_Pi0","",100,-0.02,0.02);
  TH1F *hist2c_Pi0 = new TH1F("hist2c_Pi0","",100,0,3);
  TH1F *hist3c_Pi0 = new TH1F("hist3c_Pi0","",100,-0.5,0.5);
  TH1F *hist4c_Pi0 = new TH1F("hist4c_Pi0","",100,-0.2,0.4);
  TH1F *hist5c_Pi0 = new TH1F("hist5c_Pi0","",100, -0.2,2);

  TH1F *hist6c_Pi0 = new TH1F("hist6c_Pi0","",100,-0.2,0.4);  
  TH1F *hist7c_Pi0 = new TH1F("hist7c_Pi0","",100,0,100);
  TH1F *hist8c_Pi0 = new TH1F("hist8c_Pi0","",100,0,2);
  TH1F *hist9c_Pi0 = new TH1F("hist9c_Pi0","",100,0,100);

  //Variables to plot
  const char *p1="mm2_eNg";
  const char *p2="mm2_eg";
  const char *p3="delta_t";
  const char *p4="delta_Phi";
  const char *p5="theta_gamma_X";

  const char *p6="delta_t";
  const char *p7="strip_El_Theta";
  const char *p8="theta_gamma_X";
  const char *p9="strip_Ph_Theta";
  
  pDVCS_Sim->Project("hist1c_Sim", p1, cut);
  pDVCS_Sim->Project("hist2c_Sim", p2, cut);
  pDVCS_Sim->Project("hist3c_Sim", p3, cut);
  pDVCS_Sim->Project("hist4c_Sim", p4, cut);
  pDVCS_Sim->Project("hist5c_Sim", p5, cut);
  //pDVCS_Sim->Project("hist6c_Sim", p6, cut);
  //pDVCS_Sim->Project("hist7c_Sim", p7, cut);
  //pDVCS_Sim->Project("hist8c_Sim", p8, cut);
  //pDVCS_Sim->Project("hist9c_Sim", p9, cut);

  pDVCS_Data->Project("hist1c_Data", p1, cut);
  pDVCS_Data->Project("hist2c_Data", p2, cut);
  pDVCS_Data->Project("hist3c_Data", p3, cut);
  pDVCS_Data->Project("hist4c_Data", p4, cut);
  pDVCS_Data->Project("hist5c_Data", p5, cut);
  //pDVCS_Data->Project("hist6c_Data", p6, cut);
  //pDVCS_Data->Project("hist7c_Data", p7, cut);
  //pDVCS_Data->Project("hist8c_Data", p8, cut);
  //pDVCS_Data->Project("hist9c_Data", p9, cut);

  pDVCS_Pi0->Project("hist1c_Pi0", p1, cut);
  pDVCS_Pi0->Project("hist2c_Pi0", p2, cut);
  pDVCS_Pi0->Project("hist3c_Pi0", p3, cut);
  pDVCS_Pi0->Project("hist4c_Pi0", p4, cut);
  pDVCS_Pi0->Project("hist5c_Pi0", p5, cut);
  //pDVCS_Pi0->Project("hist6c_Pi0", p6, cut);
  //pDVCS_Pi0->Project("hist7c_Pi0", p7, cut);
  //pDVCS_Pi0->Project("hist8c_Pi0", p8, cut);
  //pDVCS_Pi0->Project("hist9c_Pi0", p9, cut);

  hs1->Add(hist1c_Sim);
  hs1->Add(hist1c_Data);
  hs1->Add(hist1c_Pi0);

  hs2->Add(hist2c_Sim);
  hs2->Add(hist2c_Data);
  hs2->Add(hist2c_Pi0);

  hs3->Add(hist3c_Sim);
  hs3->Add(hist3c_Data);
  hs3->Add(hist3c_Pi0);

  hs4->Add(hist4c_Sim);
  hs4->Add(hist4c_Data);
  hs4->Add(hist4c_Pi0);

  hs5->Add(hist5c_Sim);
  hs5->Add(hist5c_Data);
  hs5->Add(hist5c_Pi0);

  double axtit_size = 0.06;
  double axlab_size = 0.04;
  double axXtit_offs = 0.7;
  double axYtit_offs = -0.4;

  //All
  double factor1=0.75;
  double factor2=0.4;

  //FD
  //double factor1=0.3;
  //double factor2=0.8;

  //FT
  //double factor1=1.0;
  //double factor2=0.4;

  hist1c_Sim->Scale(hist1c_Data->Integral()*factor1/hist1c_Sim->Integral());
  hist1c_Pi0->Scale(hist1c_Data->Integral()*factor2/hist1c_Pi0->Integral());

  hist2c_Sim->Scale(hist2c_Data->Integral()*factor1/hist2c_Sim->Integral());
  hist2c_Pi0->Scale(hist2c_Data->Integral()*factor2/hist2c_Pi0->Integral());

  hist3c_Sim->Scale(hist3c_Data->Integral()*factor1/hist3c_Sim->Integral());
  hist3c_Pi0->Scale(hist3c_Data->Integral()*factor2/hist3c_Pi0->Integral());

  hist4c_Sim->Scale(hist4c_Data->Integral()*factor1/hist4c_Sim->Integral());
  hist4c_Pi0->Scale(hist4c_Data->Integral()*factor2/hist4c_Pi0->Integral());

  hist5c_Sim->Scale(hist5c_Data->Integral()*factor1/hist5c_Sim->Integral());
  hist5c_Pi0->Scale(hist5c_Data->Integral()*factor2/hist5c_Pi0->Integral());

  c1->cd(1);

  hist1c_Sim->SetLineColor(kBlack);
  hist1c_Data->SetLineColor(kRed);
  hist1c_Pi0->SetLineColor(kBlue);

  hs1->Draw("nostack,hist");
    
  c1->cd(2);
  hist2c_Sim->SetLineColor(kBlack);
  hist2c_Data->SetLineColor(kRed);
  hist2c_Pi0->SetLineColor(kBlue);

  hs2->Draw("nostack,hist");


  c1->cd(3);
  hist3c_Sim->SetLineColor(kBlack);
  hist3c_Data->SetLineColor(kRed);
  hist3c_Pi0->SetLineColor(kBlue);

  hs3->Draw("nostack,hist");

   

  c1->cd(4);
  hist4c_Sim->SetLineColor(kBlack);
  hist4c_Data->SetLineColor(kRed);
  hist4c_Pi0->SetLineColor(kBlue);

  hs4->Draw("nostack,hist");

  
  c1->cd(5);
  hist5c_Data->SetLineColor(kRed);
  hist5c_Sim->SetLineColor(kBlack);
  hist5c_Pi0->SetLineColor(kBlue);

  hs5->Draw("nostack,hist");
  
  c1->cd(6);
  auto legend = new TLegend(0.1,0.1,0.9,0.9);
  legend->AddEntry(hist1c_Data,"Data","lp");
  legend->AddEntry(hist1c_Sim,"MC DVCS","lp");
  legend->AddEntry(hist1c_Pi0,"MC Pi0","lp");
  legend->Draw();

  hs1->GetXaxis()->SetTitle(t1);
  hs1->GetYaxis()->SetTitle("");
  hs1->GetXaxis()->SetTitleSize(axtit_size);
  hs1->GetYaxis()->SetTitleSize(axtit_size);
  hs1->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs1->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs1->GetXaxis()->SetLabelSize(axlab_size);
  hs1->GetYaxis()->SetLabelSize(axlab_size);

  hs2->GetXaxis()->SetTitle(t2);
  hs2->GetYaxis()->SetTitle("");
  hs2->GetXaxis()->SetTitleSize(axtit_size);
  hs2->GetYaxis()->SetTitleSize(axtit_size);
  hs2->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs2->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs2->GetXaxis()->SetLabelSize(axlab_size);
  hs2->GetYaxis()->SetLabelSize(axlab_size);

  hs3->GetXaxis()->SetTitle(t3);
  hs3->GetYaxis()->SetTitle("");
  hs3->GetXaxis()->SetTitleSize(axtit_size);
  hs3->GetYaxis()->SetTitleSize(axtit_size);
  hs3->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs3->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs3->GetXaxis()->SetLabelSize(axlab_size);
  hs3->GetYaxis()->SetLabelSize(axlab_size);

  hs4->GetXaxis()->SetTitle(t4);
  hs4->GetYaxis()->SetTitle("");
  hs4->GetXaxis()->SetTitleSize(axtit_size);
  hs4->GetYaxis()->SetTitleSize(axtit_size);
  hs4->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs4->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs4->GetXaxis()->SetLabelSize(axlab_size);
  hs4->GetYaxis()->SetLabelSize(axlab_size);

  hs5->GetXaxis()->SetTitle(t5);
  hs5->GetYaxis()->SetTitle("");
  hs5->GetXaxis()->SetTitleSize(axtit_size);
  hs5->GetYaxis()->SetTitleSize(axtit_size);
  hs5->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs5->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs5->GetXaxis()->SetLabelSize(axlab_size);
  hs5->GetYaxis()->SetLabelSize(axlab_size);

/*
  c1->cd(6);
  hist6c_Sim->Scale(factor/hist6c_Sim->Integral());
  hist6c_Sim->SetLineColor(kBlack);
  hist6c_Data->Scale(factor/hist6c_Data->Integral());
  hist6c_Data->SetLineColor(kRed);
  hist6c_Pi0->Scale(factor/hist6c_Pi0->Integral());
  hist6c_Pi0->SetLineColor(kBlue);

  hs6->Add(hist6c_Sim);
  hs6->Add(hist6c_Data);
  hs6->Add(hist6c_Pi0);
  hs6->Draw("nostack,hist");
  
  c1->cd(7);
  hist7c_Sim->Scale(factor/hist7c_Sim->Integral());
  hist7c_Sim->SetLineColor(kBlack);
  hist7c_Data->Scale(factor/hist7c_Data->Integral());
  hist7c_Data->SetLineColor(kRed);
  hist7c_Pi0->Scale(factor/hist7c_Pi0->Integral());
  hist7c_Pi0->SetLineColor(kBlue);

  hs7->Add(hist7c_Sim);
  hs7->Add(hist7c_Data);
  hs7->Add(hist7c_Pi0);
  hs7->Draw("nostack,hist");

  c1->cd(8);
  hist8c_Sim->Scale(factor/hist8c_Sim->Integral());
  hist8c_Sim->SetLineColor(kBlack);
  hist8c_Data->Scale(factor/hist8c_Data->Integral());
  hist8c_Data->SetLineColor(kRed);
  hist8c_Pi0->Scale(factor/hist8c_Pi0->Integral());
  hist8c_Pi0->SetLineColor(kBlue);

  hs8->Add(hist8c_Sim);
  hs8->Add(hist8c_Data);
  hs8->Add(hist8c_Pi0);
  hs8->Draw("nostack,hist");

  c1->cd(9);
  hist9c_Sim->Scale(factor/hist9c_Sim->Integral());
  hist9c_Sim->SetLineColor(kBlack);
  hist9c_Data->Scale(factor/hist9c_Data->Integral());
  hist9c_Data->SetLineColor(kRed);
  hist9c_Pi0->Scale(factor/hist9c_Pi0->Integral());
  hist9c_Pi0->SetLineColor(kBlue);

  hs9->Add(hist9c_Sim);
  hs9->Add(hist9c_Data);
  hs9->Add(hist9c_Pi0);
  hs9->Draw("nostack,hist");
  
  hs6->GetXaxis()->SetTitle(t6);
  hs6->GetYaxis()->SetTitle("");
  hs6->GetXaxis()->SetTitleSize(axtit_size);
  hs6->GetYaxis()->SetTitleSize(axtit_size);
  hs6->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs6->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs6->GetXaxis()->SetLabelSize(axlab_size);
  hs6->GetYaxis()->SetLabelSize(axlab_size);

  hs7->GetXaxis()->SetTitle(t7);
  hs7->GetYaxis()->SetTitle("");
  hs7->GetXaxis()->SetTitleSize(axtit_size);
  hs7->GetYaxis()->SetTitleSize(axtit_size);
  hs7->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs7->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs7->GetXaxis()->SetLabelSize(axlab_size);
  hs7->GetYaxis()->SetLabelSize(axlab_size);

  hs8->GetXaxis()->SetTitle(t8);
  hs8->GetYaxis()->SetTitle("");
  hs8->GetXaxis()->SetTitleSize(axtit_size);
  hs8->GetYaxis()->SetTitleSize(axtit_size);
  hs8->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs8->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs8->GetXaxis()->SetLabelSize(axlab_size);
  hs8->GetYaxis()->SetLabelSize(axlab_size);

  hs9->GetXaxis()->SetTitle(t9);
  hs9->GetYaxis()->SetTitle("");
  hs9->GetXaxis()->SetTitleSize(axtit_size);
  hs9->GetYaxis()->SetTitleSize(axtit_size);
  hs9->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs9->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs9->GetXaxis()->SetLabelSize(axlab_size);
  hs9->GetYaxis()->SetLabelSize(axlab_size);
  
  
  
  // input->Close();
  */
  c1->Print(Folder + TString("Training_vars.pdf"));
  
  delete hist1c_Sim;
  delete hist1c_Data;
  delete hist1c_Pi0;
  delete hist2c_Sim;
  delete hist2c_Data;
  delete hist2c_Pi0;
  delete hist3c_Sim;
  delete hist3c_Data;
  delete hist3c_Pi0;
  delete hist4c_Sim;
  delete hist4c_Data;
  delete hist4c_Pi0;
  delete hist5c_Sim;
  delete hist5c_Data;
  delete hist5c_Pi0;
  delete hist6c_Sim;
  delete hist6c_Data;
  delete hist6c_Pi0;
  delete hist7c_Sim;
  delete hist7c_Data;
  delete hist7c_Pi0;
  delete hist8c_Sim;
  delete hist8c_Data;
  delete hist8c_Pi0;
  delete hist9c_Sim;
  delete hist9c_Data;
  delete hist9c_Pi0;
  delete c1;
  delete hs1;
  delete hs2;
  delete hs3;
  delete hs4;
  delete hs5;
  delete hs6;
  delete hs7;
  delete hs8;
  delete hs9;

  delete pDVCS_Data;
  delete pDVCS_Sim;
  delete pDVCS_Pi0;

  input2->Close();
  input3->Close();
  input4->Close();

  delete input2;
  delete input3;
  delete input4;

}

