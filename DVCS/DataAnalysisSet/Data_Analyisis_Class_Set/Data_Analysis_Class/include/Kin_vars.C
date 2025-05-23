void BDT::Kin_vars(TString Data1, TString Data2, TString Data3, TCut cut){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","Histograms");
  c1->Divide(3,3);

  TFile *input3 = new TFile(Folder + Data3,"READ");
  TTree *tree3  = (TTree*)input3->Get("pDVCS");

  TFile *input2 = new TFile(Folder + Data2,"READ");
  TTree *tree2  = (TTree*)input2->Get("pDVCS");

  TFile *input1 = new TFile(Folder + Data1,"READ");
  TTree *tree1  = (TTree*)input1->Get("pDVCS");
  

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

  TH1F* hist1c_3;
  TH1F* hist2c_3;
  TH1F* hist3c_3;
  TH1F* hist4c_3;
  TH1F* hist5c_3;
  TH1F* hist6c_3;
  TH1F* hist7c_3;
  TH1F* hist8c_3;
  TH1F* hist9c_3;

  const char *p1;
  const char *p2;
  const char *p3;
  const char *p4;
  const char *p5;
  const char *p6;
  const char *p7;
  const char *p8;
  const char *p9;


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

  hist1c_3 = new TH1F("hist1c_3","",100,2,10);
  hist2c_3 = new TH1F("hist2c_3","",100,0,7);
  hist3c_3 = new TH1F("hist3c_3","",100,0,3.0);
  hist4c_3 = new TH1F("hist4c_3","",100,0,35);
  hist5c_3 = new TH1F("hist5c_3","",100,4,35);
  hist6c_3 = new TH1F("hist6c_3","",100,4,70);  
  hist7c_3 = new TH1F("hist7c_3","",100,-180,180);
  hist8c_3 = new TH1F("hist8c_3","",100,-180,180);
  hist9c_3 = new TH1F("hist9c_3","",100,-180,180);

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


  tree3->Project("hist1c_3", p1, cut);
  tree3->Project("hist2c_3", p2, cut);
  tree3->Project("hist3c_3", p3, cut);
  tree3->Project("hist4c_3", p4, cut);
  tree3->Project("hist5c_3", p5, cut);
  tree3->Project("hist6c_3", p6, cut);
  tree3->Project("hist7c_3", p7, cut);
  tree3->Project("hist8c_3", p8, cut);
  tree3->Project("hist9c_3", p9, cut);

  tree2->Project("hist1c_2", p1, cut);
  tree2->Project("hist2c_2", p2, cut);
  tree2->Project("hist3c_2", p3, cut);
  tree2->Project("hist4c_2", p4, cut);
  tree2->Project("hist5c_2", p5, cut);
  tree2->Project("hist6c_2", p6, cut);
  tree2->Project("hist7c_2", p7, cut);
  tree2->Project("hist8c_2", p8, cut);
  tree2->Project("hist9c_2", p9, cut);

  tree1->Project("hist1c_1", p1, cut);
  tree1->Project("hist2c_1", p2, cut);
  tree1->Project("hist3c_1", p3, cut);
  tree1->Project("hist4c_1", p4, cut);
  tree1->Project("hist5c_1", p5, cut);
  tree1->Project("hist6c_1", p6, cut);
  tree1->Project("hist7c_1", p7, cut);
  tree1->Project("hist8c_1", p8, cut);
  tree1->Project("hist9c_1", p9, cut);

  hist1c_1->SetLineColor(kRed);
  hist2c_1->SetLineColor(kRed);
  hist3c_1->SetLineColor(kRed);
  hist4c_1->SetLineColor(kRed);
  hist5c_1->SetLineColor(kRed);
  hist6c_1->SetLineColor(kRed);
  hist7c_1->SetLineColor(kRed);
  hist8c_1->SetLineColor(kRed);
  hist9c_1->SetLineColor(kRed);
  hist1c_2->SetLineColor(kBlack);
  hist2c_2->SetLineColor(kBlack);
  hist3c_2->SetLineColor(kBlack);
  hist4c_2->SetLineColor(kBlack);
  hist5c_2->SetLineColor(kBlack);
  hist6c_2->SetLineColor(kBlack);
  hist7c_2->SetLineColor(kBlack);
  hist8c_2->SetLineColor(kBlack);
  hist9c_2->SetLineColor(kBlack);
  hist1c_3->SetLineColor(kBlue);
  hist2c_3->SetLineColor(kBlue);
  hist3c_3->SetLineColor(kBlue);
  hist4c_3->SetLineColor(kBlue);
  hist5c_3->SetLineColor(kBlue);
  hist6c_3->SetLineColor(kBlue);
  hist7c_3->SetLineColor(kBlue);
  hist8c_3->SetLineColor(kBlue);
  hist9c_3->SetLineColor(kBlue);

  hs1->Add(hist1c_1);
  hs1->Add(hist1c_2);
  hs1->Add(hist1c_3);
  hs2->Add(hist2c_1);
  hs2->Add(hist2c_2);
  hs2->Add(hist2c_3);
  hs3->Add(hist3c_1);
  hs3->Add(hist3c_2);
  hs3->Add(hist3c_3);
  hs4->Add(hist4c_1);
  hs4->Add(hist4c_2);
  hs4->Add(hist4c_3);
  hs5->Add(hist5c_1);
  hs5->Add(hist5c_2);
  hs5->Add(hist5c_3);
  hs6->Add(hist6c_1);
  hs6->Add(hist6c_2);
  hs6->Add(hist6c_3);
  hs7->Add(hist7c_1);
  hs7->Add(hist7c_2);
  hs7->Add(hist7c_3);
  hs8->Add(hist8c_1);
  hs8->Add(hist8c_2);
  hs8->Add(hist8c_3);
  hs9->Add(hist9c_1);
  hs9->Add(hist9c_2);
  hs9->Add(hist9c_3);

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

  hist1c_3->Scale(factor/hist1c_3->Integral());
  hist2c_3->Scale(factor/hist2c_3->Integral());
  hist3c_3->Scale(factor/hist3c_3->Integral());
  hist4c_3->Scale(factor/hist4c_3->Integral());
  hist5c_3->Scale(factor/hist5c_3->Integral());
  hist6c_3->Scale(factor/hist6c_3->Integral());
  hist7c_3->Scale(factor/hist7c_3->Integral());
  hist8c_3->Scale(factor/hist8c_3->Integral());
  hist9c_3->Scale(factor/hist9c_3->Integral());

  c1->cd(1);
  hs1->Draw("hist, nostack");  
  
  c1->cd(2);
  hs2->Draw("hist, nostack");  

  c1->cd(3);
  hs3->Draw("hist, nostack");  
  auto legend = new TLegend(0.6,0.6,0.88,0.88);
  legend->AddEntry(hist1c_1,"Data","lp");
  legend->AddEntry(hist1c_2,"MC DVCS","lp");
  legend->AddEntry(hist1c_3,"MC Pi0","lp");
  legend->Draw();

  c1->cd(4);
  hs4->Draw("hist, nostack");  

  c1->cd(5);
  hs5->Draw("hist, nostack");  

  c1->cd(6);
  hs6->Draw("hist, nostack");  
  
  c1->cd(7);
  hs7->Draw("hist, nostack");  
  
  c1->cd(8);
  hs8->Draw("hist, nostack");  
  
  c1->cd(9);
  hs9->Draw("hist, nostack");  

  double axtit_size = 0.06;
  double axlab_size = 0.04;
  double axXtit_offs = 0.7;
  double axYtit_offs = 0.7;

  hs1->GetXaxis()->SetTitle(hs1->GetTitle());
  hs1->GetYaxis()->SetTitle("counts/total events");
  hs1->SetTitle("Photon");
  hs1->GetXaxis()->SetTitleSize(axtit_size);
  hs1->GetYaxis()->SetTitleSize(axtit_size);
  hs1->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs1->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs1->GetXaxis()->SetLabelSize(axlab_size);
  hs1->GetYaxis()->SetLabelSize(axlab_size);
  hs1->GetYaxis()->SetNdivisions(6);
  hs1->GetYaxis()->SetMaxDigits(1);

  hs2->GetXaxis()->SetTitle(hs2->GetTitle());
  hs2->GetYaxis()->SetTitle("counts/total events");
  hs2->SetTitle("Electron");
  hs2->GetXaxis()->SetTitleSize(axtit_size);
  hs2->GetYaxis()->SetTitleSize(axtit_size);
  hs2->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs2->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs2->GetXaxis()->SetLabelSize(axlab_size);
  hs2->GetYaxis()->SetLabelSize(axlab_size);
  hs2->GetYaxis()->SetNdivisions(6);
  hs2->GetYaxis()->SetMaxDigits(1);

  hs3->GetXaxis()->SetTitle(hs3->GetTitle());
  hs3->GetYaxis()->SetTitle("counts/total events");
  hs3->SetTitle("Proton");
  hs3->GetXaxis()->SetTitleSize(axtit_size);
  hs3->GetYaxis()->SetTitleSize(axtit_size);
  hs3->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs3->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs3->GetXaxis()->SetLabelSize(axlab_size);
  hs3->GetYaxis()->SetLabelSize(axlab_size);
  hs3->GetYaxis()->SetNdivisions(6);
  hs3->GetYaxis()->SetMaxDigits(1);

  hs4->GetXaxis()->SetTitle(hs4->GetTitle());
  hs4->GetYaxis()->SetTitle("counts/total events");
  hs4->SetTitle("Photon");
  hs4->GetXaxis()->SetTitleSize(axtit_size);
  hs4->GetYaxis()->SetTitleSize(axtit_size);
  hs4->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs4->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs4->GetXaxis()->SetLabelSize(axlab_size);
  hs4->GetYaxis()->SetLabelSize(axlab_size);
  hs4->GetYaxis()->SetNdivisions(6);
  hs4->GetYaxis()->SetMaxDigits(1);

  hs5->GetXaxis()->SetTitle(hs5->GetTitle());
  hs5->GetYaxis()->SetTitle("counts/total events");
  hs5->SetTitle("Electron");
  hs5->GetXaxis()->SetTitleSize(axtit_size);
  hs5->GetYaxis()->SetTitleSize(axtit_size);
  hs5->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs5->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs5->GetXaxis()->SetLabelSize(axlab_size);
  hs5->GetYaxis()->SetLabelSize(axlab_size);
  hs5->GetYaxis()->SetNdivisions(6);
  hs5->GetYaxis()->SetMaxDigits(1);

  hs6->GetXaxis()->SetTitle(hs6->GetTitle());
  hs6->GetYaxis()->SetTitle("counts/total events");
  hs6->SetTitle("Proton");
  hs6->GetXaxis()->SetTitleSize(axtit_size);
  hs6->GetYaxis()->SetTitleSize(axtit_size);
  hs6->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs6->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs6->GetXaxis()->SetLabelSize(axlab_size);
  hs6->GetYaxis()->SetLabelSize(axlab_size);
  hs6->GetYaxis()->SetNdivisions(6);
  hs6->GetYaxis()->SetMaxDigits(1);

  hs7->GetXaxis()->SetTitle(hs7->GetTitle());
  hs7->GetYaxis()->SetTitle("counts/total events");
  hs7->SetTitle("Photon");
  hs7->GetXaxis()->SetTitleSize(axtit_size);
  hs7->GetYaxis()->SetTitleSize(axtit_size);
  hs7->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs7->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs7->GetXaxis()->SetLabelSize(axlab_size);
  hs7->GetYaxis()->SetLabelSize(axlab_size);
  hs7->GetYaxis()->SetNdivisions(6);
  hs7->GetYaxis()->SetMaxDigits(1);

  hs8->GetXaxis()->SetTitle(hs8->GetTitle());
  hs8->GetYaxis()->SetTitle("counts/total events");
  hs8->SetTitle("Electron");
  hs8->GetXaxis()->SetTitleSize(axtit_size);
  hs8->GetYaxis()->SetTitleSize(axtit_size);
  hs8->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs8->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs8->GetXaxis()->SetLabelSize(axlab_size);
  hs8->GetYaxis()->SetLabelSize(axlab_size);
  hs8->GetYaxis()->SetNdivisions(6);
  hs8->GetYaxis()->SetMaxDigits(1);

  hs9->GetXaxis()->SetTitle(hs9->GetTitle());
  hs9->GetYaxis()->SetTitle("counts/total events");
  hs9->SetTitle("Proton");
  hs9->GetXaxis()->SetTitleSize(axtit_size);
  hs9->GetYaxis()->SetTitleSize(axtit_size);
  hs9->GetXaxis()->SetTitleOffset(axXtit_offs);
  hs9->GetYaxis()->SetTitleOffset(axYtit_offs);
  hs9->GetXaxis()->SetLabelSize(axlab_size);
  hs9->GetYaxis()->SetLabelSize(axlab_size);
  hs9->GetYaxis()->SetNdivisions(6);
  hs9->GetYaxis()->SetMaxDigits(1);

  c1->Print(Folder + TString("Kin_Vars.pdf"));


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
  
  delete hist1c_3;
  delete hist2c_3;
  delete hist3c_3;
  delete hist4c_3;
  delete hist5c_3;
  delete hist6c_3;
  delete hist7c_3;
  delete hist8c_3;
  delete hist9c_3;

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
  delete tree3;
  input1->Close();
  input2->Close();
  input3->Close();
  delete input1;
  delete input2;
  delete input3;
}

