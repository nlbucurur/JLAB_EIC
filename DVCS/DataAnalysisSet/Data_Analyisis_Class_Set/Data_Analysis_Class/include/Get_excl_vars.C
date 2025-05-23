void BDT::Get_excl_vars(TCut cut)
{
   TCanvas *canvas2 = new TCanvas("canvas2", "Exclusivity variables");
   canvas2->Divide(3,2);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

  TChain *pDVCS= new TChain("pDVCS");
  pDVCS->Add(Folder + TData);

  TChain *cont= new TChain("pDVCS");
  cont->Add(Folder + TString("TMostafa_pi0.root"));

  TChain *simu= new TChain("pDVCS");
  simu->Add(Folder + TDVCS);

  auto hs1 = new THStack("hs1","M^{2}_{e#gamma p X} (GeV^{2})");
  auto hs2 = new THStack("hs2","M^{2}_{e#gamma X} (GeV^{2})");
  auto hs3 = new THStack("hs3","#Delta t (GeV^{2})");
  auto hs4 = new THStack("hs4","#Delta #phi");
  auto hs5 = new THStack("hs5","#theta_{#gamma X} (deg)");
  auto hs6 = new THStack("hs6","P_{miss} (GeV)");

  TH1F *hist1c_Data = new TH1F("hist1c_Data","",100,-0.02,0.02);
  TH1F *hist2c_Data = new TH1F("hist2c_Data","",100,0,2);
  TH1F *hist3c_Data = new TH1F("hist3c_Data","",100,-0.5,0.5);
  TH1F *hist4c_Data = new TH1F("hist4c_Data","",100,-1,1);
  TH1F *hist5c_Data = new TH1F("hist5c_Data","",100, 0,2);
  TH1F *hist6c_Data = new TH1F("hist6c_Data","",100, 0,1);

  TH1F *hist1c_Cont = new TH1F("hist1c_Cont","",100,-0.02,0.02);
  TH1F *hist2c_Cont = new TH1F("hist2c_Cont","",100, 0,2);
  TH1F *hist3c_Cont = new TH1F("hist3c_Cont","",100,-0.5,0.5);
  TH1F *hist4c_Cont = new TH1F("hist4c_Cont","",100,-1,1);
  TH1F *hist5c_Cont = new TH1F("hist5c_Cont","",100, 0,2);
  TH1F *hist6c_Cont = new TH1F("hist6c_Cont","",100, 0,1);

  TH1F *hist1c_Sim = new TH1F("hist1c_Sim","",100,-0.02,0.02);
  TH1F *hist2c_Sim = new TH1F("hist2c_Sim","",100,0,2);
  TH1F *hist3c_Sim = new TH1F("hist3c_Sim","",100,-0.5,0.5);
  TH1F *hist4c_Sim = new TH1F("hist4c_Sim","",100,-1,1);
  TH1F *hist5c_Sim = new TH1F("hist5c_Sim","",100, 0,2);
  TH1F *hist6c_Sim = new TH1F("hist6c_Sim","",100, 0,1);

  //Variables to plot
  const char *p1="mm2_eNg";
  const char *p2="mm2_eg";
  const char *p3="delta_t";
  const char *p4="delta_Phi";
  const char *p5="theta_gamma_X";
  const char *p6="miss_mom_eNg";

  cont->Project("hist1c_Cont", p1, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value)))* TCut("Weight"));
  cont->Project("hist2c_Cont", p2, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))) * TCut("Weight"));
  cont->Project("hist3c_Cont", p3, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))) * TCut("Weight"));
  cont->Project("hist4c_Cont", p4, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))) * TCut("Weight"));
  cont->Project("hist5c_Cont", p5, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))) * TCut("Weight"));
  cont->Project("hist6c_Cont", p6, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))) * TCut("Weight"));

  pDVCS->Project("hist1c_Data", p1, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  pDVCS->Project("hist2c_Data", p2, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  pDVCS->Project("hist3c_Data", p3, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  pDVCS->Project("hist4c_Data", p4, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  pDVCS->Project("hist5c_Data", p5, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  pDVCS->Project("hist6c_Data", p6, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));

  simu->Project("hist1c_Sim", p1, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  simu->Project("hist2c_Sim", p2, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  simu->Project("hist3c_Sim", p3, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  simu->Project("hist4c_Sim", p4, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  simu->Project("hist5c_Sim", p5, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  simu->Project("hist6c_Sim", p6, (cut + cut_ref + TCut(Form("_strip_Nuc_BDT > %f", BDT_value))));
  
  hist1c_Data->Add(hist1c_Data, hist1c_Cont, 1, -1);
  hist2c_Data->Add(hist2c_Data, hist2c_Cont, 1, -1);
  hist3c_Data->Add(hist3c_Data, hist3c_Cont, 1, -1);
  hist4c_Data->Add(hist4c_Data, hist4c_Cont, 1, -1);
  hist5c_Data->Add(hist5c_Data, hist5c_Cont, 1, -1);
  hist6c_Data->Add(hist6c_Data, hist6c_Cont, 1, -1);

  hs1->Add(hist1c_Sim);
  hs1->Add(hist1c_Data);
  hs2->Add(hist2c_Sim);
  hs2->Add(hist2c_Data);
  hs3->Add(hist3c_Sim);
  hs3->Add(hist3c_Data);
  hs4->Add(hist4c_Sim);
  hs4->Add(hist4c_Data);
  hs5->Add(hist5c_Sim);
  hs5->Add(hist5c_Data);
  hs6->Add(hist6c_Sim);
  hs6->Add(hist6c_Data);

  hs1->SetHistogram(new TH1F("hstot1","",100,-0.02,0.02));
  hs1->GetHistogram()->GetXaxis()->SetTitle(hs1->GetTitle());
  hs1->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs1->GetHistogram()->SetTitle("");
  hs1->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs1->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs1->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs1->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs1->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs1->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs1->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs1->GetHistogram()->GetXaxis()->SetNdivisions(4);
  hs1->GetHistogram()->GetYaxis()->SetMaxDigits(2);
  
  hs2->SetHistogram(new TH1F("hstot2","",100,0,2));
  hs2->GetHistogram()->GetXaxis()->SetTitle(hs2->GetTitle());
  hs2->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs2->GetHistogram()->SetTitle("");
  hs2->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs2->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs2->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs2->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs2->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs2->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs2->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs2->GetHistogram()->GetXaxis()->SetNdivisions(8);
  hs2->GetHistogram()->GetYaxis()->SetMaxDigits(2);

  hs3->SetHistogram(new TH1F("hstot3","",100,-0.5,0.5));
  hs3->GetHistogram()->GetXaxis()->SetTitle(hs3->GetTitle());
  hs3->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs3->GetHistogram()->SetTitle("");
  hs3->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs3->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs3->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs3->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs3->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs3->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs3->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs3->GetHistogram()->GetXaxis()->SetNdivisions(8);
  hs3->GetHistogram()->GetYaxis()->SetMaxDigits(2);

  hs4->SetHistogram(new TH1F("hstot4","",100,-1,1));
  hs4->GetHistogram()->GetXaxis()->SetTitle(hs4->GetTitle());
  hs4->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs4->GetHistogram()->SetTitle("");
  hs4->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs4->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs4->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs4->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs4->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs4->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs4->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs4->GetHistogram()->GetXaxis()->SetNdivisions(8);
  hs4->GetHistogram()->GetYaxis()->SetMaxDigits(2);

  hs5->SetHistogram(new TH1F("hstot5","",100, 0,2));
  hs5->GetHistogram()->GetXaxis()->SetTitle(hs5->GetTitle());
  hs5->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs5->GetHistogram()->SetTitle("");
  hs5->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs5->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs5->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs5->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs5->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs5->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs5->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs5->GetHistogram()->GetXaxis()->SetNdivisions(8);
  hs5->GetHistogram()->GetYaxis()->SetMaxDigits(2);

  hs6->SetHistogram(new TH1F("hstot6","",100, 0,1));
  hs6->GetHistogram()->GetXaxis()->SetTitle(hs6->GetTitle());
  hs6->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
  hs6->GetHistogram()->SetTitle("");
  hs6->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  hs6->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  hs6->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  hs6->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  hs6->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
  hs6->GetHistogram()->GetXaxis()->SetTitleOffset(0.7);
  hs6->GetHistogram()->GetYaxis()->SetNdivisions(8);
  hs6->GetHistogram()->GetXaxis()->SetNdivisions(8);
  hs6->GetHistogram()->GetYaxis()->SetMaxDigits(2);
  
  canvas2->cd(1);
  hist1c_Sim->Scale(hist1c_Data->Integral()/hist1c_Sim->Integral());
  hist1c_Sim->SetLineColor(kBlack);
  hist1c_Data->SetLineColor(kRed);
  hs1->Draw("nostack,hist");
    
  canvas2->cd(2);
  hist2c_Sim->Scale(hist2c_Data->Integral()/hist2c_Sim->Integral());
  hist2c_Sim->SetLineColor(kBlack);
  hist2c_Data->SetLineColor(kRed);
  hs2->Draw("nostack,hist");


  canvas2->cd(3);
  hist3c_Sim->Scale(hist3c_Data->Integral()/hist3c_Sim->Integral());
  hist3c_Sim->SetLineColor(kBlack);
  hist3c_Data->SetLineColor(kRed);
  hs3->Draw("nostack,hist");

   

  canvas2->cd(4);
  hist4c_Sim->Scale(hist4c_Data->Integral()/hist4c_Sim->Integral());
  hist4c_Sim->SetLineColor(kBlack);
  hist4c_Data->SetLineColor(kRed);
  hs4->Draw("nostack,hist");

  
  canvas2->cd(5);
  hist5c_Data->SetLineColor(kRed);
  hist5c_Sim->Scale(hist5c_Data->Integral()/hist5c_Sim->Integral());
  hist5c_Sim->SetLineColor(kBlack);
  hs5->Draw("nostack,hist");

  canvas2->cd(6);
  hist6c_Data->SetLineColor(kRed);
  hist6c_Sim->Scale(hist6c_Data->Integral()/hist6c_Sim->Integral());
  hist6c_Sim->SetLineColor(kBlack);
  hs6->Draw("nostack,hist");

  canvas2->Print(Folder + TString("Excl_vars.pdf"));
  
  
  delete hist1c_Sim;
  delete hist1c_Data;
  delete hist1c_Cont;
  delete hist2c_Sim;
  delete hist2c_Data;
  delete hist2c_Cont;
  delete hist3c_Sim;
  delete hist3c_Data;
  delete hist3c_Cont;
  delete hist4c_Sim;
  delete hist4c_Data;
  delete hist4c_Cont;
  delete hist5c_Sim;
  delete hist5c_Data;
  delete hist5c_Cont;
  delete hist6c_Sim;
  delete hist6c_Data;
  delete hist6c_Cont;

  delete canvas2;
  delete hs1;
  delete hs2;
  delete hs3;
  delete hs4;
  delete hs5;
  delete hs6;

  delete pDVCS;
  delete simu;
  delete cont;
}
