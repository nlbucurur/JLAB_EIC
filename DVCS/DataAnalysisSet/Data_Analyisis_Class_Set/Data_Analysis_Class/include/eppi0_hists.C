void BDT::eppi0_hists(){
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","Histograms",1000,500);
  c1->Divide(2);

  TFile *input2 = new TFile( Folder + eppi0_name,"READ");
  TTree *pDVCS_Data = (TTree*)input2->Get("eppi0");
  
  TFile *input3 = new TFile(sim_eppi0,"READ");
  TTree *pDVCS_Sim = (TTree*)input3->Get("eppi0");

  double factor=1;
  
  TCut cut =TCut("bestCandidateFlag==1 && abs(mm2_eNgg)<0.01 && strip_Pi0_2DChi2 < 0.04 && strip_Pi0_E > 2 && strip_Pi0_E<9 && strip_Pi0_Theta > 5 && strip_Pi0_Theta < 33 && N_Ph<3 && abs(t_Pi0/strip_Q2)<1");//&& strip_Pi0_IM2 > 0.012 && strip_Pi0_IM2 < 0.022");

  //TCut cut =TCut("bestCandidateFlag==1 && strip_Q2 > 1.0 && strip_W > 2 && strip_Nuc_P > 0.35 && strip_El_P > 1.0  && strip_El_vz < 10 && strip_El_vz > -12 && TMath::Abs(Phi_Nuc - Phi_Pi0) < 2 && TMath::Abs(t_Nuc - t_Pi0) < 2 && TMath::Sqrt(Xbal * Xbal + Ybal*Ybal + Zbal*Zbal) <1");
  
  TH1F *hist3c_Sim = new TH1F("hist3c_Sim","Histogram",100,0,3);
  TH1F *hist9c_Sim = new TH1F("hist9c_Sim","Histogram",100,0.0,0.03);

  TH1F *hist3c_Data = new TH1F("hist3c_Data","Histogram",100,0,3);
  TH1F *hist9c_Data = new TH1F("hist9c_Data","Histogram",100,0.0,0.03);


  //Variables to plot

  const char *p3="mm2_egg";
  const char *p9="strip_Pi0_IM2";
  
  pDVCS_Sim->Project("hist3c_Sim", p3, cut);
  pDVCS_Sim->Project("hist9c_Sim", p9, cut);

  pDVCS_Data->Project("hist3c_Data", p3, cut);
  pDVCS_Data->Project("hist9c_Data", p9, cut);

  auto hs1 = new THStack("hs1","Missing mass of ep #rightarrow e#gamma#gamma (GeV^{2})");
  auto hs2 = new THStack("hs2","M_{#pi^{0}}^{2} (GeV^{2})");

  c1->cd(1);
  hist3c_Sim->SetTitle("Missing mass of ep #rightarrow e#gamma#gamma (GeV^{2})");
  hist3c_Data->SetTitle("Missing mass of ep #rightarrow e#gamma#gamma (GeV^{2})");
  hist3c_Sim->Scale(factor/hist3c_Sim->GetEntries());
  hist3c_Sim->SetLineColor(kBlack);
  hist3c_Data->Scale(factor/hist3c_Data->GetEntries());
  hist3c_Data->SetLineColor(kRed);

  hs1->Add(hist3c_Sim);
  hs1->Add(hist3c_Data);
  hs1->Draw("nostack,hist");


  
  c1->cd(2);
  hist9c_Sim->SetTitle("M_{#pi^{0}}^{2} (GeV^{2})");
  hist9c_Data->SetTitle("M_{#pi^{0}}^{2} (GeV^{2})");
  hist9c_Sim->Scale(factor/hist9c_Sim->GetEntries());
  hist9c_Sim->SetLineColor(kBlack);
  hist9c_Data->Scale(factor/hist9c_Data->GetEntries());
  hist9c_Data->SetLineColor(kRed);

  hs2->Add(hist9c_Sim);
  hs2->Add(hist9c_Data);
  hs2->Draw("nostack,hist");

  
  // input->Close();
  c1->Print(Folder + TString("eppi0_Data_vs_MC.pdf"));

  delete hs1;
  delete hs2;

}

