void BDT::Compare_three(TCut cut, TString Data1, TString Data2){
  gStyle->SetOptStat(0);
  TString Data3 = TString("Tested_1gamma.root");
  TCanvas *c1 = new TCanvas("c1","Histograms");
  c1->Divide(3,3);

  TFile *input2 = new TFile(Folder + Data1,"READ");
  TTree *tree1 = (TTree*)input2->Get("pDVCS");
  
  TFile *input3 = new TFile(Folder + Data2,"READ");
  TTree *tree2 = (TTree*)input3->Get("pDVCS");

  TFile *input4 = new TFile(Folder + Data3,"READ");
  TTree *tree3 = (TTree*)input4->Get("pDVCS");

  double factor=1;
  TH1F *hist1c_3 = new TH1F("hist1c_3","#bf{q}'(GeV)",100,2,10);
  TH1F *hist2c_3 = new TH1F("hist2c_3","#bf{k}'(GeV)",100,0,7);
  TH1F *hist3c_3 = new TH1F("hist3c_3","M^{2}_{ep#gamma X} (GeV^{2})",100,-0.1,0.1);
  TH1F *hist4c_3 = new TH1F("hist4c_3","M^{2}_{e#gamma X} (GeV^{2})",100,0,3);
  TH1F *hist5c_3 = new TH1F("hist5c_3","#Delta #phi (deg)",100,-0.5,0.5);
  TH1F *hist6c_3 = new TH1F("hist6c_3","#Delta t (GeV^{2})",100,-0.2,0.4);  
  TH1F *hist7c_3 = new TH1F("hist7c_3","P_{miss}(ep#gamma) (GeV)",100,0,0.5);
  TH1F *hist8c_3 = new TH1F("hist8c_3","#theta_{#gamma X} (deg)",100,0,2);
  TH1F *hist9c_3 = new TH1F("hist9c_3","P_{#perp} (GeV)",100,0,0.06);

  TH1F *hist1c_2 = new TH1F("hist1c_2","#bf{q}'(GeV)",100,2,10);
  TH1F *hist2c_2 = new TH1F("hist2c_2","#bf{k}'(GeV)",100,0,7);
  TH1F *hist3c_2 = new TH1F("hist3c_2","M^{2}_{ep#gamma X} (GeV^{2})",100,-0.1,0.1);
  TH1F *hist4c_2 = new TH1F("hist4c_2","M^{2}_{e#gamma X} (GeV^{2})",100,0,3);
  TH1F *hist5c_2 = new TH1F("hist5c_2","#Delta #phi (deg)",100,-0.5,0.5);
  TH1F *hist6c_2 = new TH1F("hist6c_2","#Delta t (GeV^{2})",100,-0.2,0.4);  
  TH1F *hist7c_2 = new TH1F("hist7c_2","P_{miss}(ep#gamma) (GeV)",100,0,0.5);
  TH1F *hist8c_2 = new TH1F("hist8c_2","#theta_{#gamma X} (deg)",100,0,2);
  TH1F *hist9c_2 = new TH1F("hist9c_2","P_{#perp} (GeV)",100,0,0.06);

  TH1F *hist1c_1 = new TH1F("hist1c_1","#bf{q}'(GeV)",100,2,10);
  TH1F *hist2c_1 = new TH1F("hist2c_1","#bf{q}'(GeV)",100,0,7);
  TH1F *hist3c_1 = new TH1F("hist3c_1","M^{2}_{ep#gamma X} (GeV^{2})",100,-0.1,0.1);
  TH1F *hist4c_1 = new TH1F("hist4c_1","M^{2}_{e#gamma X} (GeV^{2})",100,0,3);
  TH1F *hist5c_1 = new TH1F("hist5c_1","#Delta #phi (deg)",100,-0.5,0.5);
  TH1F *hist6c_1 = new TH1F("hist6c_1","#Delta t (GeV^{2})",100,-0.2,0.4);  
  TH1F *hist7c_1 = new TH1F("hist7c_1","P_{miss}(ep#gamma) (GeV)",100,0,0.5);
  TH1F *hist8c_1 = new TH1F("hist8c_1","#theta_{#gamma X} (deg)",100,0,2);
  TH1F *hist9c_1 = new TH1F("hist9c_1","P_{#perp} (GeV)",100,0,0.06);

  //Variables to plot
  const char *p1="strip_Ph_P";
  const char *p2="strip_El_P";
  const char *p3="mm2_eNg";
  const char *p4="mm2_eg";
  const char *p5="delta_Phi";
  const char *p6="delta_t";
  const char *p7="miss_mom_eNg";
  const char *p8="theta_gamma_X";
  const char *p9="p_perp";

  tree3->Project("hist1c_3", p1, cut);
  tree3->Project("hist2c_3", p2, cut);
  tree3->Project("hist3c_3", p3, cut);
  tree3->Project("hist4c_3", p4, cut);
  tree3->Project("hist5c_3", p5, cut);
  tree3->Project("hist6c_3", p6, cut);
  tree3->Project("hist7c_3", p7, cut);
  tree3->Project("hist8c_3", p8, cut);
  tree3->Project("hist9c_3", p9, cut);

  tree2->Project("hist1c_2", p1, TCut("Weight")*cut);
  tree2->Project("hist2c_2", p2, TCut("Weight")*cut);
  tree2->Project("hist3c_2", p3, TCut("Weight")*cut);
  tree2->Project("hist4c_2", p4, TCut("Weight")*cut);
  tree2->Project("hist5c_2", p5, TCut("Weight")*cut);
  tree2->Project("hist6c_2", p6, TCut("Weight")*cut);
  tree2->Project("hist7c_2", p7, TCut("Weight")*cut);
  tree2->Project("hist8c_2", p8, TCut("Weight")*cut);
  tree2->Project("hist9c_2", p9, TCut("Weight")*cut);

  tree1->Project("hist1c_1", p1, TCut("Weight")*cut);
  tree1->Project("hist2c_1", p2, TCut("Weight")*cut);
  tree1->Project("hist3c_1", p3, TCut("Weight")*cut);
  tree1->Project("hist4c_1", p4, TCut("Weight")*cut);
  tree1->Project("hist5c_1", p5, TCut("Weight")*cut);
  tree1->Project("hist6c_1", p6, TCut("Weight")*cut);
  tree1->Project("hist7c_1", p7, TCut("Weight")*cut);
  tree1->Project("hist8c_1", p8, TCut("Weight")*cut);
  tree1->Project("hist9c_1", p9, TCut("Weight")*cut);


  hist1c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist1c_1->GetXaxis()->SetTitle(hist1c_1->GetTitle());
  hist1c_1->SetTitle(""); 
  hist1c_1->GetXaxis()->SetTitleSize(0.05);
  hist1c_1->GetYaxis()->SetTitleSize(0.06);
  hist1c_1->GetYaxis()->SetTitleOffset(0.8);
  hist1c_1->GetXaxis()->SetTitleOffset(0.75);
  hist1c_1->GetXaxis()->SetLabelSize(0.04);
  hist1c_1->GetYaxis()->SetLabelSize(0.04);
  hist1c_1->GetXaxis()->SetNdivisions(6);
  hist1c_1->GetYaxis()->SetNdivisions(6);
  hist1c_1->GetYaxis()->SetMaxDigits(2);

  hist2c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist2c_1->GetXaxis()->SetTitle(hist2c_1->GetTitle());
  hist2c_1->SetTitle(""); 
  hist2c_1->GetXaxis()->SetTitleSize(0.05);
  hist2c_1->GetYaxis()->SetTitleSize(0.06);
  hist2c_1->GetYaxis()->SetTitleOffset(0.8);
  hist2c_1->GetXaxis()->SetTitleOffset(0.75);
  hist2c_1->GetXaxis()->SetLabelSize(0.04);
  hist2c_1->GetYaxis()->SetLabelSize(0.04);
  hist2c_1->GetXaxis()->SetNdivisions(6);
  hist2c_1->GetYaxis()->SetNdivisions(6);
  hist2c_1->GetYaxis()->SetMaxDigits(2);

  hist3c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist3c_1->GetXaxis()->SetTitle(hist3c_1->GetTitle());
  hist3c_1->SetTitle(""); 
  hist3c_1->GetXaxis()->SetTitleSize(0.05);
  hist3c_1->GetYaxis()->SetTitleSize(0.06);
  hist3c_1->GetYaxis()->SetTitleOffset(0.8);
  hist3c_1->GetXaxis()->SetTitleOffset(0.75);
  hist3c_1->GetXaxis()->SetLabelSize(0.04);
  hist3c_1->GetYaxis()->SetLabelSize(0.04);
  hist3c_1->GetXaxis()->SetNdivisions(6);
  hist3c_1->GetYaxis()->SetNdivisions(6);
  hist3c_1->GetYaxis()->SetMaxDigits(2);

  hist4c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist4c_1->GetXaxis()->SetTitle(hist4c_1->GetTitle());
  hist4c_1->SetTitle(""); 
  hist4c_1->GetXaxis()->SetTitleSize(0.05);
  hist4c_1->GetYaxis()->SetTitleSize(0.06);
  hist4c_1->GetYaxis()->SetTitleOffset(0.8);
  hist4c_1->GetXaxis()->SetTitleOffset(0.75);
  hist4c_1->GetXaxis()->SetLabelSize(0.04);
  hist4c_1->GetYaxis()->SetLabelSize(0.04);
  hist4c_1->GetXaxis()->SetNdivisions(6);
  hist4c_1->GetYaxis()->SetNdivisions(6);
  hist4c_1->GetYaxis()->SetMaxDigits(2);

  hist5c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist5c_1->GetXaxis()->SetTitle(hist5c_1->GetTitle());
  hist5c_1->SetTitle(""); 
  hist5c_1->GetXaxis()->SetTitleSize(0.05);
  hist5c_1->GetYaxis()->SetTitleSize(0.06);
  hist5c_1->GetYaxis()->SetTitleOffset(0.8);
  hist5c_1->GetXaxis()->SetTitleOffset(0.75);
  hist5c_1->GetXaxis()->SetLabelSize(0.04);
  hist5c_1->GetYaxis()->SetLabelSize(0.04);
  hist5c_1->GetXaxis()->SetNdivisions(6);
  hist5c_1->GetYaxis()->SetNdivisions(6);
  hist5c_1->GetYaxis()->SetMaxDigits(2);

  hist6c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist6c_1->GetXaxis()->SetTitle(hist6c_1->GetTitle());
  hist6c_1->SetTitle(""); 
  hist6c_1->GetXaxis()->SetTitleSize(0.05);
  hist6c_1->GetYaxis()->SetTitleSize(0.06);
  hist6c_1->GetYaxis()->SetTitleOffset(0.8);
  hist6c_1->GetXaxis()->SetTitleOffset(0.75);
  hist6c_1->GetXaxis()->SetLabelSize(0.04);
  hist6c_1->GetYaxis()->SetLabelSize(0.04);
  hist6c_1->GetXaxis()->SetNdivisions(6);
  hist6c_1->GetYaxis()->SetNdivisions(6);
  hist6c_1->GetYaxis()->SetMaxDigits(2);

  hist7c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist7c_1->GetXaxis()->SetTitle(hist7c_1->GetTitle());
  hist7c_1->SetTitle(""); 
  hist7c_1->GetXaxis()->SetTitleSize(0.05);
  hist7c_1->GetYaxis()->SetTitleSize(0.06);
  hist7c_1->GetYaxis()->SetTitleOffset(0.8);
  hist7c_1->GetXaxis()->SetTitleOffset(0.75);
  hist7c_1->GetXaxis()->SetLabelSize(0.04);
  hist7c_1->GetYaxis()->SetLabelSize(0.04);
  hist7c_1->GetXaxis()->SetNdivisions(6);
  hist7c_1->GetYaxis()->SetNdivisions(6);
  hist7c_1->GetYaxis()->SetMaxDigits(2);

  hist8c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist8c_1->GetXaxis()->SetTitle(hist8c_1->GetTitle());
  hist8c_1->SetTitle(""); 
  hist8c_1->GetXaxis()->SetTitleSize(0.05);
  hist8c_1->GetYaxis()->SetTitleSize(0.06);
  hist8c_1->GetYaxis()->SetTitleOffset(0.8);
  hist8c_1->GetXaxis()->SetTitleOffset(0.75);
  hist8c_1->GetXaxis()->SetLabelSize(0.04);
  hist8c_1->GetYaxis()->SetLabelSize(0.04);
  hist8c_1->GetXaxis()->SetNdivisions(6);
  hist8c_1->GetYaxis()->SetNdivisions(6);
  hist8c_1->GetYaxis()->SetMaxDigits(2);

  hist9c_1->GetYaxis()->SetTitle("Counts/total events"); 
  hist9c_1->GetXaxis()->SetTitle(hist9c_1->GetTitle());
  hist9c_1->SetTitle(""); 
  hist9c_1->GetXaxis()->SetTitleSize(0.05);
  hist9c_1->GetYaxis()->SetTitleSize(0.06);
  hist9c_1->GetYaxis()->SetTitleOffset(0.8);
  hist9c_1->GetXaxis()->SetTitleOffset(0.75);
  hist9c_1->GetXaxis()->SetLabelSize(0.04);
  hist9c_1->GetYaxis()->SetLabelSize(0.04);
  hist9c_1->GetXaxis()->SetNdivisions(6);
  hist9c_1->GetYaxis()->SetNdivisions(6);
  hist9c_1->GetYaxis()->SetMaxDigits(2);


  hist1c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist1c_2->GetXaxis()->SetTitle(hist1c_2->GetTitle());
  hist1c_2->SetTitle(""); 
  hist1c_2->GetXaxis()->SetTitleSize(0.05);
  hist1c_2->GetYaxis()->SetTitleSize(0.06);
  hist1c_2->GetYaxis()->SetTitleOffset(0.8);
  hist1c_2->GetXaxis()->SetTitleOffset(0.75);
  hist1c_2->GetXaxis()->SetLabelSize(0.04);
  hist1c_2->GetYaxis()->SetLabelSize(0.04);
  hist1c_2->GetXaxis()->SetNdivisions(6);
  hist1c_2->GetYaxis()->SetNdivisions(6);
  hist1c_2->GetYaxis()->SetMaxDigits(2);

  hist2c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist2c_2->GetXaxis()->SetTitle(hist2c_2->GetTitle());
  hist2c_2->SetTitle(""); 
  hist2c_2->GetXaxis()->SetTitleSize(0.05);
  hist2c_2->GetYaxis()->SetTitleSize(0.06);
  hist2c_2->GetYaxis()->SetTitleOffset(0.8);
  hist2c_2->GetXaxis()->SetTitleOffset(0.75);
  hist2c_2->GetXaxis()->SetLabelSize(0.04);
  hist2c_2->GetYaxis()->SetLabelSize(0.04);
  hist2c_2->GetXaxis()->SetNdivisions(6);
  hist2c_2->GetYaxis()->SetNdivisions(6);
  hist2c_2->GetYaxis()->SetMaxDigits(2);

  hist3c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist3c_2->GetXaxis()->SetTitle(hist3c_2->GetTitle());
  hist3c_2->SetTitle(""); 
  hist3c_2->GetXaxis()->SetTitleSize(0.05);
  hist3c_2->GetYaxis()->SetTitleSize(0.06);
  hist3c_2->GetYaxis()->SetTitleOffset(0.8);
  hist3c_2->GetXaxis()->SetTitleOffset(0.75);
  hist3c_2->GetXaxis()->SetLabelSize(0.04);
  hist3c_2->GetYaxis()->SetLabelSize(0.04);
  hist3c_2->GetXaxis()->SetNdivisions(6);
  hist3c_2->GetYaxis()->SetNdivisions(6);
  hist3c_2->GetYaxis()->SetMaxDigits(2);

  hist4c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist4c_2->GetXaxis()->SetTitle(hist4c_2->GetTitle());
  hist4c_2->SetTitle(""); 
  hist4c_2->GetXaxis()->SetTitleSize(0.05);
  hist4c_2->GetYaxis()->SetTitleSize(0.06);
  hist4c_2->GetYaxis()->SetTitleOffset(0.8);
  hist4c_2->GetXaxis()->SetTitleOffset(0.75);
  hist4c_2->GetXaxis()->SetLabelSize(0.04);
  hist4c_2->GetYaxis()->SetLabelSize(0.04);
  hist4c_2->GetXaxis()->SetNdivisions(6);
  hist4c_2->GetYaxis()->SetNdivisions(6);
  hist4c_2->GetYaxis()->SetMaxDigits(2);

  hist5c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist5c_2->GetXaxis()->SetTitle(hist5c_2->GetTitle());
  hist5c_2->SetTitle(""); 
  hist5c_2->GetXaxis()->SetTitleSize(0.05);
  hist5c_2->GetYaxis()->SetTitleSize(0.06);
  hist5c_2->GetYaxis()->SetTitleOffset(0.8);
  hist5c_2->GetXaxis()->SetTitleOffset(0.75);
  hist5c_2->GetXaxis()->SetLabelSize(0.04);
  hist5c_2->GetYaxis()->SetLabelSize(0.04);
  hist5c_2->GetXaxis()->SetNdivisions(6);
  hist5c_2->GetYaxis()->SetNdivisions(6);
  hist5c_2->GetYaxis()->SetMaxDigits(2);

  hist6c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist6c_2->GetXaxis()->SetTitle(hist6c_2->GetTitle());
  hist6c_2->SetTitle(""); 
  hist6c_2->GetXaxis()->SetTitleSize(0.05);
  hist6c_2->GetYaxis()->SetTitleSize(0.06);
  hist6c_2->GetYaxis()->SetTitleOffset(0.8);
  hist6c_2->GetXaxis()->SetTitleOffset(0.75);
  hist6c_2->GetXaxis()->SetLabelSize(0.04);
  hist6c_2->GetYaxis()->SetLabelSize(0.04);
  hist6c_2->GetXaxis()->SetNdivisions(6);
  hist6c_2->GetYaxis()->SetNdivisions(6);
  hist6c_2->GetYaxis()->SetMaxDigits(2);

  hist7c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist7c_2->GetXaxis()->SetTitle(hist7c_2->GetTitle());
  hist7c_2->SetTitle(""); 
  hist7c_2->GetXaxis()->SetTitleSize(0.05);
  hist7c_2->GetYaxis()->SetTitleSize(0.06);
  hist7c_2->GetYaxis()->SetTitleOffset(0.8);
  hist7c_2->GetXaxis()->SetTitleOffset(0.75);
  hist7c_2->GetXaxis()->SetLabelSize(0.04);
  hist7c_2->GetYaxis()->SetLabelSize(0.04);
  hist7c_2->GetXaxis()->SetNdivisions(6);
  hist7c_2->GetYaxis()->SetNdivisions(6);
  hist7c_2->GetYaxis()->SetMaxDigits(2);

  hist8c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist8c_2->GetXaxis()->SetTitle(hist8c_2->GetTitle());
  hist8c_2->SetTitle(""); 
  hist8c_2->GetXaxis()->SetTitleSize(0.05);
  hist8c_2->GetYaxis()->SetTitleSize(0.06);
  hist8c_2->GetYaxis()->SetTitleOffset(0.8);
  hist8c_2->GetXaxis()->SetTitleOffset(0.75);
  hist8c_2->GetXaxis()->SetLabelSize(0.04);
  hist8c_2->GetYaxis()->SetLabelSize(0.04);
  hist8c_2->GetXaxis()->SetNdivisions(6);
  hist8c_2->GetYaxis()->SetNdivisions(6);
  hist8c_2->GetYaxis()->SetMaxDigits(2);

  hist9c_2->GetYaxis()->SetTitle("Counts/total events"); 
  hist9c_2->GetXaxis()->SetTitle(hist9c_2->GetTitle());
  hist9c_2->SetTitle(""); 
  hist9c_2->GetXaxis()->SetTitleSize(0.05);
  hist9c_2->GetYaxis()->SetTitleSize(0.06);
  hist9c_2->GetYaxis()->SetTitleOffset(0.8);
  hist9c_2->GetXaxis()->SetTitleOffset(0.75);
  hist9c_2->GetXaxis()->SetLabelSize(0.04);
  hist9c_2->GetYaxis()->SetLabelSize(0.04);
  hist9c_2->GetXaxis()->SetNdivisions(6);
  hist9c_2->GetYaxis()->SetNdivisions(6);
  hist9c_2->GetYaxis()->SetMaxDigits(2);


  hist1c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist1c_3->GetXaxis()->SetTitle(hist1c_3->GetTitle());
  hist1c_3->SetTitle(""); 
  hist1c_3->GetXaxis()->SetTitleSize(0.05);
  hist1c_3->GetYaxis()->SetTitleSize(0.06);
  hist1c_3->GetYaxis()->SetTitleOffset(0.8);
  hist1c_3->GetXaxis()->SetTitleOffset(0.75);
  hist1c_3->GetXaxis()->SetLabelSize(0.04);
  hist1c_3->GetYaxis()->SetLabelSize(0.04);
  hist1c_3->GetXaxis()->SetNdivisions(6);
  hist1c_3->GetYaxis()->SetNdivisions(6);
  hist1c_3->GetYaxis()->SetMaxDigits(2);

  hist2c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist2c_3->GetXaxis()->SetTitle(hist2c_3->GetTitle());
  hist2c_3->SetTitle(""); 
  hist2c_3->GetXaxis()->SetTitleSize(0.05);
  hist2c_3->GetYaxis()->SetTitleSize(0.06);
  hist2c_3->GetYaxis()->SetTitleOffset(0.8);
  hist2c_3->GetXaxis()->SetTitleOffset(0.75);
  hist2c_3->GetXaxis()->SetLabelSize(0.04);
  hist2c_3->GetYaxis()->SetLabelSize(0.04);
  hist2c_3->GetXaxis()->SetNdivisions(6);
  hist2c_3->GetYaxis()->SetNdivisions(6);
  hist2c_3->GetYaxis()->SetMaxDigits(2);

  hist3c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist3c_3->GetXaxis()->SetTitle(hist3c_3->GetTitle());
  hist3c_3->SetTitle(""); 
  hist3c_3->GetXaxis()->SetTitleSize(0.05);
  hist3c_3->GetYaxis()->SetTitleSize(0.06);
  hist3c_3->GetYaxis()->SetTitleOffset(0.8);
  hist3c_3->GetXaxis()->SetTitleOffset(0.75);
  hist3c_3->GetXaxis()->SetLabelSize(0.04);
  hist3c_3->GetYaxis()->SetLabelSize(0.04);
  hist3c_3->GetXaxis()->SetNdivisions(6);
  hist3c_3->GetYaxis()->SetNdivisions(6);
  hist3c_3->GetYaxis()->SetMaxDigits(2);

  hist4c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist4c_3->GetXaxis()->SetTitle(hist4c_3->GetTitle());
  hist4c_3->SetTitle(""); 
  hist4c_3->GetXaxis()->SetTitleSize(0.05);
  hist4c_3->GetYaxis()->SetTitleSize(0.06);
  hist4c_3->GetYaxis()->SetTitleOffset(0.8);
  hist4c_3->GetXaxis()->SetTitleOffset(0.75);
  hist4c_3->GetXaxis()->SetLabelSize(0.04);
  hist4c_3->GetYaxis()->SetLabelSize(0.04);
  hist4c_3->GetXaxis()->SetNdivisions(6);
  hist4c_3->GetYaxis()->SetNdivisions(6);
  hist4c_3->GetYaxis()->SetMaxDigits(2);

  hist5c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist5c_3->GetXaxis()->SetTitle(hist5c_3->GetTitle());
  hist5c_3->SetTitle(""); 
  hist5c_3->GetXaxis()->SetTitleSize(0.05);
  hist5c_3->GetYaxis()->SetTitleSize(0.06);
  hist5c_3->GetYaxis()->SetTitleOffset(0.8);
  hist5c_3->GetXaxis()->SetTitleOffset(0.75);
  hist5c_3->GetXaxis()->SetLabelSize(0.04);
  hist5c_3->GetYaxis()->SetLabelSize(0.04);
  hist5c_3->GetXaxis()->SetNdivisions(6);
  hist5c_3->GetYaxis()->SetNdivisions(6);
  hist5c_3->GetYaxis()->SetMaxDigits(2);

  hist6c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist6c_3->GetXaxis()->SetTitle(hist6c_3->GetTitle());
  hist6c_3->SetTitle(""); 
  hist6c_3->GetXaxis()->SetTitleSize(0.05);
  hist6c_3->GetYaxis()->SetTitleSize(0.06);
  hist6c_3->GetYaxis()->SetTitleOffset(0.8);
  hist6c_3->GetXaxis()->SetTitleOffset(0.75);
  hist6c_3->GetXaxis()->SetLabelSize(0.04);
  hist6c_3->GetYaxis()->SetLabelSize(0.04);
  hist6c_3->GetXaxis()->SetNdivisions(6);
  hist6c_3->GetYaxis()->SetNdivisions(6);
  hist6c_3->GetYaxis()->SetMaxDigits(2);

  hist7c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist7c_3->GetXaxis()->SetTitle(hist7c_3->GetTitle());
  hist7c_3->SetTitle(""); 
  hist7c_3->GetXaxis()->SetTitleSize(0.05);
  hist7c_3->GetYaxis()->SetTitleSize(0.06);
  hist7c_3->GetYaxis()->SetTitleOffset(0.8);
  hist7c_3->GetXaxis()->SetTitleOffset(0.75);
  hist7c_3->GetXaxis()->SetLabelSize(0.04);
  hist7c_3->GetYaxis()->SetLabelSize(0.04);
  hist7c_3->GetXaxis()->SetNdivisions(6);
  hist7c_3->GetYaxis()->SetNdivisions(6);
  hist7c_3->GetYaxis()->SetMaxDigits(2);

  hist8c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist8c_3->GetXaxis()->SetTitle(hist8c_3->GetTitle());
  hist8c_3->SetTitle(""); 
  hist8c_3->GetXaxis()->SetTitleSize(0.05);
  hist8c_3->GetYaxis()->SetTitleSize(0.06);
  hist8c_3->GetYaxis()->SetTitleOffset(0.8);
  hist8c_3->GetXaxis()->SetTitleOffset(0.75);
  hist8c_3->GetXaxis()->SetLabelSize(0.04);
  hist8c_3->GetYaxis()->SetLabelSize(0.04);
  hist8c_3->GetXaxis()->SetNdivisions(6);
  hist8c_3->GetYaxis()->SetNdivisions(6);
  hist8c_3->GetYaxis()->SetMaxDigits(2);

  hist9c_3->GetYaxis()->SetTitle("Counts/total events"); 
  hist9c_3->GetXaxis()->SetTitle(hist9c_3->GetTitle());
  hist9c_3->SetTitle(""); 
  hist9c_3->GetXaxis()->SetTitleSize(0.05);
  hist9c_3->GetYaxis()->SetTitleSize(0.06);
  hist9c_3->GetYaxis()->SetTitleOffset(0.8);
  hist9c_3->GetXaxis()->SetTitleOffset(0.75);
  hist9c_3->GetXaxis()->SetLabelSize(0.04);
  hist9c_3->GetYaxis()->SetLabelSize(0.04);
  hist9c_3->GetXaxis()->SetNdivisions(6);
  hist9c_3->GetYaxis()->SetNdivisions(6);
  hist9c_3->GetYaxis()->SetMaxDigits(2);


  c1->cd(1);
  hist1c_3->Scale(factor/hist1c_3->Integral());
  hist1c_3->SetLineColor(kBlue);
  hist1c_2->Scale(factor/hist1c_2->Integral());
  hist1c_2->SetLineColor(kBlack);
  hist1c_1->Scale(factor/hist1c_1->Integral());
  hist1c_1->SetLineColor(kRed);

  hist1c_3->Draw("Hist,Same");
  hist1c_2->Draw("Hist,Same");
  hist1c_1->Draw("Hist,Same");
  
  c1->cd(2);
  hist2c_3->Scale(factor/hist2c_3->Integral());
  hist2c_3->SetLineColor(kBlue);
  hist2c_2->Scale(factor/hist2c_2->Integral());
  hist2c_2->SetLineColor(kBlack);
  hist2c_1->Scale(factor/hist2c_1->Integral());
  hist2c_1->SetLineColor(kRed);

  hist2c_3->Draw("Hist,Same");
  hist2c_2->Draw("Hist,Same");
  hist2c_1->Draw("Hist,Same");


  c1->cd(3);
  hist3c_3->Scale(factor/hist3c_3->Integral());
  hist3c_3->SetLineColor(kBlue);
  hist3c_2->Scale(factor/hist3c_2->Integral());
  hist3c_2->SetLineColor(kBlack);
  hist3c_1->Scale(factor/hist3c_1->Integral());
  hist3c_1->SetLineColor(kRed);

  hist3c_3->Draw("Hist,Same");
  hist3c_2->Draw("Hist,Same");
  hist3c_1->Draw("Hist,Same");

  c1->cd(4);
  hist4c_3->Scale(factor/hist4c_3->Integral());
  hist4c_3->SetLineColor(kBlue);
  hist4c_2->Scale(factor/hist4c_2->Integral());
  hist4c_2->SetLineColor(kBlack);
  hist4c_1->Scale(factor/hist4c_1->Integral());
  hist4c_1->SetLineColor(kRed);

  hist4c_3->Draw("Hist,Same");
  hist4c_2->Draw("Hist,Same");
  hist4c_1->Draw("Hist,Same");


  c1->cd(5);
  hist5c_3->Scale(factor/hist5c_3->Integral());
  hist5c_3->SetLineColor(kBlue);
  hist5c_1->Scale(factor/hist5c_1->Integral());
  hist5c_1->SetLineColor(kRed);
  hist5c_2->Scale(factor/hist5c_2->Integral());
  hist5c_2->SetLineColor(kBlack);

  hist5c_3->Draw("Hist,Same");
  hist5c_2->Draw("Hist,Same");
  hist5c_1->Draw("Hist,Same");


  c1->cd(6);
  hist6c_3->Scale(factor/hist6c_3->Integral());
  hist6c_3->SetLineColor(kBlue);
  hist6c_2->Scale(factor/hist6c_2->Integral());
  hist6c_2->SetLineColor(kBlack);
  hist6c_1->Scale(factor/hist6c_1->Integral());
  hist6c_1->SetLineColor(kRed);

  hist6c_3->Draw("Hist,Same");
  hist6c_2->Draw("Hist,Same");
  hist6c_1->Draw("Hist,Same");

  
  c1->cd(7);
  hist7c_3->Scale(factor/hist7c_3->Integral());
  hist7c_3->SetLineColor(kBlue);
  hist7c_1->Scale(factor/hist7c_1->Integral());
  hist7c_1->SetLineColor(kRed);
  hist7c_2->Scale(factor/hist7c_2->Integral());
  hist7c_2->SetLineColor(kBlack);

  hist7c_3->Draw("Hist,Same");
  hist7c_2->Draw("Hist,Same");
  hist7c_1->Draw("Hist,Same");

  
  c1->cd(8);
  hist8c_3->Scale(factor/hist8c_3->Integral());
  hist8c_3->SetLineColor(kBlue);
  hist8c_2->Scale(factor/hist8c_2->Integral());
  hist8c_2->SetLineColor(kBlack);
  hist8c_1->Scale(factor/hist8c_1->Integral());
  hist8c_1->SetLineColor(kRed);

  hist8c_3->Draw("Hist,Same");
  hist8c_2->Draw("Hist,Same");
  hist8c_1->Draw("Hist,Same");

  
  c1->cd(9);
  hist9c_3->Scale(factor/hist9c_3->Integral());
  hist9c_3->SetLineColor(kBlue);
  hist9c_2->Scale(factor/hist9c_2->Integral());
  hist9c_2->SetLineColor(kBlack);
  hist9c_1->Scale(factor/hist9c_1->Integral());
  hist9c_1->SetLineColor(kRed);

  hist9c_3->Draw("Hist,Same");
  hist9c_2->Draw("Hist,Same");
  hist9c_1->Draw("Hist,Same");

  c1->Print(Folder + TString("three_comparison.pdf"));

  delete c1;
  /*
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

  delete p1;
  delete p2;
  delete p3;
  delete p4;
  delete p5;
  delete p6;
  delete p7;
  delete p8;
  delete p9;
  */
  delete tree1;
  delete tree2;
  delete tree3;
  input2->Close();
  input3->Close();
  input4->Close();
  delete input2;
  delete input3;
  delete input4;

}

