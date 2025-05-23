
void BDT::Get_Contamination(TCut cut, double BDT_cut)
{
  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Folder + TDVCS);

  TFile *input4 = new TFile(Folder + TPi0,"READ");
  TTree *ch2;
  if( input4->GetListOfKeys()->Contains("pDVCS") )
    ch2     = (TTree*)input4->Get("pDVCS");
  else
    ch2     = (TTree*)input4->Get("eppi0");

  TChain *ch3= new TChain("pDVCS");
  ch3->Add(Folder + TData);

  double N_DVCS_Train, N_Pi0_Train, N_before, N_before_FT, N_before_FD;
  double N_DVCS_after_BDT, N_Pi0_after_BDT, N_after;
  double N0_hDVCS, N0_hPi0;
  double sig_eff, bkg_eff;

  TH1F *hDVCS_0 = new TH1F("hDVCS_0","hDVCS_0",100,0,4);
  ch1->Project("hDVCS_0","mm2_eg",cut);
  N_DVCS_Train=hDVCS_0->GetEntries();

  TH1F *hPi0_0 = new TH1F("hPi0_0","hPi0_0",100,0,4);
  ch2->Project("hPi0_0","mm2_eg",cut);
  N_Pi0_Train=hPi0_0->GetEntries();
  
  TH1F *hData_0 = new TH1F("hData_0","hData_0",100,0,4);
  ch3->Project("hData_0","mm2_eg",cut);
  N_before=hData_0->GetEntries();

  TH1F *hData_0_FT = new TH1F("hData_0_FT","hData_0_FT",100,0,4);
  ch3->Project("hData_0_FT","mm2_eg",cut + TCut("strip_Ph_Theta < 5"));
  N_before_FT=hData_0_FT->GetEntries();

  TH1F *hData_0_FD = new TH1F("hData_0_FD","hData_0",100,0,4);
  ch3->Project("hData_0_FD","mm2_eg",cut + TCut("strip_Ph_Theta > 5"));
  N_before_FD=hData_0_FD->GetEntries();
  std::cout<<"Number of events on bin (All/FT/FD) "<<N_before<<" "<<N_before_FT<<" "<<N_before_FD<<endl;
  std::cout<<"*******************************************"<<endl;
  std::cout<<"*******************************************"<<endl;
  entries_bef_BDT=N_before;
  entries_bef_BDT_FT=N_before_FT;
  entries_bef_BDT_FD=N_before_FD;
  
  if(N_DVCS_Train==0 || N_Pi0_Train==0 || N_before==0)
    {
      std::cout<<"ERROR: Unable to get contamination from BDT"<<endl;
      std::cout<<N_DVCS_Train <<" "<< N_Pi0_Train <<" "<< N_before<<endl;
      delete ch1;
      delete ch2;
      delete ch3;
      
      delete hDVCS_0;
      delete hPi0_0;
      delete hData_0;
      boundaries.push_back(-1);
      boundaries.push_back(-1);

      return;
    }
  TH1F *hDVCS = new TH1F("hDVCS","hDVCS",100,0,4);
  ch1->Project("hDVCS","mm2_eg",cut && TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
  N_DVCS_after_BDT=hDVCS->GetEntries();

  TH1F *hPi0 = new TH1F("hPi0","hPi0",100,0,4);
  ch2->Project("hPi0","mm2_eg",cut && TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
  N_Pi0_after_BDT=hPi0->GetEntries();
  
  TH1F *hData = new TH1F("hData","hData",100,0,4);
  ch3->Project("hData","mm2_eg",cut && TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
  N_after=hData->GetEntries();

  sig_eff=N_DVCS_after_BDT/N_DVCS_Train;
  bkg_eff=N_Pi0_after_BDT/N_Pi0_Train;

  N0_hDVCS=-(-N_after + bkg_eff*N_before)/(sig_eff - bkg_eff);
  N0_hPi0 = (-N_after + sig_eff*N_before)/(sig_eff - bkg_eff);

  std::cout<<"BDT Estimation "<<N0_hPi0/N_before<<" "<<N0_hPi0*bkg_eff/N_after<<endl;
  boundaries.push_back(N0_hPi0/N_before);
  boundaries.push_back(N0_hPi0*bkg_eff/N_after);

  
  TCanvas* c1 = new TCanvas("c1","Histograms",1000,500);
  c1->Divide(2,1);

  hDVCS_0->Scale(N0_hDVCS/N_DVCS_Train);
  hPi0_0->Scale(N0_hPi0/N_Pi0_Train);
  hDVCS_0->Add(hPi0_0);
  TH1F *Mixed_MC_before=hDVCS_0;

  hDVCS->Scale(N0_hDVCS*sig_eff/N_DVCS_after_BDT);
  hPi0->Scale(N0_hPi0*bkg_eff/N_Pi0_after_BDT);
  hDVCS->Add(hPi0);
  TH1F *Mixed_MC_after=hDVCS;
  
  c1->cd(1);
  hData_0->SetLineColor(kRed);
  hData_0->SetTitle("Before");
  Mixed_MC_before->SetTitle("Before");
  Mixed_MC_before->GetXaxis()->SetTitle("M_{e#gamma X}^{2} GeV^{2}");
  Mixed_MC_before->GetYaxis()->SetTitle("counts");
  Mixed_MC_before->GetXaxis()->SetTitleSize(0.05);
  Mixed_MC_before->GetYaxis()->SetTitleSize(0.05);
  Mixed_MC_before->GetXaxis()->SetLabelSize(0.04);
  Mixed_MC_before->GetYaxis()->SetLabelSize(0.04);
  Mixed_MC_before->GetYaxis()->SetTitleOffset(0.9);
  Mixed_MC_before->GetXaxis()->SetTitleOffset(0.7);
  Mixed_MC_before->GetYaxis()->SetNdivisions(6);
  Mixed_MC_before->GetXaxis()->SetNdivisions(4);
  Mixed_MC_before->GetYaxis()->SetMaxDigits(2);
  Mixed_MC_before->Draw("HIST");
  hData_0->Draw("HIST,SAME");

  c1->cd(2);
  hData->SetLineColor(kRed);
  hData->SetTitle("After");
  Mixed_MC_after->SetTitle("After");
  Mixed_MC_after->GetXaxis()->SetTitle("M_{e#gamma X}^{2} GeV^{2}");
  Mixed_MC_after->GetYaxis()->SetTitle("counts");
  Mixed_MC_after->GetXaxis()->SetTitleSize(0.05);
  Mixed_MC_after->GetYaxis()->SetTitleSize(0.05);
  Mixed_MC_after->GetXaxis()->SetLabelSize(0.04);
  Mixed_MC_after->GetYaxis()->SetLabelSize(0.04);
  Mixed_MC_after->GetYaxis()->SetTitleOffset(0.9);
  Mixed_MC_after->GetXaxis()->SetTitleOffset(0.7);
  Mixed_MC_after->GetYaxis()->SetNdivisions(6);
  Mixed_MC_after->GetXaxis()->SetNdivisions(4);
  Mixed_MC_after->GetYaxis()->SetMaxDigits(2);
  Mixed_MC_after->Draw("HIST");
  hData->Draw("HIST,SAME");

  c1->Print(Folder + TString("Contamination.pdf"));
    
  delete hData;
  delete Mixed_MC_before;
  delete Mixed_MC_after;
  delete ch1;
  delete ch2;
  delete ch3;
  delete c1;
  
  delete hPi0_0;
  delete hData_0;
  delete hPi0;

  input4->Close();
  delete input4;
  return;
}
