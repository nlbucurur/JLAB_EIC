TH1F* BDT::No_Rec_correction(TString Data, TString MCData, double BDT_cut, int Nphi){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString String_cut = cut_bin.GetTitle();
  String_cut.ReplaceAll("strip_Q2", "strip_Q2_MC");
  String_cut.ReplaceAll("strip_Xbj", "strip_Xbj_MC");
  String_cut.ReplaceAll("t_Ph", "t_Ph_MC");
  TCut cut_bin_gen = TCut(String_cut);

  TChain *Data_Tree= new TChain("pDVCS");
  Data_Tree->Add(Data);
  
  TChain *MCData_Tree= new TChain("pDVCS");
  MCData_Tree->Add(MCData);

  TH1F *Phi_rec = new TH1F("Phi_rec","",Nphi,0,360);
  TH1F *Phi_mc = new TH1F("Phi_mc","",Nphi,0,360);
  TH1F *ratio = new TH1F("Phi_norec","F_noRec",Nphi,0,360);

  Data_Tree  ->Project("Phi_rec", "Phi_Ph", (cut + cut_bin_gen + cut_ref + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut))));
  MCData_Tree->Project("Phi_mc" , "Phi_Ph", TCut("bestCandidateFlag==1") + cut_bin);

  ratio->Divide(Phi_rec, Phi_mc,1.0,1.0);
  
  std::cout<<"Acc correction"<<endl;
  for(int k=1; k<=Nphi; k++)
  {
    //std::cout<<Phi_rec->GetBinContent(k)<<" "<<Phi_mc->GetBinContent(k)<<" "<<ratio->GetBinContent(k)<<std::endl;
    if(ratio->GetBinContent(k)==0)
      ratio->SetBinError(k,0);
    else
      ratio->SetBinError(k,ratio->GetBinContent(k)*sqrt(1.0/Phi_rec->GetBinContent(k) + 1.0/Phi_mc->GetBinContent(k)));
    
      std::cout<<ratio->GetBinContent(k)<<" "<<ratio->GetBinError(k)<<std::endl;
  }

  delete Phi_rec;
  delete Phi_mc;
  delete Data_Tree;
  delete MCData_Tree;

  return ratio;
}
