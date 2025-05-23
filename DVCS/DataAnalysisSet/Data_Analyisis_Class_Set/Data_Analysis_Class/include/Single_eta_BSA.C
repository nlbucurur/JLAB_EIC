
TH1* BDT::Single_eta_BSA(TString Data, int P_bins, TCut cut)
{
  //Convert cut to string
  TString String_cut = cut.GetTitle();

  //replace Ph by eta
  String_cut.ReplaceAll("Ph ", "eta ");
  String_cut.ReplaceAll("Ph>", "eta>");
  String_cut.ReplaceAll("Ph<", "eta<");
  String_cut.ReplaceAll("strip_Ph_P", "strip_W");
  String_cut.ReplaceAll("gamma", "eta");

  //Create cut for 2gamma case
  TCut cut2g = TCut(String_cut) + TCut("mm2_egg<1.5 && strip_eta_2DChi2 < 0.04 ");
  
  
  TChain *ch1= new TChain("epeta");
  ch1->Add(Data);

  TCanvas* c1 = new TCanvas("c1","Histograms");
		  
  TCut cut_p = TCut("bestCandidateFlag==1 && Helicity>0");
  TCut cut_m = TCut("bestCandidateFlag==1 && Helicity<0");

  TH1F *histp = new TH1F("histp","1dgraph",P_bins,0,360);
  TH1F *histm = new TH1F("histm","1dgraph",P_bins,0,360);
  histp->Sumw2();
  histm->Sumw2();
	  
  ch1->Project("histp","Phi_Nuc",cut2g + cut_p);
  ch1->Project("histm","Phi_Nuc",cut2g + cut_m);
	  
  TH1 * BA=histm->GetAsymmetry(histp);
  BA->SetLineColor(kBlack);
  BA->SetMarkerColor(kBlack);
  BA->SetAxisRange(-1.0, 1.0,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi(deg)");

  BA->SetTitle("");
  BA->GetXaxis()->SetTitleSize(0.06);
  BA->GetYaxis()->SetTitle("BSA");
  BA->GetXaxis()->SetTitleSize(0.07);
  BA->GetYaxis()->SetTitleSize(0.07);
  BA->GetXaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetTitleOffset(0.3);
  BA->GetXaxis()->SetTitleOffset(0.3);
  BA->GetYaxis()->SetNdivisions(4);
  BA->GetXaxis()->SetNdivisions(4);

  BA->DrawClone();
	  
  BSA_Amplitude=BA->GetBinContent(BA->GetMaximumBin());
  std::cout<<"Amplitude: "<<BSA_Amplitude<<endl;
	  
  c1->Print(Folder + "BSA_eta.png");  
  delete histp;
  delete histm;
  delete ch1;
  delete c1;

  return BA;
}
