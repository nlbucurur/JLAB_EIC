
void BDT::Single_BSA_2(TString BDT_Data, TString Bin_Data, vector<double> boundaries, int P_bins=12)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Folder + Bin_Data);

  TChain *ch2= new TChain("pDVCS");
  ch2->Add(Folder + BDT_Data);

  double meanQ;
  double meant;
  double meanx;

  TF1 *fit;
  
  TCanvas* c1 = new TCanvas("c1","Histograms");
		  
  TCut cut_p = TCut("(bestCandidateFlag==1 && Helicity>0)");
  TCut cut_m = TCut("(bestCandidateFlag==1 && Helicity<0)");

  //--------------------------------------
  //Get Means

  TCut cutM = TCut("bestCandidateFlag==1");
	  
  TH1F *tmean = new TH1F("tmean","t_temp",100,0,0);
  ch2->Project("tmean","t_Ph",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meant=tmean->GetMean();

  TH1F *Qmean = new TH1F("Qmean","Q_temp",100,0,0);
  ch2->Project("Qmean","strip_Q2",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanQ=Qmean->GetMean();

  TH1F *xmean = new TH1F("xmean","x_temp",100,0,0);
  ch2->Project("xmean","strip_Xbj",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanx=xmean->GetMean();

  TH1F *histp = new TH1F("histp","1dgraph",P_bins,0,360);
  TH1F *histm = new TH1F("histm","1dgraph",P_bins,0,360);

  TH1F *histp_F = new TH1F("histp_F","1dgraph",P_bins,0,360);
  TH1F *histm_F = new TH1F("histm_F","1dgraph",P_bins,0,360);

  ch1->Project("histp","Phi_Nuc",cut_p);
  ch1->Project("histm","Phi_Nuc",cut_m);
	  
  TH1 * BA=histm->GetAsymmetry(histp);
  BA->Scale(1.0/Bpol);
  //BA->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  delete tmean;
  delete Qmean;
  delete xmean;
  //----------------------------------------

  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf = new TF1("fitf","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf->SetParameter(0,0.1);
  fitf->SetParameter(1,-0.3);
  fitf->SetParLimits(0,0,0.3);
  fitf->SetParLimits(1,-1.,1.);
  fitf->SetLineColor(kWhite);
  //fitf->SetLineWidth(2);

	  
  //BA->Fit("fitf");
  BA->SetLineColor(kBlack);
  BA->SetLineWidth(2);
  BA->SetMarkerColor(kBlack);
  BA->SetAxisRange(-1.0, 1.0,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
	  
  //Display means
  TLatex *box = new TLatex(.6,.8,Form("#splitline{#LT t #GT= %.3f}{#splitline{#LT Q^{2} #GT= %.3f}{#LT x_{B} #GT= %.3f}}",meant,meanQ,meanx));
  box->SetNDC(kTRUE);

  ch2->Project("histp_F","Phi_Nuc",cut_p + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  ch2->Project("histm_F","Phi_Nuc",cut_m + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
	  
  TH1 * BA2=histm_F->GetAsymmetry(histp_F);
  BA2->Scale(1.0/Bpol);
  //BA2->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf2 = new TF1("fitf2","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf2->SetParameter(0,0.1);
  fitf2->SetParameter(1,-0.3);
  fitf2->SetParLimits(0,0,0.3);
  fitf2->SetParLimits(1,-1.,1.);
  fitf2->SetLineColor(kWhite);
	  
  //BA2->Fit("fitf2");
  BA2->SetLineColor(kRed);
  BA2->SetLineWidth(2);
  BA2->SetMarkerColor(kRed);
  BA2->SetAxisRange(-1.0, 1.0,"Y");
  BA2->SetAxisRange( 0. ,360.,"X");

  THStack *hist = new THStack("hist","");
  //hist->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  BA->SetTitle("Binned");
  BA2->SetTitle("Unbinned");

  hist->Add(BA);
  hist->Add(BA2);
  hist->Draw("nostack");
  c1->BuildLegend();	  
 
  hist->GetXaxis()->SetTitle("#phi(deg)");
  hist->GetYaxis()->SetTitle("BSA");
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.5);
  hist->GetXaxis()->SetTitleOffset(0.5);
  hist->GetYaxis()->SetNdivisions(6);
  hist->GetXaxis()->SetNdivisions(4);

  //box->DrawClone("SAME");
  c1->Print(Folder + "BSA_Global_vs_Bin.png");  
  delete histp;
  delete histm;
  delete BA;
  delete ch1;
  delete histp_F;
  delete histm_F;
  delete BA2;
  delete ch2;
  delete c1;
  return;
}
