
void BDT::Single_BSA_Fit(TString BDT_Data, TString Fit_Data, vector<double> boundaries, int P_bins=12)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Folder + Fit_Data);

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
  ch2->Project("tmean","t_Ph",cutM);
  meant=tmean->GetMean();

  TH1F *Qmean = new TH1F("Qmean","Q_temp",100,0,0);
  ch2->Project("Qmean","strip_Q2",cutM);
  meanQ=Qmean->GetMean();

  TH1F *xmean = new TH1F("xmean","x_temp",100,0,0);
  ch2->Project("xmean","strip_Xbj",cutM);
  meanx=xmean->GetMean();

  TH1F *histp = new TH1F("histp","1dgraph",P_bins,0,360);
  TH1F *histm = new TH1F("histm","1dgraph",P_bins,0,360);

  TH1F *histp_F = new TH1F("histp_F","1dgraph",P_bins,0,360);
  TH1F *histm_F = new TH1F("histm_F","1dgraph",P_bins,0,360);

  ch1->Project("histp","Phi_Nuc",TCut(cut_p.GetTitle() + TString("*(SyieldWeight)")));
  ch1->Project("histm","Phi_Nuc",TCut(cut_m.GetTitle() + TString("*(SyieldWeight)")));
	  
  TH1 * BA=histm->GetAsymmetry(histp);
  BA->Scale(1.0/Bpol);
  BA->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

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
  fitf->SetLineColor(kBlack);

	  
  BA->Fit("fitf");
  BA->SetLineColor(kBlack);
  BA->SetMarkerColor(kBlack);
  BA->SetAxisRange(-0.3, 0.3,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi (deg)");
  //BA->DrawClone("SAME");
	  
  //Display means
  TLatex *box = new TLatex(.6,.8,Form("#splitline{#LT t #GT= %.3f}{#splitline{#LT Q^{2} #GT= %.3f}{#LT x_{B} #GT= %.3f}}",meant,meanQ,meanx));
  box->SetNDC(kTRUE);

  ch2->Project("histp_F","Phi_Nuc",cut_p);
  ch2->Project("histm_F","Phi_Nuc",cut_m);
	  
  TH1 * BA2=histm_F->GetAsymmetry(histp_F);
  BA2->Scale(1.0/Bpol);
  BA2->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf2 = new TF1("fitf2","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf2->SetParameter(0,0.1);
  fitf2->SetParameter(1,-0.3);
  fitf2->SetParLimits(0,0,0.3);
  fitf2->SetParLimits(1,-1.,1.);
  fitf2->SetLineColor(kBlue);
	  
  BA2->Fit("fitf2");
  BA2->SetLineColor(kBlue);
  BA2->SetMarkerColor(kBlue);
  BA2->SetAxisRange(-0.3, 0.3,"Y");
  BA2->SetAxisRange( 0. ,360.,"X");
  BA2->GetXaxis()->SetTitle("#phi (deg)");
  //BA2->DrawClone("SAME");

  THStack *hist = new THStack("hist","");
  hist->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  hist->Add(BA);
  hist->Add(BA2);
  hist->Draw("nostack");
  box->DrawClone("SAME");
	  
  c1->Print(Folder + "BSA_Fit.png");  
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
