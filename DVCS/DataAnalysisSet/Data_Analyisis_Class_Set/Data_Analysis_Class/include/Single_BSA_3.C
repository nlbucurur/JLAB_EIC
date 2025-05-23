
void BDT::Single_BSA_3(TString BDT_Data, TString Bin_Data, vector<double> boundaries, int P_bins=12)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Folder + Bin_Data);

  TChain *ch2= new TChain("pDVCS");
  ch2->Add(Folder + BDT_Data);

  TChain *ch3= new TChain("pDVCS");
  ch3->Add("/home/munoz/Documents/RG-A:2/NP/FT/Data_NP_Theta_g_5.root");
  ch3->Add("/home/munoz/Documents/RG-A:2/NP/FD/Data_NP_Theta_g_5.root");

  double meanQ;
  double meant;
  double meanx;
  double meanQ_NP;
  double meant_NP;
  double meanx_NP;

  TF1 *fit;
  
  TCanvas* c1 = new TCanvas("c1","Histograms");
		  
  TCut cut_p = TCut("(bestCandidateFlag==1 && Helicity>0)");
  TCut cut_m = TCut("(bestCandidateFlag==1 && Helicity<0)");

  //--------------------------------------
  //Get Means

  TCut cutM = TCut("bestCandidateFlag==1");
	  
  TH1F *tmean = new TH1F("tmean","t_temp",100,0,0);
  ch1->Project("tmean","t_Ph",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meant=tmean->GetMean();

  TH1F *Qmean = new TH1F("Qmean","Q_temp",100,0,0);
  ch1->Project("Qmean","strip_Q2",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanQ=Qmean->GetMean();

  TH1F *xmean = new TH1F("xmean","x_temp",100,0,0);
  ch1->Project("xmean","strip_Xbj",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanx=xmean->GetMean();


  TH1F *tmean_NP = new TH1F("tmean_NP","t_temp_NP",100,0,0);
  ch3->Project("tmean_NP","t_Ph",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meant_NP=tmean_NP->GetMean();

  TH1F *Qmean_NP = new TH1F("Qmean_NP","Q_temp_NP",100,0,0);
  ch3->Project("Qmean_NP","strip_Q2",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanQ_NP=Qmean_NP->GetMean();

  TH1F *xmean_NP = new TH1F("xmean_NP","x_temp_NP",100,0,0);
  ch3->Project("xmean_NP","strip_Xbj",cutM + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  meanx_NP=xmean_NP->GetMean();
  
  delete tmean;
  delete Qmean;
  delete xmean;

  delete tmean_NP;
  delete Qmean_NP;
  delete xmean_NP;
//----------------------------------------

  TH1F *histp = new TH1F("histp","1dgraph",P_bins,0,360);
  TH1F *histm = new TH1F("histm","1dgraph",P_bins,0,360);

  TH1F *histp_F = new TH1F("histp_F","1dgraph",P_bins,0,360);
  TH1F *histm_F = new TH1F("histm_F","1dgraph",P_bins,0,360);

  TH1F *histp_NP = new TH1F("histp_NP","1dgraph",P_bins,0,360);
  TH1F *histm_NP = new TH1F("histm_NP","1dgraph",P_bins,0,360);

  ch1->Project("histp","Phi_Nuc",cut_p);
  ch1->Project("histm","Phi_Nuc",cut_m);
	  
  TH1 * BA=histm->GetAsymmetry(histp);
  BA->Scale(1.0/Bpol);
  BA->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));


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
  BA->SetMarkerStyle(20);
  BA->SetMarkerSize(1);
  BA->SetLineWidth(2);

  BA->SetAxisRange(-0.3, 0.3,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi (deg)");
  //BA->DrawClone("SAME");
	  
  //Display means
  TLatex *box = new TLatex(.6,.8,Form("#splitline{P(NP)}{#splitline{#LT t #GT= %.3f(%.3f)}{#splitline{#LT Q^{2} #GT= %.3f(%.3f)}{#LT x_{B} #GT= %.3f(%.3f)}}}",meant,meant_NP,meanQ,meanQ_NP,meanx,meanx_NP));
  box->SetNDC(kTRUE);
  //------------------------------------------------------------------------------------------------------------------------
  ch2->Project("histp_F","Phi_Nuc",cut_p + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
  ch2->Project("histm_F","Phi_Nuc",cut_m + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
	  
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
  BA2->SetMarkerStyle(20);
  BA2->SetMarkerSize(1);
  BA2->SetLineWidth(2);
  BA2->SetAxisRange(-0.3, 0.3,"Y");
  BA2->SetAxisRange( 0. ,360.,"X");
  BA2->GetXaxis()->SetTitle("#phi (deg)");
  //BA2->DrawClone("SAME");


  //------------------------------------------------------------------------------------------------------------------------
    ch3->Project("histp_NP","Phi_Nuc",cut_p + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
    ch3->Project("histm_NP","Phi_Nuc",cut_m + TCut(Form("t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5))));
	  
  TH1 * BA3=histm_NP->GetAsymmetry(histp_NP);
  BA3->Scale(1.0/Bpol);
  BA3->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf3 = new TF1("fitf3","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf3->SetParameter(0,0.1);
  fitf3->SetParameter(1,-0.3);
  fitf3->SetParLimits(0,0,0.3);
  fitf3->SetParLimits(1,-1.,1.);
  fitf3->SetLineColor(kRed);
	  
  BA3->Fit("fitf3");
  BA3->SetLineColor(kRed);
  BA3->SetMarkerColor(kRed);
  BA3->SetMarkerStyle(20);
  BA3->SetMarkerSize(1);
  BA3->SetLineWidth(2);

  BA3->SetAxisRange(-0.3, 0.3,"Y");
  BA3->SetAxisRange( 0. ,360.,"X");
  BA3->GetXaxis()->SetTitle("#phi (deg)");
  //------------------------------------------------------------------------------------------------------------------------
  THStack *hist = new THStack("hist","");
  hist->SetTitle(Form("%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5)));

  hist->Add(BA);
  //hist->Add(BA2);
  hist->Add(BA3);
  hist->Draw("nostack");
  box->DrawClone("SAME");
	  
  c1->Print(Folder + "P_vs_NP_vs_bin.png");  
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
