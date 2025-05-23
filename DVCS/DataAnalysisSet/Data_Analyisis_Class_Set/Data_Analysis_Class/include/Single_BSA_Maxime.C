
TH1* BDT::Single_BSA_Maxime(TString Data, vector<double> Pbins)
{
  TChain *ch1= new TChain("pDVCS");
  ch1->Add(Folder + Data);

  double meanQ;
  double meant;
  double meanx;

  TF1 *fit;
  
  TCanvas* c1 = new TCanvas("c1","Histograms");
		  
  TCut cut_p = TCut("bestCandidateFlag==1 && Helicity>0");
  TCut cut_m = TCut("bestCandidateFlag==1 && Helicity<0");

  //--------------------------------------
  //Get Means

  TCut cutM = TCut("bestCandidateFlag==1");
	  
  TH1F *tmean = new TH1F("tmean","t_temp",100,0,0);
  ch1->Project("tmean","t_Ph",cutM);
  meant=tmean->GetMean();

  TH1F *Qmean = new TH1F("Qmean","Q_temp",100,0,0);
  ch1->Project("Qmean","strip_Q2",cutM);
  meanQ=Qmean->GetMean();

  TH1F *xmean = new TH1F("xmean","x_temp",100,0,0);
  ch1->Project("xmean","strip_Xbj",cutM);
  meanx=xmean->GetMean();

  TH1F *histp = new TH1F("histp","1dgraph",Pbins.size()-1,0,360);
  TH1F *histm = new TH1F("histm","1dgraph",Pbins.size()-1,0,360);
  histp->Sumw2();
  histm->Sumw2();
	  
  ch1->Project("histp","Phi_Nuc",cut_p);
  ch1->Project("histm","Phi_Nuc",cut_m);
	  
  TH1 * BA=histm->GetAsymmetry(histp);
  BA->Scale(1.0/Bpol);
  BA->SetTitle(Form("#splitline{%.3f<t<%.3f , %.3f<Q^{2}<%.3f , %.3f<x_{B}<%.3f}{#splitline{BDT_bef %.3f, BDT_aft %.3f}{ #splitline{#frac{1#gamma}{2#gamma}_bef ratio %.3f,  #frac{1#gamma}{2#gamma}_aft ratio %.3f}{Maxime_bef=%.3f, Maxime_aft=%.3f}}}",boundaries.at(0),boundaries.at(1),boundaries.at(2),boundaries.at(3),boundaries.at(4),boundaries.at(5),boundaries.at(6),boundaries.at(7),boundaries.at(8),boundaries.at(9),boundaries.at(10),boundaries.at(11)));

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
	  
  BA->Fit("fitf","Q");
  BA->SetLineColor(kBlack);
  BA->SetMarkerColor(kBlack);
  BA->SetAxisRange(-0.3, 0.3,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi (deg)");
  BA->DrawClone();
	  

  //Display means
  TLatex *box = new TLatex(.6,.6,Form("#splitline{#LT t #GT= %.3f}{#splitline{#LT Q^{2} #GT= %.3f}{#LT x_{B} #GT= %.3f}}",meant,meanQ,meanx));
  box->SetNDC(kTRUE);
  box->DrawClone("SAME");

  //auto *g = new TGraph(Form("Triple%i/Triple%i.txt",k+1,NBINS_Q*l + m+1));
  //g->DrawClone("SAME");

	  
  c1->Print(Folder + "BSA.png");  
  delete histp;
  delete histm;
  delete ch1;
  delete c1;

  return BA;
}
