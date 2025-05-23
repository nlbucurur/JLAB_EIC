
void BDT::Bin_View(TString Data, double *bins_t, TTree*& ch1, int NumEv)
{
  
  int NBINS_t=3;
  int NBINS_x=3;
  //int NumEv=29000;
  
  TH1F *temp3 = new TH1F("temp3","Histogram",10,-20.0,0.0);
  ch1->Project("temp3", "t_Ph",cut);
  int DATASZ=temp3->GetEntries();
  delete temp3;
  
  int NBINS=DATASZ/(NumEv*NBINS_t);
  int NBINS_Q;
  if(NBINS%NBINS_x >0) {NBINS_Q =(NBINS/NBINS_x) + 1;}
  if(NBINS%NBINS_x==0) {NBINS_Q =(NBINS/NBINS_x);}

  int NumEv_Corr=DATASZ/(NBINS*NBINS_t);  


  std::vector<vector<double>> bins;
  std::vector<double> binst;
  std::vector<double> binsQ;
  std::vector<double> binsx;

  double *bins_Q;
  double *bins_x;

  double meanQ;
  double meant;
  double meanx;
  
  TCanvas* c1 = new TCanvas("c1","Histograms",1500,500);
  c1->Divide(NBINS_t,1);

  for(int k=0;k<NBINS_t;k++)
    {
      bins=Grid_Bins(k, bins_t, NBINS_t, NumEv,NBINS_x, ch1);
      c1->cd(k+1);
      bins_Q=bins.at(0).data();
	    
      TH2F *Bin2D = new TH2F("Bin2D","Bin2D",40,bins_Q[0],12,40,0,0.9);
      TCut cutt = TCut(Form("bestCandidateFlag==1  && t_Ph>%f && t_Ph<%f ",bins_t[k],bins_t[k+1]));
      ch1->Project("Bin2D","strip_Xbj:strip_Q2",cut + cutt);
      Bin2D->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
      Bin2D->GetYaxis()->SetTitle("x_{B}");
      Bin2D->SetTitle(Form(" %.3f <t(GeV^{2})< %.3f",bins_t[k],bins_t[k+1]));
      Bin2D->DrawClone("COLZ");
	      
      for(int m=1;m<NBINS_Q;m++)
	{
	  bins_x=bins.at(m).data();
	  TLine *l1 = new TLine(bins_Q[m],0,bins_Q[m],0.9);
	  l1->SetLineWidth(1);

	  l1->DrawClone();
	  for(int n=1;n<NBINS_x;n++)
	    {
	      TLine *l2 = new TLine(bins_Q[m-1],bins_x[n],bins_Q[m],bins_x[n]);
	      l2->SetLineWidth(1);
	      l2->DrawClone();
	      delete l2;
	    }
	  delete l1;
	}
      delete Bin2D;

      bins_x=bins.at(NBINS_Q).data();
      int last_X_bins;
      if(NBINS%NBINS_x > 0){last_X_bins= NBINS%NBINS_x;}
      if(NBINS%NBINS_x== 0){last_X_bins= NBINS_x;}
      
      for(int n=1;n<last_X_bins;n++)
	{
	  TLine *l2 = new TLine(bins_Q[NBINS_Q-1],bins_x[n],12,bins_x[n]);
	  l2->SetLineWidth(1);
	  l2->DrawClone();
	  delete l2;
	}
      
      bins.clear();
    }

  c1->Print(Folder + "Bins.pdf");

  delete c1;  

  return;
}
