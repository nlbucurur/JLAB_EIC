
void BDT::To_bin_or_not_to_bin(TString Data, int NumEv)//, TCut cut, TString DVCS, TString Pi0, vector<TString> vars)
{
  TFile *File = new TFile(Data);
  TTree *ch1 = (TTree *)(File->Get("pDVCS"));
  
  int NBINS_t=3;
  int NBINS_x=3;
  //int NumEv=29000;
  
  TH1F *temp2 = new TH1F("temp2","Histogram",10,0,0);
  ch1->Project("temp2", "t_Ph","bestCandidateFlag==1");
  int DATASZ=ch1->GetEntries();
  int NBINS=DATASZ/(NumEv*NBINS_t);
  int NBINS_Q = NBINS%NBINS_x >0 ? (NBINS/NBINS_x) + 1 : (NBINS/NBINS_x);
  
  int NumEv_Corr=DATASZ/(NBINS*NBINS_t);  

  std::vector<vector<double>> bins;
  std::vector<double> binst;

  double *bins_t;
  double *bins_Q;
  double *bins_x;
  int bin_number=1;
  TString Folder_old=Folder;

  binst=bin_C(2,NBINS_t,ch1);
  bins_t=binst.data();
  std::vector<double> boundaries;
  //Bin_View(Data, bins_t, NumEv);
  for(int k=0;k<NBINS_t;k++)
    {
      bins=Grid_Bins(k, bins_t, NBINS_t, NumEv,NBINS_x, ch1);
      bins_Q=bins.at(0).data();
      for(int m=0;m<NBINS_Q;m++)
	{
	  bins_x=bins.at(m+1).data();
	  for(int n=0;n<min(NBINS_x,NBINS - m*NBINS_Q);n++)
	    {
	      TCut cut_bin = TCut(Form("bestCandidateFlag==1  && t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1]));
	      std::cout<<Form("t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f && strip_Xbj>%f && strip_Xbj<%f",bins_t[k],bins_t[k+1],bins_Q[m],bins_Q[m+1],bins_x[n],bins_x[n+1])<<endl;
	      boundaries.clear();
	      boundaries.push_back(bins_t[k]);
	      boundaries.push_back(bins_t[k+1]);
	      boundaries.push_back(bins_Q[m]);
	      boundaries.push_back(bins_Q[m+1]);
	      boundaries.push_back(bins_x[n]);
	      boundaries.push_back(bins_x[n+1]);
	      Folder = Folder_old + TString("bin_")+Form("%i/",bin_number);
	      gSystem->Exec(TString("mkdir -p ") + Folder);
	      //Make BSA global vs bin_BDT
	      Single_BSA_2("../Data_NP_Theta_g_5.root","Data_NP_Theta_g_5.root", boundaries, 12);

	      bin_number+=1;	      
	    }
	}

      bins.clear();
    }

  Folder = Folder_old;

  delete ch1;
  File->Close();
  delete File;
  delete temp2;
  return;
}
