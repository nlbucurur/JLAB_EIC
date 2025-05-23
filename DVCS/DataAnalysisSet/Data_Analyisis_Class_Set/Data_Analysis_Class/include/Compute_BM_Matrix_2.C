int BDT::get_index(double t, double Q, double x, double p)
{
  int index=0;
  for(int k=0; k<64; k++)
    {
      if(k>0)
	index+=Nphibins[k-1];
	  
      if(t > bins[k][0] && t < bins[k][1] && Q > bins[k][2] && Q < bins[k][3] && x > bins[k][4] && x < bins[k][5])
	{
	  for(int l=0; l<Nphibins[k]; l++)
	    {
	      if(p > (360./Nphibins[k])*l && p < (360./Nphibins[k])*(l+1) )
		{
		  index=index + l;
		  return index;
		}
	    }
	}
    }
  std::cout<<"Error"<<endl;
  return -1;
}

void BDT::Compute_BM_Matrix() 
{
  TChain *Data= new TChain("pDVCS");
  Data->Add(Folder + TString("Tested_BM_Sim.root"));

  const int N=1186;
  int Hel;
  static vector<double> *t;
  static vector<double> *Q;
  static vector<double> *x;
  static vector<double> *p;

  static vector<double> *t_MC;
  static vector<double> *Q_MC;
  static vector<double> *x_MC;
  static vector<double> *p_MC;

  static double entriesvec[N][N]={0};
//  static double matrix[N][N]={0};
 
  static vector<int> *flag;
  static vector<double> *xsec;
  static double BDTvar;

  Data->SetBranchAddress("t_Ph",&t);
  Data->SetBranchAddress("strip_Q2",&Q);
  Data->SetBranchAddress("strip_Xbj",&x);
  Data->SetBranchAddress("Phi_Ph",&p);

  Data->SetBranchAddress("t_Ph_MC",&t_MC);
  Data->SetBranchAddress("strip_Q2_MC",&Q_MC);
  Data->SetBranchAddress("strip_Xbj_MC",&x_MC);
  Data->SetBranchAddress("Phi_Ph_MC",&p_MC);

  Data->SetBranchAddress("xsec_rc",&xsec);

  Data->SetBranchAddress("Helicity",&Hel);
  Data->SetBranchAddress("bestCandidateFlag",&flag);
  Data->SetBranchAddress("_strip_Nuc_BDT",&BDTvar);
  
  int col, row;
  int entries=Data->GetEntries();
  TH2D* hMatrix = new TH2D("hMatrix", "hMatrix",N,1,N,N,1,N);
  TH1D* input_BSAs = new TH1D("input_BSAs", "input_BSAs",N,1,N);
  
  int input_index=1;
  double BSA, Err, nentriesbin;
  char comma;
  std::unordered_set<int> to_remove;

  //Fill input BSAs
  for(int i=1; i<=64; i++)
  {
	  std::ifstream inputFile(Folder + Form("bin_%i/BSA_Most_Values.txt",i));
	  if (!inputFile.is_open()) {
		  std::cerr << "Error: Unable to open "<<Folder + Form("bin_%i/BSA_Most_Values.txt",i)<< std::endl;
		  }
	 
	 while (inputFile >> nentriesbin >> comma >> BSA >> comma >> Err) 
    {
		input_BSAs->SetBinContent(input_index, nentriesbin);
		input_BSAs->SetBinError(input_index, sqrt(nentriesbin));
		if(Err==0) to_remove.insert(input_index);
		input_index++;
    }
    inputFile.close();	  
  }

  //Fill Matrix
  for(int i=0; i<entries; i++)
    {
      printProgress(i*1.0/entries);
      Data->GetEntry(i);
      for(int j=0; j<flag->size(); j++)
	{
	  if(flag->at(j)==1 && BDTvar > BDT_value )
	    {
	      //Get index of the event, rows are generated, columns are reconstructed
	      //Generated has to be in columns (left to right), reconstructed in rows (top to bottom)
	      row=get_index(t->at(j), Q->at(j), x->at(j), p->at(j));
	      col=get_index(t_MC->at(j), Q_MC->at(j), x_MC->at(j), p_MC->at(j));  
		
		entriesvec[row][col]++;
		//Verify the definition needed for TSVDUnfold
		hMatrix->AddBinContent(row+1, col+1);
	    }
	}
    }

  //Fill MC true and recon BSAs
  TH1D* true_MC = new TH1D("true_MC", "true_MC",N,1,N);
  TH1D* recon_MC= new TH1D("recon_MC", "recon_MC",N,1,N);
  true_MC->Sumw2();
  recon_MC->Sumw2();

  input_index=1;
  for(int i=0; i<64; i++)
    {
      printProgress(i*1.0/64);
		TH1D* aux1_p= new TH1D("aux1_p", "aux1_p",Nphibins[i],0,360);
		TH1D* aux1_m= new TH1D("aux1_m", "aux1_m",Nphibins[i],0,360);
		TH1D* aux2_p= new TH1D("aux2_p", "aux2_p",Nphibins[i],0,360);
		TH1D* aux2_m= new TH1D("aux2_m", "aux2_m",Nphibins[i],0,360);
		
		TCut cut = TCut(Form("bestCandidateFlag==1 && _strip_Nuc_BDT > %f && t_Ph > %f && t_Ph < %f && strip_Q2 > %f && strip_Q2 < %f && strip_Xbj > %f && strip_Xbj < %f", BDT_value, bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]));
		TCut cut_MC = TCut(Form("bestCandidateFlag==1 && _strip_Nuc_BDT > %f && t_Ph_MC > %f && t_Ph_MC < %f && strip_Q2_MC > %f && strip_Q2_MC < %f && strip_Xbj_MC > %f && strip_Xbj_MC < %f", BDT_value, bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]));
		Data->Project("aux1_p","Phi_Ph", cut + TCut("Helicity>0"));
		Data->Project("aux1_m","Phi_Ph", cut + TCut("Helicity<0"));
		Data->Project("aux2_p","Phi_Ph_MC", cut_MC + TCut("Helicity>0"));
		Data->Project("aux2_m","Phi_Ph_MC", cut_MC + TCut("Helicity<0"));
		
		TH1 *BA1= aux1_p->GetAsymmetry(aux1_m);
		TH1 *BA2= aux2_p->GetAsymmetry(aux2_m);

	for(int j=1; j<=Nphibins[i];j++)
	{
		if(aux1_p->GetBinContent(j) + aux1_m->GetBinContent(j)==0) 
		{
			BA1->SetBinContent(j,0);
			BA1->SetBinError(j,0);
		}
		if(aux2_p->GetBinContent(j) + aux2_m->GetBinContent(j)==0) 
		{
			BA2->SetBinContent(j,0);
			BA2->SetBinError(j,0);
		}
	}

	for(int j=1; j<=Nphibins[i];j++)
	{
		recon_MC->SetBinContent(input_index,aux1_p->GetBinContent(j));
		true_MC ->SetBinContent(input_index,aux2_p->GetBinContent(j));
		recon_MC->SetBinError(input_index,aux1_p->GetBinError(j));
		true_MC ->SetBinError(input_index,aux2_p->GetBinError(j));
		input_index++;
	}
	delete aux1_p;
	delete aux1_m;
	delete aux2_p;
	delete aux2_m;
	delete BA1;
	delete BA2;
	}	

  //Normalize the matrix

  double sum, esum;
  int lele=0;
  double val;
  /*
  for(int i=1; i<=N; i++)
    {
      sum=0;
      esum=0;
      //Get normalization
      //Sum_j of f^{i}_{j} should be 1 
      for(int j=1; j<=N; j++)
	{
	  esum+=entriesvec[i][j];
	  //sum+=matrix[i][j];
	}

      //The (i,j) entry must give the probability of being reconstructed in the i bin, given that it was generated on j (nrec_i/ngen_j).
      //The sum of entries of a given row must add to 1, as given a generated event it has to be reconstructed somewhere
      for(int j=1; j<=N; j++)
	{
	  //Add diagonal element for bins with no generated events
	  if(esum<20)// || true_MC->GetBincontent(j)==0 || recon_MC->GetBincontent(j)==0 || input_BSAs->GetBincontent(j)==0)
	  {
		lele++;
		val = (i==j) ? 1 : 0;
	    hMatrix->SetBinContent(i,j, val);
	  }
	  else
	    hMatrix->SetBinContent(i,j, hMatrix->GetBinContent(i,j)/esum);
	 }
    }	
		  std::cout<<"\n Bins with less than 20 events: "<<lele/entries<<endl;

 */
 
  //Remove bins with null BSA measurement
  //The covariance matrix is built with the error of the measurements
  //A null measurement is a null eigenvalue in the matrix, so it cannot be inverted
  int N2=N - to_remove.size();
  TH2D* hMatrix_filtered = new TH2D("hMatrix_filtered", "hMatrix_filtered",N2,1,N2,N2,1,N2);
  TH1D* input_BSAs_filtered = new TH1D("input_BSAs_filtered", "input_BSAs_filtered",N2,1,N2);
  TH1D* true_MC_filtered = new TH1D("true_MC_filtered", "true_MC_filtered",N2,1,N2);
  TH1D* recon_MC_filtered= new TH1D("recon_MC_filtered", "recon_MC_filtered",N2,1,N2);

  int input_i=1;
  int input_j=1;

  for(int i=1; i<=N; i++)
    {
    if (to_remove.find(i) == to_remove.end())
    {
		input_BSAs_filtered->SetBinContent(input_i,input_BSAs->GetBinContent(i));
		true_MC_filtered->SetBinContent(input_i,true_MC->GetBinContent(i));
		recon_MC_filtered->SetBinContent(input_i,recon_MC->GetBinContent(i));

		input_BSAs_filtered->SetBinError(input_i,input_BSAs->GetBinError(i));
		true_MC_filtered->SetBinError(input_i,true_MC->GetBinError(i));
		recon_MC_filtered->SetBinError(input_i,recon_MC->GetBinError(i));
		
		input_j=1;
		for(int j=1; j<=N; j++)
		{
			if (to_remove.find(j) == to_remove.end())
			{
				hMatrix_filtered->SetBinContent(input_i,input_j, hMatrix->GetBinContent(i,j));
				hMatrix_filtered->SetBinError(input_i,input_j, hMatrix->GetBinError(i,j));
				input_j++;
			}
		}
		input_i++;
	}
	}


  TSVDUnfold* BMunfold = new TSVDUnfold(input_BSAs_filtered, recon_MC_filtered, true_MC_filtered, hMatrix_filtered);
  BMunfold->SetNormalize( kFALSE ); // no normalisation here
  TH1D* output_BSAs = BMunfold->Unfold( 13 );
  
/*
  TUnfold BMunfold(hMatrix_filtered,TUnfold::kHistMapOutputVert);
  BMunfold.SetInput(input_BSAs_filtered);
  
  // scan L curve and find best point
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  iBest=BMunfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

  //BMunfold.DoUnfold(1.0);
  BMunfold.GetOutput(output_BSAs);
*/
  
  //Write output
  std::ofstream outFile;
  outFile.open(Folder + TString("factor_mig_most.csv"));
  input_index=1;
  for(int i=1; i<=N;i++)
  {
    if (to_remove.find(i) == to_remove.end())
    {
	  outFile<<output_BSAs->GetBinContent(input_index)/input_BSAs_filtered->GetBinContent(input_index)<<" "<<hMatrix_filtered->GetBinContent(input_index,input_index)<<endl;
	  input_index++;
	}
	else
	  outFile<<1.0<<endl;	
  }


}
