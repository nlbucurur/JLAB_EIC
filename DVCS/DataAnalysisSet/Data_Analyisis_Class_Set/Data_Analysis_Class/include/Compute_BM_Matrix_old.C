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

  std::ofstream outFile;
  outFile.open(Folder + TString("BM_Matrix.csv"));

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
  static double matrix[N][N]={0};
 
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
		matrix[row][col]++;//=xsec->at(j);

	    }
	}
    }

  double sum, esum;
  int lele=0;
  
  for(int i=0; i<N; i++)
    {
      sum=0;
      esum=0;
      //Get normalization
      //Sum_j of f^{i}_{j} should be 1 
      for(int j=0; j<N; j++)
	{
	  esum+=entriesvec[i][j];
	  sum+=matrix[i][j];
	}

      //The (i,j) entry must give the probability of being reconstructed in the i bin, given that it was generated on j (nrec_i/ngen_j).
      //The sum of entries of a given row must add to 1, as given a generated event it has to be reconstructed somewhere
      for(int j=0; j<N; j++)
	{
	  //Add diagonal element for bins with no generated events
	  if(esum<20)
	  {
		  lele++;
	    matrix[i][j] = (i==j) ? 1 : 0;
	  }
	  else
	    matrix[i][j]=matrix[i][j]/esum;
	 }
    }	
		  std::cout<<"\n Bins with less than 20 events: "<<lele<<endl;


//Write output
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  if(j==N-1)
	      outFile<<matrix[i][j];
	  else
	      outFile<<matrix[i][j]<<",";
	}
      outFile<<endl;
	}
  
  outFile.close();
  //Notice that before the inversion, Sum P(R|G)N_G = N_R  , the matrix needs to have the generated over the columns, not the rows.
  
  gSystem->Exec(TString("cp include/Bin_Migrate.m ") + Folder);
  gSystem->Exec(TString("cp include/collect_data.sh ") + Folder);
  gSystem->cd(Folder);
  gSystem->Exec("rm factor_maxi_mig.csv");
  gSystem->Exec("rm factor_most_mig.csv");
  gSystem->Exec("chmod 755 collect_data.sh && ./collect_data.sh");
  gSystem->Exec("math -run < Bin_Migrate.m");  
  gSystem->cd("..");
}
