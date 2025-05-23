void BDT::Check_GPU(bool Pi0_flag=false){

  static vector<int> *Flag1;
  static vector<int> *Flag2;

  static int Ev1;
  static int Ev2;
  
  TString tree;
  if(Pi0_flag) tree="eppi0";
  else  tree="pDVCS";

  TFile *input = new TFile(Folder + TString("Data_P_Theta_g_5.root"),"READ");
  TTree *pDVCS = (TTree*)input->Get(tree);
  TFile *input2 = new TFile(Folder + TString("Data_NP_Theta_g_5.root"),"READ");
  TTree *pDVCS_NP = (TTree*)input2->Get("pDVCS");

  pDVCS->SetBranchStatus("*", 0); // disable all branches
  pDVCS->SetBranchStatus("EventNumber", 1); 
  pDVCS->SetBranchStatus("bestCandidateFlag", 1); 

  pDVCS_NP->SetBranchStatus("*", 0); // disable all branches
  pDVCS_NP->SetBranchStatus("EventNumber", 1);
  pDVCS_NP->SetBranchStatus("bestCandidateFlag", 1);

  pDVCS_NP->SetBranchStatus("EventNumber", 1);
  pDVCS->SetBranchAddress("EventNumber",&Ev1);
  pDVCS_NP->SetBranchAddress("EventNumber",&Ev2);
  
  pDVCS->SetBranchAddress("bestCandidateFlag",&Flag1);
  pDVCS_NP->SetBranchAddress("bestCandidateFlag",&Flag2);

  std::ofstream myfile;
  myfile.open ("Data_P.csv"); //Change file name each time
  for(int m=0; m<pDVCS->GetEntries();m++)
    {     
      pDVCS->GetEntry(m);
      //std::cout<<m<<" "<<Flag1->size()<<endl;
      for(int i=0; i<Flag1->size();i++)
	{
	  if(Flag1->at(i)==1)
	    {
	      myfile<<Ev1<<endl;
	    }
	}
    }
  myfile.close();

  std::ofstream myfile2;
  myfile2.open ("Data_NP.csv"); //Change file name each time
  for(int m=0; m<pDVCS_NP->GetEntries();m++)
    {     
      pDVCS_NP->GetEntry(m);
      for(int i=0; i<Flag2->size();i++)
	{
	  if(Flag2->at(i)==1)
	    {
	      myfile2<<Ev2<<endl;
	    }
	}
    }
  myfile2.close();
  gSystem->Exec("math -run < Sort.m");
  if(!Pi0_flag)
    {
      gSystem->Exec("rm Pi0_flag_in_Data.csv");
    }
  gSystem->Exec("rm Data_P.csv");
  gSystem->Exec("rm Data_NP.csv");

  delete Flag1;
  delete Flag2;

  
  delete input;
  delete input2;
  delete pDVCS;
  delete pDVCS_NP;
  
  
}
