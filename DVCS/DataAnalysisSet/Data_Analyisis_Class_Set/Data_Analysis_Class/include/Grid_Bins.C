#include <fstream>

std::vector<double> BDT::bin_C(int select, int NBINS, TTree*& ch1 )
{
  //1 for Q2, 2 for t, 3 for x_B
  static std::vector<double> *Q2, *t_P, *xB;
  static std::vector<double> bins_temp;
  static std::vector<int> *flag;
  bins_temp.clear();
  std::vector<double> smp;
  int entries=ch1->GetEntries(), DATASZ;

  ch1->SetBranchAddress("strip_Q2",&Q2);
  ch1->SetBranchAddress("strip_Xbj",&xB);
  ch1->SetBranchAddress("t_Ph",&t_P);
  ch1->SetBranchAddress("bestCandidateFlag",&flag);

  //Read Data

  DATASZ=0;
  std::ofstream myfile;
  myfile.open ("Data.csv"); //Change file name each time
  for(int m=0; m<entries;m++)
    {
      ch1->GetEntry(m);
      for (int i=0; i<Q2->size();i++)//any branch has the same size so Q2 here is not special
	{
	  if(select==1 && flag->at(i)==1){myfile<<Q2->at(i)<<endl;DATASZ+=1;}
	  if(select==2 && flag->at(i)==1){myfile<<t_P->at(i)<<endl;DATASZ+=1;}
	  if(select==3 && flag->at(i)==1){myfile<<xB->at(i)<<endl;DATASZ+=1;}
	  
	}
    }
  myfile.close();
  smp.clear();
  int delta=DATASZ/NBINS;
   //Organize data
  gSystem->Exec("math -run < Sort_BSA.m");  

  //Now we create the bins
  std::string filename("Data1.csv");
  std::ifstream input_file(filename);
  string line;
  int m=0,i=0;
  while (getline(input_file, line))
    {
      if(i==delta*m)
	{
	  bins_temp.push_back(std::stod(line));
	  m+=1;
	}
      //std::stod(line);
      //smp.push_back(std::stod(line));
      i+=1;
      if(bins_temp.size()<NBINS+1 && i==DATASZ)
	{
	  bins_temp.push_back(std::stod(line));
	}
    }

  //bins_temp.push_back(std::stod(line));
  
  gSystem->Exec("rm Data.csv");   
  gSystem->Exec("rm Data1.csv");

  return bins_temp;
 
}


std::vector<double> BDT::bin_C2(int select, int delta, TTree*& ch1 )
{
  //1 for Q2, 2 for t, 3 for x_B

  static std::vector<double> *Q2, *t_P, *xB;
  static std::vector<double> bins_temp;
  static std::vector<int> *flag;
  bins_temp.clear();
  std::vector<double> smp;
  int entries=ch1->GetEntries(), DATASZ;

  ch1->SetBranchAddress("strip_Q2",&Q2);
  ch1->SetBranchAddress("strip_Xbj",&xB);
  ch1->SetBranchAddress("t_Ph",&t_P);
  ch1->SetBranchAddress("bestCandidateFlag",&flag);

  //Read Data

  DATASZ=0;
  std::ofstream myfile;
  myfile.open ("Data.csv"); //Change file name each time
  for(int m=0; m<entries;m++)
    {     
      ch1->GetEntry(m);
      for (int i=0; i<Q2->size();i++)//any branch has the same size so Q2 here is not special
	{
	  if(select==1 && flag->at(i)==1){myfile<<Q2->at(i)<<endl;DATASZ+=1;}
	  if(select==2 && flag->at(i)==1){myfile<<t_P->at(i)<<endl;DATASZ+=1;}
	  if(select==3 && flag->at(i)==1){myfile<<xB->at(i)<<endl;DATASZ+=1;}
	  
	}
    }
  myfile.close();
  smp.clear();
  int NBINS=DATASZ/delta +1; //the +1 fixes the error
   //Organize data
  gSystem->Exec("math -run < Sort_BSA.m");  

  //Now we create the bins
  std::string filename("Data1.csv");
  std::ifstream input_file(filename);
  string line;
  int m=0,i=0;
  while (getline(input_file, line))
    {
      if(i==delta*m)
	{
	  bins_temp.push_back(std::stod(line));
	  m+=1;
	}
      //std::stod(line);
      //smp.push_back(std::stod(line));
      i+=1;
      if(bins_temp.size()<NBINS+1 && i==DATASZ)
	{
	  bins_temp.push_back(std::stod(line));
	}
    }

  //bins_temp.push_back(std::stod(line));
  
  gSystem->Exec("rm Data.csv");   
  gSystem->Exec("rm Data1.csv");

  return bins_temp;
 
}


vector<vector<double>> BDT::Grid_Bins(int ii, double *bins_t, int NBINS_t, int NumEv, int NBINS_x, TTree*& Tree)
{
  //i = i-th t bin
  //  int i=0;

  std::vector<vector<double>> bins;
  std::vector<double> binsQ;
  std::vector<double> binst;
  std::vector<double> binsx;
  
  double *bins_Q;
  double *bins_x;
  
  
  TH1F *temp = new TH1F("temp","Histogram",10,-20.0,0.0);
  Tree->Project("temp", "t_Ph",cut);
  int DATASZ=temp->GetEntries();
  delete temp;

  int NBINS=DATASZ/(NumEv*NBINS_t);
  int NBINS_Q;
  if(NBINS%NBINS_x >0) {NBINS_Q =(NBINS/NBINS_x) + 1;}
  if(NBINS%NBINS_x==0) {NBINS_Q =(NBINS/NBINS_x);}
  std::cout<<NBINS_Q<<endl;
  int NumEv_Corr=DATASZ/(NBINS*NBINS_t);  
  
  TCut cutT=TCut(Form("bestCandidateFlag==1 && t_Ph>%f && t_Ph<%f",bins_t[ii],bins_t[ii+1]));
  TFile *file1 = new TFile("file1.root", "RECREATE");
  TTree *TreeT = Tree->CopyTree(cutT);
 
  //Create Q2 bins
  //bin_C2 requires number of events rather than number of bins
  binsQ=bin_C2(1,(NumEv_Corr+1)*NBINS_x,TreeT);
  bins.push_back(binsQ);
  bins_Q=binsQ.data();
  std::cout<<"WARNING: Number of bins in Q2 is "<<binsQ.size() -1<<". Expected: "<<NBINS_Q<<endl;
  std::cout<<"Total number of bins: "<<NBINS<<". Bins in x: "<<NBINS_x<<endl;
  delete TreeT;
  file1->Close();
  delete file1;
  gSystem->Exec("rm file1.root");

  for(int k=0; k<NBINS_Q;k++) //-1 to not include last bins
    {
      TCut cutQ=TCut(Form("bestCandidateFlag==1 && t_Ph>%f && t_Ph<%f && strip_Q2>%f && strip_Q2<%f",bins_t[ii],bins_t[ii+1],bins_Q[k],bins_Q[k+1]));
      TFile *file2 = new TFile("file2.root", "RECREATE");
      TTree *TreeTQ = Tree->CopyTree(cutQ);
      std::cout<<k<<" "<<Tree->GetEntries()<<" "<<TreeTQ->GetEntries()<<endl;
      //Create Xbj bins
      if(k==NBINS_Q-1 && NBINS%NBINS_x>0)
	{
	  binsx=bin_C(3,NBINS%NBINS_x,TreeTQ);
	}
      else
	{
	  binsx=bin_C(3,NBINS_x,TreeTQ);	  
	}
      bins.push_back(binsx);
      gSystem->Exec("rm file2.root");
      delete TreeTQ;
      file2->Close();
      delete file2;
    }


  return bins ;
  
}
