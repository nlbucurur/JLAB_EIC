#include <vector>

void BDT::Add_BDT_var_Formula(TCut cutSB, TString Data, TString output, vector<TString> vars)
{

  cout << "-Start TMVA-" << endl;

  cout << " Adding Variables to Reader" << endl;
  // --- Create the Reader object
  TMVA::Reader *reader1 = new TMVA::Reader("!Color:!Silent");
  
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  static Float_t _t_Ph1;
  static Float_t _mm2_eg1;
  static Float_t _delta_Phi1;
  static Float_t _delta_t1;
  static Float_t _theta_gamma_X1;
  
  static vector<double>* _t_Ph;
  static vector<double>* _mm2_eg;
  static vector<double>* _delta_Phi;
  static vector<double>* _delta_t;
  static vector<double>* _theta_gamma_X;
  static vector<int>* flag;


  reader1->AddVariable("t_Ph", &_t_Ph1);
  reader1->AddVariable("abs(mm2_eg-0.88)", &_mm2_eg1);
  reader1->AddVariable("abs(delta_Phi)", &_delta_Phi1);
  reader1->AddVariable("abs(delta_t)", &_delta_t1);
  reader1->AddVariable("theta_gamma_X", &_theta_gamma_X1);
  
  cout << " Booking Method" << endl;
  reader1->BookMVA("BDT method", Folder + TString("dataset/weights/TMVAClassification_BDT.weights.xml"));

  TFile *input4 = new TFile(Data,"READ");
  TTree *signal;
  if( input4->GetListOfKeys()->Contains("pDVCS") )
    signal     = (TTree*)input4->Get("pDVCS");
  else
    signal     = (TTree*)input4->Get("eppi0");
 
  Long64_t nentries = signal->GetEntries();
  cout << "Reading File with " << nentries << " to Add MVA response" << endl;
  cout << " Reading Variables" << endl;


  // for DVCS channels
  TFile f(Folder + output, "RECREATE");
  TTree *Tree = signal->CopyTree(cutSB);

  Tree->SetBranchStatus("*", 1); // enable all branches

  Tree->SetBranchAddress("t_Ph", &_t_Ph);
  Tree->SetBranchAddress("mm2_eg", &_mm2_eg);
  Tree->SetBranchAddress("delta_Phi", &_delta_Phi);
  Tree->SetBranchAddress("delta_t", &_delta_t);
  Tree->SetBranchAddress("theta_gamma_X", &_theta_gamma_X);
  Tree->SetBranchAddress("bestCandidateFlag", &flag);

  cout << " Loop on all events " << signal->GetEntries() << " to compute MVA output" << endl;
  
  double _strip_Nuc_BDT;

  TBranch *newBranch = Tree->Branch("_strip_Nuc_BDT", &_strip_Nuc_BDT, "_strip_Nuc_BDT/D");

  ////////////////////////////


  for (Long64_t ievt = 0; ievt < Tree->GetEntries(); ievt++)
  {
    Tree->GetEntry(ievt);
    printProgress(ievt*1.0/Tree->GetEntries());
    for(int i=0; i<flag->size();i++)
      {
	if(flag->at(i)==1)
	  {
	    //_strip_Ph_P1=_strip_Ph_P->at(0);
	    //    std::cout<<_strip_Ph_P1<<endl;
	    
	    _t_Ph1 = _t_Ph->at(i);
	    //std::cout<<_t_Ph1<<endl;
	    
	    _mm2_eg1 = abs(_mm2_eg->at(i) - 0.88);
	    //std::cout<<_mm2_eg1<<endl;
	    
	    _delta_Phi1 = abs(_delta_Phi->at(i));
	    //std::cout<<_delta_Phi1<<endl;
	    
	    _delta_t1 = abs(_delta_t->at(i));
	    //std::cout<<_delta_t1<<endl;
	    
	    _theta_gamma_X1 = _theta_gamma_X->at(i);
	    //std::cout<<_theta_gamma_X1<<endl;
	    
	    // Get TMVA output
	    _strip_Nuc_BDT = reader1->EvaluateMVA("BDT method");
	    
	    newBranch->Fill();
	  }
      }
  }
  cout << "Writing file" << endl;
  Tree->Print();
  Tree->Write();
  cout << "closing file" << endl;

  cout << "Added " << Tree->GetEntries() << " MVA variables to file" << endl;
  delete reader1;
  delete signal;
  input4->Close();
  delete input4;
  
  delete Tree;
  f.Close();

 }
