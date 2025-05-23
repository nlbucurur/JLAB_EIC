#include <vector>

void BDT::Add_BDT_var_float_Formula(TCut cutSB, TString Data, TString output)
{

  cout << "-Start TMVA-" << endl;

  cout << " Adding Variables to Reader" << endl;
  // --- Create the Reader object
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  TChain *signal= new TChain("pDVCS");
  signal->Add(Data);


  static Float_t _t_Ph1;
  static Float_t _mm2_eg1;
  static Float_t _delta_Phi1;
  static Float_t _delta_t1;
  static Float_t _theta_gamma_X1;
  
  static Double_t _t_Ph;
  static Double_t _mm2_eg;
  static Double_t _delta_Phi;
  static Double_t _delta_t;
  static Double_t _theta_gamma_X;

  reader->AddVariable("t_Ph", &_t_Ph1);
  reader->AddVariable("abs(mm2_eg-0.88)", &_mm2_eg1);
  reader->AddVariable("abs(delta_Phi)", &_delta_Phi1);
  reader->AddVariable("abs(delta_t)", &_delta_t1);
  reader->AddVariable("theta_gamma_X", &_theta_gamma_X1);

  
  cout << " Booking Method" << endl;

  reader->BookMVA("BDT method", Folder + TString("dataset/weights/TMVAClassification_BDT.weights.xml"));
  
  Long64_t nentries = signal->GetEntries();
  cout << "Reading File with " << nentries << " to Add MVA response" << endl;
  cout << " Reading Variables" << endl;


  //========================================
  //signal->SetBranchStatus("*", 1); // enable all branches
  
  TFile f(Folder + output, "RECREATE");

  TTree *Tree = signal->CopyTree(cutSB);

  Tree->SetBranchStatus("*", 1); // enable all branches

  Tree->SetBranchAddress("t_Ph", &_t_Ph);
  Tree->SetBranchAddress("mm2_eg", &_mm2_eg);
  Tree->SetBranchAddress("delta_Phi", &_delta_Phi);
  Tree->SetBranchAddress("delta_t", &_delta_t);
  Tree->SetBranchAddress("theta_gamma_X", &_theta_gamma_X);

  cout << " Loop on all events " << Tree->GetEntries() << " to compute MVA output" << endl;
  
  double _strip_Nuc_BDT;
  

  TBranch *newBranch = Tree->Branch("_strip_Nuc_BDT", &_strip_Nuc_BDT, "_strip_Nuc_BDT/D");

  ////////////////////////////


  for (Long64_t ievt = 0; ievt < Tree->GetEntries(); ievt++)
  {
    Tree->GetEntry(ievt);
    printProgress(ievt*1.0/signal->GetEntries());
	    _t_Ph1 = _t_Ph;
	    _mm2_eg1 = abs(_mm2_eg - 0.88);
	    _delta_Phi1 = abs(_delta_Phi);
	    _delta_t1 = abs(_delta_t);
	    _theta_gamma_X1 = _theta_gamma_X;
    
    // Get TMVA output
    _strip_Nuc_BDT = reader->EvaluateMVA("BDT method");
    
    newBranch->Fill();
  }
  cout << "Writing file" << endl;
  Tree->Write();
  cout << "closing file" << endl;

  cout << "Added " << Tree->GetEntries() << " MVA variables to file" << endl;

  delete reader;
  delete signal;
  delete Tree;
  f.Close();
  
}
