#include <vector>

void BDT::Add_BDT_var_float(TCut cutSB, TString Data, TString output)
{

  cout << "-Start TMVA-" << endl;

  cout << " Adding Variables to Reader" << endl;
  // --- Create the Reader object
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  static vector<Float_t> _vars1(Vars.size());
  static vector<Double_t> _vars(Vars.size());  
  static Float_t Ph_Th;
  static Float_t Nuc_Th;
  double Ph_Th_V;
  double Nuc_Th_V;

  TChain *signal= new TChain("pDVCS");
  signal->Add(Data);

  //==================================

    for(int l=0; l<Vars.size(); l++)
    {
 	reader->AddVariable(Vars.at(l), &_vars1.at(l));
    }
   reader->AddSpectator( "strip_Ph_Theta", &Ph_Th);
   reader->AddSpectator( "strip_Nuc_Theta", &Nuc_Th);


  cout << " Booking Method" << endl;
  //reader->BookMVA("BDT method", Folder + TString("dataset/weights/TMVAClassification_BDT.weights.xml"));
  if(!nodata)
  {
  if(categories)
  	reader->BookMVA("BDT method", Folder + TString("dataset/weights/TMVAClassificationCategory_Ph_Topology.weights.xml"));
  else
  	reader->BookMVA("BDT method", Folder + TString("dataset/weights/TMVAClassificationCategory_BDT.weights.xml"));
  }


  //========================================
  //signal->SetBranchStatus("*", 1); // enable all branches
  
  TFile f(Folder + output, "RECREATE");
  TTree *Tree = signal->CopyTree(cutSB);

  if (!Tree) {
    std::cerr << "No events...exiting" << std::endl;
    
delete reader;
delete signal;
 delete Tree;
 f.Close();
    return;
}

Long64_t nentries = Tree->GetEntries();
  cout << "Reading File with " << nentries << " to Add MVA response" << endl;
  cout << " Reading Variables" << endl;

    for(int l=0; l<Vars.size(); l++)
    {
      Tree->SetBranchAddress(Vars.at(l), &_vars.at(l));
    }
  Tree->SetBranchAddress("strip_Ph_Theta", &Ph_Th_V);
  Tree->SetBranchAddress("strip_Nuc_Theta", &Nuc_Th_V);

  cout << " Loop on all events " << Tree->GetEntries() << " to compute MVA output" << endl;
  
  double _strip_Nuc_BDT;
  

  TBranch *newBranch = Tree->Branch("_strip_Nuc_BDT", &_strip_Nuc_BDT, "_strip_Nuc_BDT/D");

  ////////////////////////////


  for (Long64_t ievt = 0; ievt < Tree->GetEntries(); ievt++)
  {
    Tree->GetEntry(ievt);
    printProgress(ievt*1.0/Tree->GetEntries());
    for(int l=0; l<Vars.size(); l++)
      {
	_vars1.at(l) = _vars.at(l);
      }
	Ph_Th = Ph_Th_V;
	Nuc_Th = Nuc_Th_V;
    
    // Get TMVA output
    if( (categories || (!categories && Ph_Th > 5)) && !nodata )
    {
    _strip_Nuc_BDT = reader->EvaluateMVA("BDT method");
    }
    else
    {
      _strip_Nuc_BDT = 2;
    }
  
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
