#include <vector>

void BDT::Add_BDT_var2(TString filename)
{
  int nthreads = 4;
  ROOT::EnableImplicitMT(nthreads);

  cout << "-Start TMVA-" << endl;

  cout << " Adding Variables to Reader" << endl;
  // --- Create the Reader object
  TMVA::Reader *reader1 = new TMVA::Reader("!Color:!Silent");
  
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  static vector<Float_t> _vars1(Vars.size());
  static vector<vector<double>*> _vars(Vars.size());  
  static vector<int>* flag;

  TFile *File = new TFile(Folder + filename,"update");
  TTree *signal = (TTree*)File->Get("pDVCS");   
 
  //==================================

    for(int l=0; l<Vars.size(); l++)
    {
      reader1->AddVariable(Vars.at(l), &_vars1.at(l));
      signal->SetBranchAddress(Vars.at(l), &_vars.at(l));

    }
  reader1->BookMVA("BDT method", "Analysis/dataset/weights/TMVAClassification_BDT.weights.xml");
  
  cout << " Booking Method" << endl;

  Long64_t nentries = signal->GetEntries();
  cout << "Reading File with " << nentries << " to Add MVA response" << endl;
  cout << " Reading Variables" << endl;

  //========================================
  signal->SetBranchStatus("*", 1); // enable all branches
  signal->SetBranchAddress("bestCandidateFlag", &flag);


  // for DVCS channels
   TCut cut =TCut("bestCandidateFlag==1");

  cout << " Loop on all events " << signal->GetEntries() << " to compute MVA output" << endl;
  
  double _strip_Nuc_BDT_2;

  TBranch *newBranch = signal->Branch("_strip_Nuc_BDT_2", &_strip_Nuc_BDT_2, "_strip_Nuc_BDT_Data_2/D");

  ////////////////////////////


for (Long64_t ievt = 0; ievt < signal->GetEntries(); ievt++)
  {
    signal->GetEntry(ievt);
    printProgress(ievt*1.0/signal->GetEntries());
    for(int i=0; i<flag->size();i++)
      {
	if(flag->at(i)==1)
	  {

	    for(int l=0; l<Vars.size(); l++)
	      {
		_vars1.at(l) = _vars.at(l)->at(i);
	      }

	    // Get TMVA output
	    _strip_Nuc_BDT_2 = reader1->EvaluateMVA("BDT method");
	    
	    newBranch->Fill();
	  }
      }
  }
    cout << "Writing file" << endl;
  signal->Print();
  signal->Write();
  cout << "closing file" << endl;

  cout << "Added " << signal->GetEntries() << " MVA variables to file" << endl;
  delete reader1;
  delete signal;
  File->Close();
  delete File;
}
