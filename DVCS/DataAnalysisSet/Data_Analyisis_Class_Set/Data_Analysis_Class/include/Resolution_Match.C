void BDT::Resolution_Match(TCut cut, TString Data, TString Sim){
  cout << "Computing matching weights..." << endl;
  //Run Python script
  gSystem->Exec(TString("mkdir -p ") + Folder + TString("Reweighting_Plots"));
  gSystem->Exec(TString("cp include/get_weights.py ") + Folder);
  gSystem->Exec(TString("cd ") + Folder + TString(" && python3 get_weights.py ") + Sim + TString(" ") + Data + TString(" && cd - && pwd"));
  
  std::ifstream inputFile(Folder + TString("rWeights.dat"));
  if (!inputFile.is_open()) {
    std::cerr << "Error: Unable to open weight file"<< std::endl;
  }

  double rw;
  
  TFile *input = new TFile(Folder + Sim,"READ");
  TTree *signal = (TTree*)input->Get("eppi0");

  signal->SetBranchStatus("*", 1); // enable all branches  
  TFile f(Folder + TString("w") + Sim, "RECREATE");
  TTree *Tree = signal->CopyTree(cut);

  double rWeight;
  TBranch *newBranch = Tree->Branch("rWeight", &rWeight, "rWeight/D");

  for (Long64_t ievt = 0; ievt < Tree->GetEntries(); ievt++)
  {
    Tree->GetEntry(ievt);
    printProgress(ievt*1.0/Tree->GetEntries());
    inputFile >> rw;
    rWeight = rw;
    //std::cout<<rw<<endl;
    newBranch->Fill();
  }
  cout << "Writing file" << endl;
  Tree->Write();
  cout << "closing file" << endl;
  cout << "Added resolution matching weights" << endl;
  inputFile.close();

  //gSystem->Exec(TString("rm ") + Folder + Sim );
  gSystem->Exec(TString("rm ") + Folder + TString("rWeights.dat") );
}
