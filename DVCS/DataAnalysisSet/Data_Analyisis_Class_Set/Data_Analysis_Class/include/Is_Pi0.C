#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


void BDT::Is_Pi0()//TString Data, TString output)
{
  TString Data=Folder + "Data_NP_Theta_g_5.root";
  TString output = Folder + "Data_Pi0flag.root";
  // Open the first file and get the tree
  TFile *file1 = TFile::Open(Data);
  TTree *tree1 = (TTree*)file1->Get("pDVCS");

  // Create the output file and tree
  TFile *outfile = new TFile(output,"RECREATE");
  TTree *outtree = tree1->CopyTree("bestCandidateFlag==1");

  // Define variables to hold the branch values
  // Assuming that the first tree has branch "event" and the second tree also has a branch "event"
  static int event1, value;
  static int N_Ph;
  int count=0;
  // Set the branches to the tree1 and tree2
  tree1->SetBranchAddress("EventNumber", &event1);
  tree1->SetBranchAddress("N_Ph", &N_Ph);

  // Set the output tree branches
  TBranch* IsPi0 = outtree->Branch("Is_Pi0", &value);

  std::string filename("Pi0_in_Data.csv");
  string line;
  int event2;

  // Loop over the events in tree1
  for (Long64_t ievent1 = 0; ievent1 < tree1->GetEntries(); ++ievent1)
    {
      value=0;
      // Get the event from tree1
      printProgress(ievent1*1.0/tree1->GetEntries());
      tree1->GetEntry(ievent1);
      if(N_Ph>1)
	{
	  std::ifstream input_file(filename);
	  // Loop over the events tagged
	  while (std::getline(input_file, line))
	    {
	      std::istringstream ss(line);
	      ss >> event2;
	      if(event1==event2)
		{
		  count++;
		  value=1;
		}
	    }
	  input_file.close();
	}
      IsPi0->Fill();
    }

  // Write the output tree to the output file
  outfile->Write();

  // Clean up
  file1->Close();
  outfile->Close();
  //gSystem->Exec("rm Pi0_in_Data.csv");
  std::cout<<"count of Pi0 events "<<count<<" , Expected "<<gSystem->Exec("cat Pi0_in_Data.csv | wc -l")<<endl;
}
