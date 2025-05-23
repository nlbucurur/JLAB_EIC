void BDT::Merge_Bins(TString data, TString output=""){
    
  if(output =="")
    output=data;
    
  TChain *Tree= new TChain("pDVCS");
  for(int i=1; i<=Nbins;i++)
  {
    Tree->Add(Folder + Form("bin_%i/",i) + data);
  }
  
    TFile f(Folder + output,"RECREATE");    
    printf("... copying tree\n");
    TTree* sTree = Tree->CopyTree("bestCandidateFlag==1");
    sTree->SetMaxTreeSize(400000000000LL);
    
    printf("... tree merged ... \n");
    sTree->Write();
    f.Close();
    
    return;
}




