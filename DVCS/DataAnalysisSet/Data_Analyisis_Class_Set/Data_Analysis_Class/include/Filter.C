
#include <sstream>
void BDT::Filter(TString TData, TCut TheCut, TString output){
    
  TChain *chain= new TChain("pDVCS");
  chain->Add(Folder+TData);
    
  
  TTree* Tree = chain ;
  TFile f(Folder + output,"RECREATE");    

  printf("... copying tree\n");
  TTree* sTree = Tree->CopyTree(TheCut);
  sTree->SetMaxTreeSize(4000000000LL);
    
  printf("... tree copied ... \n");
  sTree->Write();
  f.Close();

  //delete sTree;
  //delete Tree;
  //delete chain;
  return;
}




