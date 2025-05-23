
#include <sstream>
void BDT::FilterP(TCut TheCut){
    
  TChain *chain= new TChain("pDVCS");
  chain->Add("Quality2_Data_P.root");
    
  
  TTree* Tree = chain ;
  TFile f(Folder + TString("Data_P_Theta_g_5.root"),"RECREATE");    

  printf("... copying tree\n");
  TTree* sTree = Tree->CopyTree(TheCut);
  sTree->SetMaxTreeSize(4000000000LL);
    
  printf("... tree copied ... \n");
  sTree->Write();
  f.Close();
    
  return;
}




