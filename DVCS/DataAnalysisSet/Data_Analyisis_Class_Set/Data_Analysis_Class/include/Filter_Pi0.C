
#include <sstream>
void BDT::Filter_Pi0(TString Pi0Data, TCut TheCut, TString output){
    
  //replace Ph by Pi0
  TString String_cut = TheCut.GetTitle();
  String_cut.ReplaceAll("t_Ph", "t_Pi0");
	
  TString l1 = Folder[0];
  TString l2 = Folder[1];
  TString l3 = Folder[2];
  TString temps = l1+l2+l3;
  if(temps == "FT/")
    String_cut.ReplaceAll("strip_Ph_Theta", "mm2_eNgg");//mm2_eNgg for FT and theta_Pi0_e for FD
  if(temps == "FD/")
    String_cut.ReplaceAll("strip_Ph_Theta", "theta_Pi0_e");//mm2_eNgg for FT and theta_Pi0_e for FD
	  
  String_cut.ReplaceAll("strip_Ph_P", "strip_Pi0_P");
  String_cut.ReplaceAll("mm2_eg", "strip_Xbj");
  String_cut.ReplaceAll("gamma", "Pi0");
  //Create cut for 2gamma case
  cut2g = TCut(String_cut) + TCut("strip_Pi0_2DChi2 < 0.04 && abs(mm2_eNgg)<0.01"); // && (strip_Ph1_Theta <5 || strip_Ph2_Theta <5)");


  TChain *chain= new TChain("eppi0");
  chain->Add(Pi0Data);
    
  
  TTree* Tree = chain ;

  TFile f(Folder + output,"RECREATE");    

  printf("... copying tree\n");
  TTree* sTree = Tree->CopyTree(cut2g);
  sTree->SetMaxTreeSize(4000000000LL);
    
  printf("... tree copied ... \n");
  sTree->Write();
  f.Close();
    
  return;
}




