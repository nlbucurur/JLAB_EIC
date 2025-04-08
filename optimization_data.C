#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TClass.h>
#include <TStyle.h>
#include <TROOT.h>

void optimization_data () {

  double P_mass = 0.938272;
  double N_mass = 0.9395654;

  TFile *file = TFile::Open("/home/lorena/Thesis/JLAB_EIC/0pDVCS_inbending_FTPhotonsCorrected_test.root");

  TTree *tree = (TTree*) file->Get("pDVCS_stripped");
 
//   tree->Print();

  std::vector<TString> branch_names = {
      "_mm2_eg", // Expecting missing proton
      "_mm2_eNg", // Expecting neutron mass
      "_mm2_eNg_N", // Expecting nothing
      "_mm2_eNX_N", // Expecting photon
      "_strip_Q2", // Q2 > M_proton for DIS
      "_strip_Xbj", // Valence region
      "_t_Nuc", // Most of the data
      "_t_Ph", // Same as above
      "_delta_t", // Almost 0 for DVCS
      "_Phi_Nuc",
      "_Phi_Ph",
      "_delta_Phi"
  };

    std::vector<TCut> cuts;

    cuts.push_back("_theta_gamma_e > 6");
    
}