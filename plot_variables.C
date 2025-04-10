#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TClass.h>
#include <TStyle.h>
#include <TROOT.h>

struct branches_cuts {
  std::string name;
  double min;
  double max;
  double cut_min;
  double cut_max;
};

struct branches_cuts_2D {
  std::string name;
  std::string x_branch;
  std::string y_branch;
  double x_min;
  double x_max;
  double y_min;
  double y_max;
};

void stats_legend (TH1D* htemp, TH1D* htemp_cut, const std::string& branch_name) {
  gPad->Update();

  TPaveStats* stats1 = (TPaveStats*)htemp->FindObject("stats");
  TPaveStats* stats2 = (TPaveStats*)htemp_cut->FindObject("stats");

  bool is_t = (branch_name == "_t_Nuc" || branch_name == "_t_Ph");	
  
  if (is_t) {
    stats1->SetX1NDC(0.15); stats1->SetX2NDC(0.33);
    stats1->SetY1NDC(0.75); stats1->SetY2NDC(0.85);
  } else {
    stats1->SetX1NDC(0.80); stats1->SetX2NDC(0.98);
    stats1->SetY1NDC(0.85); stats1->SetY2NDC(0.95);
  }
  stats1->SetTextColor(kBlack);
  
  if (is_t) {
    stats2->SetX1NDC(0.15); stats2->SetX2NDC(0.33);
    stats2->SetY1NDC(0.6); stats2->SetY2NDC(0.7);
  } else {
    stats2->SetX1NDC(0.80); stats2->SetX2NDC(0.98);
    stats2->SetY1NDC(0.7); stats2->SetY2NDC(0.8);
  }
  stats2->SetTextColor(kRed);

  TLegend* legend;
  if (is_t) {
    legend = new TLegend(0.15, 0.45, 0.33, 0.55);
  } else {
    legend = new TLegend(0.8, 0.55, 0.98, 0.65);
  }

  legend->AddEntry(htemp, "No cuts", "l");
  legend->AddEntry(htemp_cut, "Cuts", "f");
  legend->Draw();
}


void plot_variables () {

  double P_mass = 0.938272;
  double N_mass = 0.9395654;

  TFile *file = TFile::Open("/home/lorena/Thesis/JLAB_EIC/data/0pDVCS_inbending_FTPhotonsCorrected_test.root");

  TTree *tree = (TTree*) file->Get("pDVCS_stripped");
 
  tree->Print();

  std::vector<branches_cuts> branch_names = {
    {"_strip_Q2", 0, 12, 1, 6}, // Q2 > M_proton for DIS
    {"_strip_Xbj", 0, 0.8, 0.1, 0.6}, // Valence region
    {"_t_Nuc", -12, 1, -2.5, 0.1}, // Most of the data
    {"_t_Ph", -12, 1, -2.5, 0.1}, // Same as above
    {"_delta_t", -2, 2, -0.5, 0.5}, // Almost 0 for DVCS
    {"_mm2_eNg", -2, 6, N_mass - 0.1, N_mass + 0.1}, // Expecting neutron mass
    {"_mm2_eNg_N", -2, 2, -0.2, 0.2}, // Expecting nothing
    {"_mm2_eg", -3, 6, P_mass - 0.1, P_mass + 0.1}, // Expecting missing proton
    {"_mm2_eNX_N", -10, 10, -0.5, 0.5}, // Expecting photon
    {"_strip_El_chi2pid", -5, 5, -1, 0.8},
    {"_strip_Ph_chi2pid", -1, 1, -1, 1},
    {"_strip_Nuc_chi2pid", -6, 6, -0.8, 1.5}
  };

  gStyle->SetOptStat(0);
  TCanvas *canvas_all = new TCanvas("canvas 1D Histograms", "canvas 1D Histograms", 1800, 1600);
  canvas_all->Divide(3, 4);

  int canvas_entry = 1;

  for (const auto& branch_name : branch_names) {
    // TCanvas *canvas = new TCanvas(("canvas_" + branch_name.name).c_str(), branch_name.name.c_str(), 800, 600);
    canvas_all->cd(canvas_entry);

    std::string hist_name = "htemp_" + std::to_string(canvas_entry);
    std::string hist_name_cut = "htemp_cut_" + std::to_string(canvas_entry);

    std::string no_cut = branch_name.name + " > " + std::to_string(branch_name.min) + " && " +
                         branch_name.name + " < " + std::to_string(branch_name.max);

    std::string cut = branch_name.name + " > " + std::to_string(branch_name.cut_min) + " && " +
                      branch_name.name + " < " + std::to_string(branch_name.cut_max);

    TH1D* htemp = new TH1D(hist_name.c_str(), ("DVCS" + branch_name.name).c_str(), 40, branch_name.min, branch_name.max);
    TH1D* htemp_cut = new TH1D(hist_name_cut.c_str(), ("DVCS" + branch_name.name).c_str(), 40, branch_name.min, branch_name.max);
    
    tree->Project(hist_name.c_str(), branch_name.name.c_str(), no_cut.c_str());
    tree->Project(hist_name_cut.c_str(), branch_name.name.c_str(), cut.c_str());

    htemp->GetXaxis()->SetTitle(("DVCS" + branch_name.name).c_str());
    htemp->GetYaxis()->SetTitle("Events");
    
    // double mean = htemp->GetMean();
    // double mid_range = (htemp->GetXaxis()->GetXmax() - htemp->GetXaxis()->GetXmin()) / 2;
    // htemp->GetXaxis()->SetRangeUser(mean - mid_range, mean + mid_range);
    
    // htemp->SetName(("hist_" + branch_name.name).c_str());
    htemp->SetLineColor(kBlack);
    htemp_cut->SetFillStyle(3004);
    htemp_cut->SetLineColor(kRed);
    htemp_cut->SetFillColor(kRed - 9);

    htemp->SetStats(true);
    htemp_cut->SetStats(true);

    gStyle->SetOptStat("emr");

    htemp->Draw("HIST");
    htemp_cut->Draw("HIST SAMES");
    
    stats_legend(htemp, htemp_cut, branch_name.name);

    canvas_entry++;

    // std::string output_file = "./plots_variables/" + branch_name.name + "_plot.pdf";
    // canvas->SaveAs(output_file.c_str());
    // delete canvas;
}


  std::string output_file = "./plots_variables/Histograms_Variables.png";
  std::string output_file_pdf = "./plots_variables/Histograms_Variables.pdf";
  canvas_all->SaveAs(output_file.c_str());
  canvas_all->SaveAs(output_file_pdf.c_str());
  delete canvas_all;

  //****//
  // 2D //
  //****//

  // std::vector<branches_cuts_2D> branch_names_2D = {
  //   {"Electron_P_vs_Phi", "_strip_El_Phi", "_strip_El_P", -180, 180, 0, 9},
  //   {"Electron_Theta_vs_Phi", "_strip_El_Phi", "_strip_El_Theta", -180, 180, 0, 50},
    
  //   {"Photon_P_vs_Phi", "_strip_Ph_Phi", "_strip_Ph_P", -180, 180, 0, 11},
  //   {"Photon_Theta_vs_Phi", "_strip_Ph_Phi", "_strip_Ph_Theta", -180, 180, 0, 40},
    
  //   {"Nucleon_P_vs_Phi", "_strip_Nuc_Phi", "_strip_Nuc_P", -180, 180, 0, 8},
  //   {"Nucleon_Theta_vs_Phi", "_strip_Nuc_Phi", "_strip_Nuc_Theta", -180, 180, 0, 100}
  // };
  
  // gStyle->SetOptStat(0);
  
  // for (const auto& branch_name_2D : branch_names_2D) {
  //   TCanvas *canvas = new TCanvas(("canvas_" + branch_name_2D.name).c_str(), branch_name_2D.name.c_str(), 1300, 600);

  //   std::string histDef = ">>htemp(60," +
  //                         std::to_string(branch_name_2D.x_min) + "," +
  //                         std::to_string(branch_name_2D.x_max) + "," + "60," +
  //                         std::to_string(branch_name_2D.y_min) + "," +
  //                         std::to_string(branch_name_2D.y_max) + ")";

  //   std::string cut = branch_name_2D.x_branch + " >= " + std::to_string(branch_name_2D.x_min) + " && " +
  //                     branch_name_2D.x_branch + " <= " + std::to_string(branch_name_2D.x_max) + " && " +
  //                     branch_name_2D.y_branch + " >= " + std::to_string(branch_name_2D.y_min) + " && " +
  //                     branch_name_2D.y_branch + " <= " + std::to_string(branch_name_2D.y_max);
    
  //   tree->Draw((branch_name_2D.y_branch + ":" + branch_name_2D.x_branch + histDef).c_str(), cut.c_str(), "COLZ");
  //   TH2D *htemp = (TH2D*) gDirectory->Get("htemp");
    
  //   htemp->GetXaxis()->SetLimits(branch_name_2D.x_min, branch_name_2D.x_max);
  //   htemp->GetYaxis()->SetLimits(branch_name_2D.y_min, branch_name_2D.y_max);
  //   htemp->SetName(("hist_" + branch_name_2D.name).c_str());
  //   htemp->SetTitle(("DVCS " + branch_name_2D.name).c_str());
  //   htemp->GetXaxis()->SetTitle(branch_name_2D.x_branch.c_str());
  //   htemp->GetYaxis()->SetTitle(branch_name_2D.y_branch.c_str());

  //   std::string output_file = "./plots/" + branch_name_2D.name + "_plot.pdf";
  //   canvas->SaveAs(output_file.c_str());
  //   delete canvas;
  //   delete htemp;
  // };

  delete file;
 
}