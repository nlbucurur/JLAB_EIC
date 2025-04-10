#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TClass.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>
#include <vector>
#include <utility>

bool should_set_logy (const std::string& branch_name) {
  std::vector<std::string> logy_branches = {
    "_t_Nuc",
    "_t_Ph",
    "_delta_t",
    "_Phi_Nuc",
    "_Phi_Ph",
    "_delta_Phi"
  };
  
  return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}

bool should_move_stats (const std::string& branch_name) {
  std::vector<std::string> move_stats = {
    "_t_Nuc",
    "_t_Ph",
    "_Phi_Nuc",
    "_Phi_Ph"
  };
  
  return std::find(move_stats.begin(), move_stats.end(), branch_name) != move_stats.end();
}

void stats_legend (TH1D* htemp, TH1D* htemp_cut, const std::string& branch_name) {
  gPad->Update();

  TPaveStats* stats1 = (TPaveStats*)htemp->FindObject("stats");
  TPaveStats* stats2 = (TPaveStats*)htemp_cut->FindObject("stats");

  bool move_stats = should_move_stats(branch_name);
  
  if (move_stats) {
    stats1->SetX1NDC(0.15); stats1->SetX2NDC(0.33);
    stats1->SetY1NDC(0.78); stats1->SetY2NDC(0.88);
  } else {
    stats1->SetX1NDC(0.80); stats1->SetX2NDC(0.98);
    stats1->SetY1NDC(0.85); stats1->SetY2NDC(0.95);
  }
  stats1->SetTextColor(kBlack);
  
  if (move_stats) {
    stats2->SetX1NDC(0.15); stats2->SetX2NDC(0.33);
    stats2->SetY1NDC(0.66); stats2->SetY2NDC(0.76);
  } else {
    stats2->SetX1NDC(0.80); stats2->SetX2NDC(0.98);
    stats2->SetY1NDC(0.7); stats2->SetY2NDC(0.8);
  }
  stats2->SetTextColor(kRed);

  TLegend* legend;
  if (move_stats) {
    // legend = new TLegend(0.15, 0.45, 0.33, 0.55);
    legend = new TLegend(0.36, 0.78, 0.59, 0.88);
  } else {
    legend = new TLegend(0.8, 0.55, 0.98, 0.65);
  }

  legend->AddEntry(htemp, "No cuts", "l");
  legend->AddEntry(htemp_cut, "Cuts", "f");
  legend->Draw();
}

void optimization_data () {

  double P_mass = 0.938272;
  double N_mass = 0.9395654;

  TFile *file = TFile::Open("/home/lorena/Thesis/JLAB_EIC/data/0pDVCS_inbending_FTPhotonsCorrected_test.root");

  TTree *tree = (TTree*) file->Get("pDVCS_stripped");

  //   tree->Print();

  std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
      {"_mm2_eg", {-3, 6}},
      {"_mm2_eNg", {-2, 6}},
      {"_mm2_eNg_N", {-2, 2}},
      {"_mm2_eNX_N", {-10, 10}},
      {"_strip_Q2", {0, 12}},
      {"_strip_Xbj", {0, 0.8}},
      {"_t_Nuc", {-12, 1}},
      {"_t_Ph", {-12, 10200}},
      {"_delta_t", {-2, 2}},
      {"_Phi_Nuc", {0, 360}},
      {"_Phi_Ph", {0, 360}},
      {"_delta_Phi", {-360, 360}},
      {"_strip_El_chi2pid", {-5, 5}},
      {"_strip_Ph_chi2pid", {-1, 1}},
      {"_strip_Nuc_chi2pid", {-6, 6}}
  };

  // printf("Number of branches: %lu\n", branch_names.size());
  // for (const auto& branch_name : branch_names) {
  //   std::cout << "Branch name: " << branch_name.first << ", Min: " << branch_name.second.first << ", Max: " << branch_name.second.second << std::endl;
  // }

  std::vector<std::pair<TString, TString>> cuts_definitions = {
    {"_theta_gamma_e", "_theta_gamma_e > 6"},
    {"_chi2pid", "_strip_El_chi2pid >= -4.9668 && _strip_El_chi2pid <= 3.8292 && "
                     "_strip_Ph_chi2pid >= -12 && _strip_Ph_chi2pid <= 10200 && "
                     "_strip_Nuc_chi2pid >= -5.246 && _strip_Nuc_chi2pid <= 6.376"},
    {"_delta_t", "_delta_t >= -1.62602 && _delta_t <= 1.61098"},
    {"_delta_phi", "abs(_delta_Phi) % 180 <= 1.5"},
    {"_mm2_eNg_neutron_expected", "_mm2_eNg >= -1.23003 && _mm2_eNg <= 3.10917"},
    {"_mm2_eNg_N_nothing_expected", "_mm2_eNg_N >= -0.5 && _mm2_eNg_N <= 0.5"},
    {"_mm2_eNX_N_photon_expected", "_mm2_eNX_N >= -0.5 && _mm2_eNX_N <= 0.5"},
    {"_mm2_eg_proton_expected", "_mm2_eg >= -0.4106 && _mm2_eg <= 2.9626"}
  };

  // for (const auto& cut_def : cuts_definitions) {
  //   std::cout << cut_def.first << ": " << cut_def.second << std::endl;
  // }

  std::vector<std::pair<TString, TCut>> cuts;

  TCut all_cuts = "";
  int cut_index = 1;
  for (const auto &def : cuts_definitions) {
    all_cuts += TCut(def.second);
    TString cut_name = Form("cut%i%s", cut_index, def.first.Data());
    cuts.emplace_back(cut_name, all_cuts);
    cut_index++;
    // std::cout << def.first << ": " << all_cuts << std::endl;
  };

  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  for (const auto& [label_cut, cut] : cuts) {

    TCanvas *canvas = new TCanvas("canvas", label_cut, 1920, 1080);
    canvas->Divide(4, 4);

    for (size_t i = 0; i < branch_names.size(); ++i) {
      const auto& [var, range] = branch_names[i];
      const auto& [min, max] = range;
      
      canvas->cd(i + 1);
      if (should_set_logy(var.Data())) gPad->SetLogy();

      TString base_hist_name = Form("h_%s_base_%s", var.Data(), label_cut.Data());
      TString cut_hist_name = Form("h_%s_cut_%s", var.Data(), label_cut.Data());

      TH1D* h_base = new TH1D(base_hist_name, Form("DVCS%s", var.Data()), 60, min, max);
      TH1D* h_cut = new TH1D(cut_hist_name, Form("DVCS%s", var.Data()), 60, min, max);

      tree->Project(base_hist_name, var, "");
      tree->Project(cut_hist_name, var, cut);

      h_base->GetXaxis()->SetTitle(Form("DVCS%s", var.Data()));
      h_base->GetYaxis()->SetTitle("Events");

      h_base->SetLineColor(kBlack);
      h_cut->SetLineColor(kRed);
      h_cut->SetFillColor(kRed - 9);
      h_cut->SetFillStyle(3004);

      h_base->SetStats(true);
      h_cut->SetStats(true);

      gStyle->SetOptStat("emr");

      h_base->Draw("HIST");
      h_cut->Draw("HIST SAMES");

      stats_legend(h_base, h_cut, var.Data());
    }

    // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.png", label_cut.Data());
    TString filename = Form("./cuts/optimization_%s.png", label_cut.Data());
    canvas->SaveAs(filename);
    delete canvas;
  }

  file->Close();
}
