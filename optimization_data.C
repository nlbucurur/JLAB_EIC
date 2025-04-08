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

void optimization_data () {

  double P_mass = 0.938272;
  double N_mass = 0.9395654;

  TFile *file = TFile::Open("/home/lorena/Thesis/JLAB_EIC/0pDVCS_inbending_FTPhotonsCorrected_test.root");

  TTree *tree = (TTree*) file->Get("pDVCS_stripped");

  std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
      {"_mm2_eg", {-3, 6}},
      {"_mm2_eNg", {-2, 6}},
      {"_mm2_eNg_N", {-2, 2}},
      {"_mm2_eNX_N", {-10, 10}},
      {"_strip_Q2", {0, 12}},
      {"_strip_Xbj", {0, 0.8}},
      {"_t_Nuc", {-12, 1}},
      {"_t_Ph", {-12, 1}},
      {"_delta_t", {-2, 2}},
      {"_Phi_Nuc", {0, 360}},
      {"_Phi_Ph", {0, 360}},
      {"_delta_Phi", {-360, 360}}
  };

  std::vector<std::pair<TString, TString>> cuts_definitions = {
    {"cut1_theta_gamma_e", "_theta_gamma_e > 6"},
    {"cut2_chi2pid", "_strip_El_chi2pid >= -1 && _strip_El_chi2pid <= 0.8 && "
                     "_strip_Ph_chi2pid >= -1 && _strip_Ph_chi2pid <= 1 && "
                     "_strip_Nuc_chi2pid >= -0.8 && _strip_Nuc_chi2pid <= 1.5"},
    {"cut3_delta_t", "_delta_t >= -0.5 && _delta_t <= 0.5"},
    {"cut4_delta_phi", "_delta_Phi >= -1.5 && _delta_Phi <= 1.5"},
    {"cut5_mm2_eNg_neutron_expected", "_mm2_eNg >= -1.23003 && _mm2_eNg <= 3.10917"},
    {"cut6_mm2_eNg_N_nothing_expected", "_mm2_eNg_N >= -0.5 && _mm2_eNg_N <= 0.5"},
    {"cut7_mm2_eNX_N_photon_expected", "_mm2_eNX_N >= -0.5 && _mm2_eNX_N <= 0.5"},
    {"cut8_mm2_eg_proton_expected", "_mm2_eg >= -0.4106 && _mm2_eg <= 2.9626"}
  };

  std::vector<std::pair<TString, TCut>> cuts;

  TCut all_cuts = "";
  for (const auto &def : cuts_definitions) {
    all_cuts += TCut(def.second);
    cuts.emplace_back(def.first, all_cuts);
  };

  gStyle->SetOptStat(0);

  for (const auto& [label_cut, cut] : cuts) {

    TCanvas *canvas = new TCanvas("canvas", label_cut, 2600, 2200);
    canvas->Divide(4, 3);

    for (size_t i = 0; i < branch_names.size(); ++i) {
      const auto& [var, range] = branch_names[i];
      const auto& [min, max] = range;

      canvas->cd(i + 1);

      TString base_hist_name = Form("h_%s_base_%s", var.Data(), label_cut.Data());
      TString cut_hist_name = Form("h_%s_cut_%s", var.Data(), label_cut.Data());

      TH1D* h_base = new TH1D(base_hist_name, Form("DVCS%s", var.Data()), 60, min, max);
      TH1D* h_cut = new TH1D(cut_hist_name, Form("DVCS%s", var.Data()), 60, min, max);

      tree->Project(base_hist_name, var, "");
      tree->Project(cut_hist_name, var, cut);

      h_base->SetLineColor(kBlack);
      h_cut->SetLineColor(kRed);
      h_cut->SetFillColor(kRed - 9);
      h_cut->SetFillStyle(3004);

      h_base->SetStats(true);
      h_cut->SetStats(true);

      gStyle->SetOptStat("emr");

      h_base->Draw("HIST");
      gPad->Update();

      TPaveStats* stats1 = (TPaveStats*)h_base->FindObject("stats");
      if (stats1) {
        stats1->SetX1NDC(0.60); stats1->SetX2NDC(0.83);
        stats1->SetY1NDC(0.70); stats1->SetY2NDC(0.80);
        stats1->SetTextColor(kBlack);
      }

      h_cut->Draw("HIST SAMES");
      gPad->Update();

      TPaveStats* stats2 = (TPaveStats*)h_cut->FindObject("stats");
      if (stats2) {
        stats2->SetX1NDC(0.60); stats2->SetX2NDC(0.83);
        stats2->SetY1NDC(0.58); stats2->SetY2NDC(0.68);
        stats2->SetTextColor(kRed);
      }

      auto legend = new TLegend(0.2, 0.7, 0.4, 0.8);
      legend->AddEntry(h_base, "No Cuts", "l");
      legend->AddEntry(h_cut, "With Cuts", "f");
      legend->Draw();
    }

    // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.pdf", label_cut.Data());
    TString filename = Form("./cuts/optimization_%s.pdf", label_cut.Data());
    canvas->SaveAs(filename);
    delete canvas;
  }

  file->Close();
}
