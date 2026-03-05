// To run: root -l -b -q 'analysis_background.C'

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

#include <map>
#include <algorithm>
#include <cmath>

#include <TLorentzVector.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TMath.h>

bool should_set_logy(const TString &branch_name)
{
  std::vector<TString> logy_branches = {
      "t_Nuc",
      "t_Ph"};

  return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}

std::vector<std::pair<TString, TCut>> generate_cuts(const std::map<TString, TH1D *> &hs_base)
{
  std::vector<std::pair<TString, TString>> cuts_definitions = {
      {"bestCand", "bestCandidateFlag==1"},
      {"theta_gamma_e", "theta_gamma_e > 6"},
      // {"chi2pid", "strip_El_chi2pid >= -4.56920 && strip_El_chi2pid <= 3.61976 && strip_Nuc_chi2pid >= -195.04711 && strip_Nuc_chi2pid <= 201.30658"},
      {"chi2pid", "strip_El_chi2pid >= -4.37898 && strip_El_chi2pid <= 3.77798 && strip_Nuc_chi2pid >= -5.57697 && strip_Nuc_chi2pid <= 6.33990"},
      // {"delta_t", "delta_t >= -0.46292 && delta_t <= 0.47175"},delta_t >= -1.03182 && delta_t <= 1.27318
      {"delta_t", "delta_t >= -1.03182 && delta_t <= 1.27318"},
      {"delta_phi", "abs(fmod(delta_Phi, 180)) <= 1.5"},
      // {"mm2_eNg_neutron_expected", "mm2_eNg >= -0.37894 && mm2_eNg <= 2.42267"},
      {"mm2_eNg_neutron_expected", "mm2_eNg >= -1.22968 && mm2_eNg <= 3.67237"},
      // {"mm2_eNg_N_nothing_expected", "mm2_eNg_N >= -0.19478 && mm2_eNg_N <= 0.15635"},
      {"mm2_eNg_N_nothing_expected", "mm2_eNg_N >= -0.55184 && mm2_eNg_N <= 0.47227"},
      // {"mm2_eNX_N_photon_expected", "mm2_eNX_N >= -3.95236 && mm2_eNX_N <= 3.74568"},
      {"mm2_eNX_N_photon_expected", "mm2_eNX_N >= -4.09789 && mm2_eNX_N <= 4.17903"},
      // {"mm2_eg_proton_expected", "mm2_eg >= -0.12854 && mm2_eg <= 2.21362"}};
      {"mm2_eg_proton_expected", "mm2_eg >= -1.35581 && mm2_eg <= 3.93882"}};

  std::vector<std::pair<TString, TCut>> cuts;

  TCut all_cuts = "";

  for (size_t i = 0; i < cuts_definitions.size(); ++i)
  {
    const auto &[cut_label, cut] = cuts_definitions[i];

    TCut this_cut = TCut(cut.Data());
    all_cuts += this_cut;

    TString cut_name = Form("cut%lu%s", i, cut_label.Data());
    cuts.emplace_back(cut_name, all_cuts);

    // std::cout << cut_name << " : " << all_cuts << std::endl;
  }

  return cuts;
}

void AdjustHistogramRange(TH1D *hist)
{
  int firstBin = -1, lastBin = -1;

  // Find the first and last non-empty bins
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    if (hist->GetBinContent(i) > 0)
    {
      if (firstBin == -1)
        firstBin = i;
      lastBin = i;
    }
  }

  // Adjust the axis range
  if (firstBin != -1 && lastBin != -1)
  {
    double xMin = hist->GetXaxis()->GetBinLowEdge(firstBin);
    double xMax = hist->GetXaxis()->GetBinUpEdge(lastBin);
    hist->GetXaxis()->SetRangeUser(xMin, xMax);
  }
}

void stats_legend(TH1D *htemp, TH1D *htemp_cut, const TString &branch_name, const std::map<TString, TString> &latex_labels)
{

  gPad->cd();

  htemp->Draw("HIST");
  htemp_cut->Draw("HIST SAMES");
  gPad->Update();

  htemp->SetFillStyle(0);

  htemp->GetXaxis()->SetTitle(Form("DVCS %s", latex_labels.at(branch_name).Data()));
  htemp->GetYaxis()->SetTitle("Events");
  // htemp->SetMinimum(10.0);

  // double max_total = std::max({htemp->GetMaximum(), htemp_cut->GetMaximum()});

  // if (should_set_logy(branch_name))
  // {
  //   gPad->SetLogy();
  //   htemp_cut->SetMaximum(100 * max_total);
  // }
  // else
  // {
  //   htemp_cut->SetMaximum(1.4 * max_total);
  // }

  // gPad->Update();

  auto stats1 = (TPaveStats *)htemp->FindObject("stats");
  auto stats2 = (TPaveStats *)htemp_cut->FindObject("stats");

  if (stats1)
  {
    stats1->SetX1NDC(0.15);
    stats1->SetX2NDC(0.33);
    stats1->SetY1NDC(0.78);
    stats1->SetY2NDC(0.88);

    stats1->SetTextColor(kBlack);
  }

  if (stats2)
  {
    stats2->SetX1NDC(0.15);
    stats2->SetX2NDC(0.33);
    stats2->SetY1NDC(0.66);
    stats2->SetY2NDC(0.76);

    stats2->SetTextColor(kRed);
  }
  
  TLegend *legend = new TLegend(0.36, 0.78, 0.46, 0.88);

  legend->AddEntry(htemp, "No cuts", "l");
  legend->AddEntry(htemp_cut, "Cuts", "f");
  legend->Draw();

  // gPad->ModifiedUpdate();
}

void analysis_background()
{

  TFile *file = TFile::Open("./data/0pDVCS_Pi0dataAsDVCS_10p2.root");
  TFile *output_file = new TFile("./output_root_hists/analysis_background.root", "RECREATE");

  TTree *tree = (TTree *)file->Get("pDVCS");

  //   tree->Print();

  double Pmass = 0.938272;
  double Nmass = 0.9395654;
  double Dmass = 1.8756;
  int RunNumber = 0;
  tree->SetBranchAddress("RunNumber", &RunNumber);

  tree->GetEntry(0);
  double Ebeam = (RunNumber > 10000) ? 10.4 : 10.2;
  // double Ebeam = 10.2;

  // if (RunNumber >= 6420)
  //   Ebeam = 10.2;

  // if (RunNumber > 10000)
  //   Ebeam = 10.4;

  TLorentzVector ElectronBeam;
  TLorentzVector NucTarget_Vec;
  TLorentzVector Target_Vec;
  TLorentzVector Ph_Vec;
  TLorentzVector Nuc_Vec;
  TLorentzVector El_Vec;

  ElectronBeam.SetXYZT(0, 0, Ebeam, Ebeam);
  Target_Vec.SetXYZT(0, 0, 0, Dmass);
  NucTarget_Vec.SetXYZT(0, 0, 0, Pmass);

  std::vector<double> *strip_El_px = nullptr;
  std::vector<double> *strip_El_py = nullptr;
  std::vector<double> *strip_El_pz = nullptr;
  std::vector<double> *strip_El_E = nullptr;

  std::vector<double> *strip_Nuc_px = nullptr;
  std::vector<double> *strip_Nuc_py = nullptr;
  std::vector<double> *strip_Nuc_pz = nullptr;
  std::vector<double> *strip_Nuc_E = nullptr;

  std::vector<double> *strip_Ph_px = nullptr;
  std::vector<double> *strip_Ph_py = nullptr;
  std::vector<double> *strip_Ph_pz = nullptr;
  std::vector<double> *strip_Ph_E = nullptr;

  tree->SetBranchAddress("strip_El_px", &strip_El_px);
  tree->SetBranchAddress("strip_El_py", &strip_El_py);
  tree->SetBranchAddress("strip_El_pz", &strip_El_pz);
  tree->SetBranchAddress("strip_El_E", &strip_El_E);

  tree->SetBranchAddress("strip_Nuc_px", &strip_Nuc_px);
  tree->SetBranchAddress("strip_Nuc_py", &strip_Nuc_py);
  tree->SetBranchAddress("strip_Nuc_pz", &strip_Nuc_pz);
  tree->SetBranchAddress("strip_Nuc_E", &strip_Nuc_E);

  tree->SetBranchAddress("strip_Ph_px", &strip_Ph_px);
  tree->SetBranchAddress("strip_Ph_py", &strip_Ph_py);
  tree->SetBranchAddress("strip_Ph_pz", &strip_Ph_pz);
  tree->SetBranchAddress("strip_Ph_E", &strip_Ph_E);

  std::map<TString, TH1D *> hs_base_MCsignal;

  TLorentzVector Pmiss;
  TLorentzVector Pmiss_Nuc;

  std::vector<double> Pmiss_mag_v, Pmiss_perp_v, Pmiss_Nuc_mag_v, Pmiss_Nuc_perp_v;

  TTree *newtree = tree->CloneTree(0);
  newtree->Branch("Pmiss_mag", &Pmiss_mag_v);
  newtree->Branch("Pmiss_perp", &Pmiss_perp_v);
  newtree->Branch("Pmiss_Nuc_mag", &Pmiss_Nuc_mag_v);
  newtree->Branch("Pmiss_Nuc_perp", &Pmiss_Nuc_perp_v);

  for (Long64_t i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);

    Pmiss_mag_v.clear();
    Pmiss_perp_v.clear();
    Pmiss_Nuc_mag_v.clear();
    Pmiss_Nuc_perp_v.clear();

    size_t n = std::min({strip_El_px->size(), strip_El_py->size(), strip_El_pz->size(), strip_El_E->size(),
                         strip_Ph_px->size(), strip_Ph_py->size(), strip_Ph_pz->size(), strip_Ph_E->size(),
                         strip_Nuc_px->size(), strip_Nuc_py->size(), strip_Nuc_pz->size(), strip_Nuc_E->size()});
    for (size_t j = 0; j < n; j++)
    {
      Ph_Vec.SetPxPyPzE(strip_Ph_px->at(j), strip_Ph_py->at(j), strip_Ph_pz->at(j), strip_Ph_E->at(j));
      Nuc_Vec.SetPxPyPzE(strip_Nuc_px->at(j), strip_Nuc_py->at(j), strip_Nuc_pz->at(j), strip_Nuc_E->at(j));
      El_Vec.SetPxPyPzE(strip_El_px->at(j), strip_El_py->at(j), strip_El_pz->at(j), strip_El_E->at(j));

      Pmiss = ElectronBeam + Target_Vec - El_Vec - Ph_Vec - Nuc_Vec;
      Pmiss_Nuc = ElectronBeam + NucTarget_Vec - El_Vec - Ph_Vec - Nuc_Vec;

      Pmiss_mag_v.push_back(Pmiss.P());
      Pmiss_Nuc_mag_v.push_back(Pmiss_Nuc.P());

      Pmiss_perp_v.push_back(Pmiss.Perp());
      Pmiss_Nuc_perp_v.push_back(Pmiss_Nuc.Perp());
    }
    newtree->Fill();
  }

  std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
      {"Pmiss_mag", {-0.1, 1.5}},
      {"Pmiss_perp", {-0.1, 0.8}},
      {"Pmiss_Nuc_mag", {-0.1, 1.5}},
      // {"_miss_mom_eNg", {0, 1.0}},
      {"Pmiss_Nuc_perp", {-0.1, 0.8}},
      {"mm2_eg", {-2, 5.5}},
      {"mm2_eNg", {-1.5, 5}},
      {"mm2_eNg_N", {-1, 1}},
      {"mm2_eNX_N", {-5, 10}},
      {"strip_Q2", {1, 8}},
      {"strip_Xbj", {0, 0.7}},
      {"t_Nuc", {-14, 1}},
      {"t_Ph", {-12, 1}},
      {"delta_t", {-2, 2}},
      {"Phi_Nuc", {0, 360}},
      {"Phi_Ph", {0, 360}},
      {"delta_Phi", {-4, 3}},
      {"strip_El_chi2pid", {-5.5, 5.5}},
      {"strip_Ph_chi2pid", {-0.2, 10100}},
      {"strip_Nuc_chi2pid", {-6, 6}}};

  // printf("Number of branches: %lu\n", branch_names.size());
  // for (const auto& branch_name : branch_names) {
  //   std::cout << "Branch name: " << branch_name.first << ", Min: " << branch_name.second.first << ", Max: " << branch_name.second.second << std::endl;
  // }

  std::map<TString, TString> latex_labels = {
      {"Pmiss_mag", "|P_{miss}| (GeV)"},
      {"Pmiss_perp", "|P_{miss}^{Perp}| (GeV)"},
      {"Pmiss_Nuc_mag", "|P_{miss} (Nuc)| (GeV)"},
      // {"_miss_mom_eNg", "|P_{miss} (Nuc Mostafa)| (GeV)"},
      {"Pmiss_Nuc_perp", "|P_{miss}^{Perp} (Nuc)| (GeV)"},
      {"mm2_eg", "MM^{2}_{P} e P#rightarrow e'#gamma(P_{miss}) (GeV^{2})"},
      {"mm2_eNg", "MM^{2}_{P} e D#rightarrow e'P'#gamma(N_{miss}) (GeV^{2})"},
      {"mm2_eNg_N", "MM^{2}_{X} e P#rightarrow e'P'#gamma (GeV^{2})"},
      {"mm2_eNX_N", "MM^{2}_{#gamma} e P#rightarrow e'P'(#gamma_{miss}) (GeV^{2})"},
      {"strip_Q2", "Q^{2}"},
      {"strip_Xbj", "x_{B}"},
      {"t_Nuc", "t_{Nuc}"},
      {"t_Ph", "t_{Ph}"},
      {"delta_t", "#Delta t"},
      {"Phi_Nuc", "#Phi_{Nuc}"},
      {"Phi_Ph", "#Phi_{Ph}"},
      {"delta_Phi", "#Delta#Phi"},
      {"strip_El_chi2pid", "#chi^{2}_{pid}^{e}"},
      {"strip_Ph_chi2pid", "#chi^{2}_{pid}^{#gamma}"},
      {"strip_Nuc_chi2pid", "#chi^{2}_{pid}^{N}"}};

  std::map<TString, TH1D *>
      hs_base_background;

  for (const auto &[var, range] : branch_names)
  {
    const auto &[min, max] = range;

    TString base_hist_name_background = Form("h%s_base_background", var.Data());

    TH1D *h_base_background = nullptr;

    if (var == "t_Nuc" || var == "t_Ph")
    {
      h_base_background = new TH1D(base_hist_name_background, Form("DVCS%s_background", var.Data()), 200, min, max);
    }
    else
    {
      h_base_background = new TH1D(base_hist_name_background, Form("DVCS%s_background", var.Data()), 200, min, max);
    }

    if (var == "Pmiss_mag" || var == "Pmiss_Nuc_mag" || var == "Pmiss_perp" || var == "Pmiss_Nuc_perp")
    {
      newtree->Project(base_hist_name_background, var, "");
    }
    else
    {
      tree->Project(base_hist_name_background, var, "");
    }

    // h_base_background->Scale(1.0 / h_base_background->Integral());

    h_base_background->SetMaximum(1.5 * h_base_background->GetMaximum());
    h_base_background->SetMinimum(10.0);

    h_base_background->SetLineColor(kBlack);
    h_base_background->SetStats(true);

    hs_base_background[var] = h_base_background;
    gStyle->SetOptStat("emr");

    output_file->cd();
    hs_base_background[var]->Write();
  }

  auto cuts = generate_cuts(hs_base_background);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  int plots_per_canvas = 4;

  for (const auto &[label_cut, cut] : cuts)
  {

    // TCanvas *canvas = new TCanvas("canvas", label_cut, 1920, 1080);
    // canvas->Divide(4, 4);

    int canvas_index = 0;
    int plot_index = 0;
    TCanvas *canvas = nullptr;

    bool canvas_has_been_created = false;

    for (size_t i = 0; i < branch_names.size(); ++i)
    {

      const auto &[var, range] = branch_names[i];
      const auto &[min, max] = range;

      bool is_missing_momentum_var = var.Contains("Pmiss");

      bool is_last_cut = (label_cut.Contains("mm2_eNX_N_photon_expected") ||
                          label_cut.Contains("mm2_eg_proton_expected"));

      if ((is_missing_momentum_var && !is_last_cut) || var.Contains("chi2pid"))
        continue;

      if (plot_index % plots_per_canvas == 0)
      {
        if (canvas_has_been_created)
        {
          canvas->SaveAs(Form("./cuts_background/optimization_%s_background_%d.png", label_cut.Data(), canvas_index));
          canvas->SaveAs(Form("./cuts_background/optimization_%s_background_%d.pdf", label_cut.Data(), canvas_index));
          delete canvas;
        }

        canvas = new TCanvas(Form("canvas_%s_%d", label_cut.Data(), canvas_index),
                             Form("%s - Part %d", label_cut.Data(), canvas_index),
                             1920, 1080);
        canvas->Divide(2, 2);
        ++canvas_index;
        canvas_has_been_created = true;
      }

      // canvas->cd(i + 1);
      canvas->cd((plot_index % plots_per_canvas) + 1);
      if (should_set_logy(var.Data()))
        gPad->SetLogy();

      TString cut_hist_name_background = Form("h%s_%s_background", var.Data(), label_cut.Data());

      TH1D *h_cut_background = nullptr;

      if (var == "t_Nuc" || var == "t_Ph")
      {
        h_cut_background = new TH1D(cut_hist_name_background, Form("DVCS%s_background", var.Data()), 200, min, max);
      }
      else
      {
        h_cut_background = new TH1D(cut_hist_name_background, Form("DVCS%s_background", var.Data()), 200, min, max);
      }

      if (is_missing_momentum_var)
      {
        newtree->Project(cut_hist_name_background, var, cut);
      }
      else
      {
        tree->Project(cut_hist_name_background, var, cut);
      }

      h_cut_background->SetMinimum(10.0);
      h_cut_background->SetLineColor(kRed);
      h_cut_background->SetFillColor(kRed - 9);
      h_cut_background->SetFillStyle(3004);
      h_cut_background->SetStats(true);
      // h_cut_background->Scale(1.0 / h_cut_background->Integral());

      gStyle->SetOptStat("emr");

      AdjustHistogramRange(hs_base_background[var]);
      AdjustHistogramRange(h_cut_background);

      // gPad->ModifiedUpdate();

      THStack *stack = new THStack(Form("stack%s", var.Data()), Form("DVCS%s", var.Data()));
      stack->Add(hs_base_background[var]);
      stack->Add(h_cut_background);
      stack->Draw("nostack");

      stats_legend(hs_base_background[var], h_cut_background, var, latex_labels);

      // gPad->Modified();

      output_file->cd();
      h_cut_background->Write();
      ++plot_index;
    }

    // // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.png", label_cut.Data());
    // TString filename = Form("./cuts_background/optimization_%s_background.png", label_cut.Data());
    // TString filename_pdf = Form("./cuts_background/optimization_%s_background.pdf", label_cut.Data());
    // canvas->SaveAs(filename);
    // canvas->SaveAs(filename_pdf);
    // delete canvas;

    if (canvas_has_been_created)
    {
      canvas->SaveAs(Form("./cuts_background/optimization_%s_background_%d.png", label_cut.Data(), canvas_index));
      canvas->SaveAs(Form("./cuts_background/optimization_%s_background_%d.pdf", label_cut.Data(), canvas_index));
      delete canvas;
    }
  }

  file->Close();
  output_file->Close();
  delete output_file;
}
