// to run: clas12root -l -b -q 'compare_MC_data_background.C'
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TPaveStats.h>

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

struct CutLabels
{
  TString title; // logical/common name
  TString data;  // actual cut label in analysis_data.root
  TString mc;    // actual cut label in analysis_MCsignal.root
  TString bkg;   // actual cut label in analysis_background.root
};

bool should_set_logy(const TString &branch_name)
{
  return (branch_name == "t_Nuc" || branch_name == "t_Ph");
}

void AdjustHistogramRange(TH1D *hist)
{
  if (!hist)
    return;

  int firstBin = -1, lastBin = -1;

  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    if (hist->GetBinContent(i) > 0)
    {
      if (firstBin == -1)
        firstBin = i;
      lastBin = i;
    }
  }

  if (firstBin != -1 && lastBin != -1)
  {
    double xMin = hist->GetXaxis()->GetBinLowEdge(firstBin);
    double xMax = hist->GetXaxis()->GetBinUpEdge(lastBin);
    hist->GetXaxis()->SetRangeUser(xMin, xMax);
  }
}

TH1D *GetHistClone(TFile *file, const TString &name)
{
  if (!file)
    return nullptr;

  TH1D *h = dynamic_cast<TH1D *>(file->Get(name));
  if (!h)
    return nullptr;

  TH1D *hc = dynamic_cast<TH1D *>(h->Clone(Form("%s_clone", name.Data())));
  if (!hc)
    return nullptr;

  hc->SetDirectory(nullptr);
  return hc;
}

TString VarNameForSample(const TString &sample, const TString &logical_var)
{
  if (sample == "background")
    return "_" + logical_var;
  return logical_var;
}

TString BaseHistName(const TString &sample, const TString &logical_var)
{
  TString v = VarNameForSample(sample, logical_var);

  if (sample == "MCsignal")
    return Form("h%s_base_MCsignal", v.Data());
  if (sample == "data")
    return Form("h%s_base_data", v.Data());
  if (sample == "background")
    return Form("h%s_base_background", v.Data());

  return "";
}

TString CutHistName(const TString &sample, const TString &logical_var, const TString &cut_label)
{
  TString v = VarNameForSample(sample, logical_var);

  if (sample == "MCsignal")
    return Form("h%s_%s_MCsignal", v.Data(), cut_label.Data());
  if (sample == "data")
    return Form("h%s_%s_data", v.Data(), cut_label.Data());
  if (sample == "background")
    return Form("h%s_%s_background", v.Data(), cut_label.Data());

  return "";
}

void stats_legend(TH1D *hMC, TH1D *hData, TH1D *hBkg,
                  const TString &branch_name,
                  const std::map<TString, TString> &latex_labels)
{
  if (!hMC || !hData || !hBkg)
    return;

  gPad->cd();
  gPad->SetLogy(0);

  hMC->SetStats(true);
  hData->SetStats(true);
  hBkg->SetStats(true);

  hData->Draw("HIST");
  gPad->Update();

  auto *stats1 = (TPaveStats *)hData->FindObject("stats");

  hMC->Draw("HIST SAMES");
  gPad->Update();

  auto *stats0 = (TPaveStats *)hMC->FindObject("stats");

  hBkg->Draw("HIST SAMES");
  gPad->Update();

  auto *stats2 = (TPaveStats *)hBkg->FindObject("stats");

  hMC->SetFillStyle(0);

  hData->GetXaxis()->SetTitle(Form("DVCS %s", latex_labels.at(branch_name).Data()));
  hData->GetYaxis()->SetTitle("Events");

  double max_total = std::max({hData->GetMaximum(), hMC->GetMaximum(), hBkg->GetMaximum()});

  if (should_set_logy(branch_name))
  {
    gPad->SetLogy();
    hData->SetMinimum(1e-6);
    hMC->SetMinimum(1e-6);
    hBkg->SetMinimum(1e-6);
    hData->SetMaximum(100.0 * max_total);
  }
  else
  {
    hData->SetMinimum(0.0);
    hMC->SetMinimum(0.0);
    hBkg->SetMinimum(0.0);
    hData->SetMaximum(1.4 * max_total);
  }

  gPad->Update();

  // auto *stats0 = (TPaveStats *)hMC->FindObject("stats");
  // auto *stats1 = (TPaveStats *)hData->FindObject("stats");
  // auto *stats2 = (TPaveStats *)hBkg->FindObject("stats");

  if (stats0)
  {
    stats0->SetX1NDC(0.15);
    stats0->SetX2NDC(0.33);
    stats0->SetY1NDC(0.78);
    stats0->SetY2NDC(0.88);
    stats0->SetTextColor(kBlack);
  }

  if (stats1)
  {
    stats1->SetX1NDC(0.36);
    stats1->SetX2NDC(0.54);
    stats1->SetY1NDC(0.78);
    stats1->SetY2NDC(0.88);
    stats1->SetTextColor(kRed);
  }

  if (stats2)
  {
    stats2->SetX1NDC(0.57);
    stats2->SetX2NDC(0.75);
    stats2->SetY1NDC(0.78);
    stats2->SetY2NDC(0.88);
    stats2->SetTextColor(kBlue);
  }

  TLegend *legend = new TLegend(0.78, 0.78, 0.90, 0.90);
  legend->AddEntry(hMC, "MC", "l");
  legend->AddEntry(hData, "Data", "f");
  legend->AddEntry(hBkg, "Bkg", "f");
  legend->Draw();

  gPad->Modified();
  gPad->Update();

  std::cout << "stats MC   = " << stats0 << std::endl;
  std::cout << "stats Data = " << stats1 << std::endl;
  std::cout << "stats Bkg  = " << stats2 << std::endl;
}

void compare_MC_data_background(bool compare_all_cuts = false)
{
  TFile *file_data = TFile::Open("/w/hallb-scshelf2102/clas12/nlbucuru/JLAB_EIC/DVCS/output_root_hists/analysis_data.root");
  TFile *file_background = TFile::Open("./output_root_hists/analysis_background.root");
  TFile *file_MC = TFile::Open("./output_root_hists/analysis_MCsignal.root");

  if (!file_data || !file_background || !file_MC)
  {
    std::cerr << "Could not open one or more input ROOT files.\n";
    return;
  }

  // Use logical/common variable names (no leading underscore).
  std::vector<TString> branch_names = {
      "Pmiss_mag",
      "Pmiss_perp",
      "Pmiss_Nuc_mag",
      "Pmiss_Nuc_perp",
      "mm2_eg",
      "mm2_eNg",
      "mm2_eNg_N",
      "mm2_eNX_N",
      "strip_Q2",
      "strip_Xbj",
      "t_Nuc",
      "t_Ph",
      "delta_t",
      "Phi_Nuc",
      "Phi_Ph",
      "delta_Phi"};

  std::map<TString, TString> latex_labels = {
      {"Pmiss_mag", "|P_{miss}| (GeV)"},
      {"Pmiss_perp", "|P_{miss}^{Perp}| (GeV)"},
      {"Pmiss_Nuc_mag", "|P_{miss} (Nuc)| (GeV)"},
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
      {"delta_Phi", "#Delta#Phi"}};

  // Common cuts only.
  // Data has an extra first cut: cut0bestCand.
  std::vector<CutLabels> cuts = {
      {"theta_gamma_e",
       "cut1theta_gamma_e",
       "cut0theta_gamma_e",
       "cut0_theta_gamma_e"},

      {"chi2pid",
       "cut2chi2pid",
       "cut1chi2pid",
       "cut1_chi2pid"},

      {"delta_t",
       "cut3delta_t",
       "cut2delta_t",
       "cut2_delta_t"},

      {"delta_phi",
       "cut4delta_phi",
       "cut3delta_phi",
       "cut3_delta_phi"},

      {"mm2_eNg_neutron_expected",
       "cut5mm2_eNg_neutron_expected",
       "cut4mm2_eNg_neutron_expected",
       "cut4_mm2_eNg_neutron_expected"},

      {"mm2_eNg_N_nothing_expected",
       "cut6mm2_eNg_N_nothing_expected",
       "cut5mm2_eNg_N_nothing_expected",
       "cut5_mm2_eNg_N_nothing_expected"},

      {"mm2_eNX_N_photon_expected",
       "cut7mm2_eNX_N_photon_expected",
       "cut6mm2_eNX_N_photon_expected",
       "cut6_mm2_eNX_N_photon_expected"},

      {"mm2_eg_proton_expected",
       "cut8mm2_eg_proton_expected",
       "cut7mm2_eg_proton_expected",
       "cut7_mm2_eg_proton_expected"}};

  std::vector<CutLabels> active_cuts;

  if (compare_all_cuts)
  {
    active_cuts = cuts;
  }
  else
  {
    active_cuts = {cuts.back()}; // only the last cut
  }

  gStyle->SetOptTitle(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  int plots_per_canvas = 4;
  int canvas_index = 0;
  TCanvas *canvas_base = nullptr;

  // -------------------------
  // Base histograms
  // -------------------------
  for (size_t i = 0; i < branch_names.size(); ++i)
  {
    if (i % plots_per_canvas == 0)
    {
      if (canvas_base)
      {
        canvas_base->SaveAs(Form("./comparison/base_all_%d.png", canvas_index));
        canvas_base->SaveAs(Form("./comparison/base_all_%d.pdf", canvas_index));
        delete canvas_base;
      }

      canvas_base = new TCanvas(Form("canvas_base_%d", canvas_index),
                                Form("MC-data-Bkg base comparison %d", canvas_index),
                                1200, 900);
      canvas_base->Divide(2, 2);
      ++canvas_index;
    }

    const TString &var = branch_names[i];

    TH1D *h_base_MC = GetHistClone(file_MC, BaseHistName("MCsignal", var));
    TH1D *h_base_data = GetHistClone(file_data, BaseHistName("data", var));
    TH1D *h_base_bkg = GetHistClone(file_background, BaseHistName("background", var));

    if (!h_base_MC || !h_base_data || !h_base_bkg)
    {
      std::cerr << "Skipping base " << var << " due to missing histstats_legend()ograms.\n";
      delete h_base_MC;
      delete h_base_data;
      delete h_base_bkg;
      continue;
    }

    AdjustHistogramRange(h_base_MC);
    AdjustHistogramRange(h_base_data);
    AdjustHistogramRange(h_base_bkg);

    double int_base_MC = h_base_MC->Integral();
    double int_base_data = h_base_data->Integral();
    double int_base_bkg = h_base_bkg->Integral();

    if (int_base_MC > 0)
      h_base_MC->Scale(1.0 / int_base_MC);
    if (int_base_data > 0)
      h_base_data->Scale(1.0 / int_base_data);
    if (int_base_bkg > 0)
      h_base_bkg->Scale(1.0 / int_base_bkg);

    h_base_MC->SetLineColor(kBlack);

    h_base_data->SetLineColor(kRed);
    h_base_data->SetFillColor(kRed - 9);
    h_base_data->SetFillStyle(3004);

    h_base_bkg->SetLineColor(kBlue);
    h_base_bkg->SetFillColor(kBlue - 9);
    h_base_bkg->SetFillStyle(3004);

    gStyle->SetOptStat("emr");
    canvas_base->cd((i % plots_per_canvas) + 1);
    stats_legend(h_base_MC, h_base_data, h_base_bkg, var, latex_labels);

    // delete h_base_MC;
    // delete h_base_data;
    // delete h_base_bkg;
  }

  if (canvas_base)
  {
    canvas_base->SaveAs(Form("./comparison/base_all_%d.png", canvas_index));
    canvas_base->SaveAs(Form("./comparison/base_all_%d.pdf", canvas_index));
    delete canvas_base;
  }

  // -------------------------
  // Cut histograms
  // -------------------------
  for (const auto &cut : active_cuts)
  {
    int canvas_cut_index = 0;
    int plot_index = 0;
    TCanvas *canvas_cut = nullptr;
    bool canvas_has_been_created = false;

    for (size_t i = 0; i < branch_names.size(); ++i)
    {
      const TString &var = branch_names[i];

      bool is_missing_momentum_var = var.Contains("Pmiss");
      bool is_last_cut = (cut.title == "mm2_eNX_N_photon_expected" ||
                          cut.title == "mm2_eg_proton_expected");

      if (is_missing_momentum_var && !is_last_cut)
        continue;

      if (plot_index % plots_per_canvas == 0)
      {
        if (canvas_has_been_created)
        {
          canvas_cut->SaveAs(Form("./comparison/%s_%d.png", cut.title.Data(), canvas_cut_index));
          canvas_cut->SaveAs(Form("./comparison/%s_%d.pdf", cut.title.Data(), canvas_cut_index));
          delete canvas_cut;
        }

        canvas_cut = new TCanvas(Form("canvas_cut_%s_%d", cut.title.Data(), canvas_cut_index),
                                 Form("MC-data-Bkg cut comparison %s %d", cut.title.Data(), canvas_cut_index),
                                 1200, 900);
        canvas_cut->Divide(2, 2);
        ++canvas_cut_index;
        canvas_has_been_created = true;
      }

      TH1D *h_cut_MC = GetHistClone(file_MC,
                                    CutHistName("MCsignal", var, cut.mc));

      TH1D *h_cut_data = GetHistClone(file_data,
                                      CutHistName("data", var, cut.data));

      TH1D *h_cut_bkg = GetHistClone(file_background,
                                     CutHistName("background", var, cut.bkg));

      if (!h_cut_MC || !h_cut_data || !h_cut_bkg)
      {
        std::cerr << "Skipping " << var << " in cut " << cut.title
                  << " due to missing histograms.\n";
        // delete h_cut_MC;
        // delete h_cut_data;
        // delete h_cut_bkg;
        continue;
      }

      AdjustHistogramRange(h_cut_MC);
      AdjustHistogramRange(h_cut_data);
      AdjustHistogramRange(h_cut_bkg);

      double int_cut_MC = h_cut_MC->Integral();
      double int_cut_data = h_cut_data->Integral();
      double int_cut_bkg = h_cut_bkg->Integral();

      if (int_cut_data > 0 && int_cut_MC > 0)
        h_cut_MC->Scale(int_cut_data / int_cut_MC);

      if (int_cut_data > 0 && int_cut_bkg > 0)
        h_cut_bkg->Scale(int_cut_data / int_cut_bkg);

      h_cut_MC->SetLineColor(kBlack);

      h_cut_data->SetLineColor(kRed);
      h_cut_data->SetFillColor(kRed - 9);
      h_cut_data->SetFillStyle(3004);

      h_cut_bkg->SetLineColor(kBlue);
      h_cut_bkg->SetFillColor(kBlue - 9);
      h_cut_bkg->SetFillStyle(3004);

      gStyle->SetOptStat("emr");
      canvas_cut->cd((plot_index % plots_per_canvas) + 1);
      stats_legend(h_cut_MC, h_cut_data, h_cut_bkg, var, latex_labels);

      // delete h_cut_MC;
      // delete h_cut_data;
      // delete h_cut_bkg;

      ++plot_index;
    }

    if (canvas_cut)
    {
      canvas_cut->SaveAs(Form("./comparison/%s_%d.png", cut.title.Data(), canvas_cut_index));
      canvas_cut->SaveAs(Form("./comparison/%s_%d.pdf", cut.title.Data(), canvas_cut_index));
      delete canvas_cut;
    }
  }

  file_MC->Close();
  file_data->Close();
  file_background->Close();
}