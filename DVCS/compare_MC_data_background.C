#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>

bool should_set_logy(const TString &branch_name)
{
  std::vector<TString> logy_branches = {
      "_t_Nuc",
      "_t_Ph"};

  return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
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

void stats_legend(TH1D *htemp_MC, TH1D *htemp_data, TH1D *htemp_bkg,
                  const TString &branch_name,
                  const std::map<TString, TString> &latex_labels) //,
                                                                  // const std::map<TString, std::pair<double, double>> &x_axis_ranges)
{

  gPad->cd();

  htemp_data->Draw("HIST");
  htemp_MC->Draw("HIST SAMES");
  htemp_bkg->Draw("HIST SAMES");

  htemp_MC->SetFillStyle(0);

  htemp_data->GetXaxis()->SetTitle(Form("DVCS %s", latex_labels.at(branch_name).Data()));
  htemp_data->GetYaxis()->SetTitle("Events");
  double max_total = std::max({htemp_data->GetMaximum(), htemp_MC->GetMaximum(), htemp_bkg->GetMaximum()});

  if (should_set_logy(branch_name.Data()))
  {
    gPad->SetLogy();
    htemp_data->SetMaximum(100 * max_total);
  }
  else
  {
    htemp_data->SetMaximum(1.4 * max_total);
  }

  gPad->ModifiedUpdate();

  TPaveStats *stats0 = (TPaveStats *)htemp_MC->FindObject("stats");
  TPaveStats *stats1 = (TPaveStats *)htemp_data->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)htemp_bkg->FindObject("stats");

  stats0->SetX1NDC(0.15);
  stats0->SetX2NDC(0.33);
  stats0->SetY1NDC(0.78);
  stats0->SetY2NDC(0.88);

  stats1->SetX1NDC(0.36);
  stats1->SetX2NDC(0.54);
  stats1->SetY1NDC(0.78);
  stats1->SetY2NDC(0.88);

  stats2->SetX1NDC(0.57);
  stats2->SetX2NDC(0.75);
  stats2->SetY1NDC(0.78);
  stats2->SetY2NDC(0.88);

  stats0->SetTextColor(kBlack);
  stats1->SetTextColor(kRed);
  stats2->SetTextColor(kBlue);

  // legend = new TLegend(0.15, 0.45, 0.33, 0.55);
  TLegend *legend = new TLegend(0.78, 0.78, 0.88, 0.88);

  legend->AddEntry(htemp_MC, "MC", "l");
  legend->AddEntry(htemp_data, "data", "f");
  legend->AddEntry(htemp_bkg, "Bkg", "f");
  legend->Draw();

  gPad->ModifiedUpdate();
}

void compare_MC_data_background()
{
  TFile *file_data = TFile::Open("./output_root_hists/analysis_data.root");
  TFile *file_background = TFile::Open("./output_root_hists/analysis_background.root");
  TFile *file_MC = TFile::Open("./output_root_hists/analysis_MCsignal.root");

  std::vector<TString> branch_names = {
      "_Pmiss_mag",
      "_Pmiss_perp",
      "_Pmiss_Nuc_mag",
      "_miss_mom_eNg",
      "_Pmiss_Nuc_perp",
      "_mm2_eg",
      "_mm2_eNg",
      "_mm2_eNg_N",
      "_mm2_eNX_N",
      "_strip_Q2",
      "_strip_Xbj",
      "_t_Nuc",
      "_t_Ph",
      "_delta_t",
      "_Phi_Nuc",
      "_Phi_Ph",
      "_delta_Phi"};

  std::map<TString, TString> latex_labels = {
      {"_Pmiss_mag", "|P_{miss}| (GeV)"},
      {"_Pmiss_perp", "|P_{miss}^{Perp}| (GeV)"},
      {"_Pmiss_Nuc_mag", "|P_{miss} (Nuc)| (GeV)"},
      {"_miss_mom_eNg", "|P_{miss} (Nuc Mostafa)| (GeV)"},
      {"_Pmiss_Nuc_perp", "|P_{miss}^{Perp} (Nuc)| (GeV)"},
      {"_mm2_eg", "MM^{2}_{P} e P#rightarrow e'#gamma(P_{miss}) (GeV^{2})"},
      {"_mm2_eNg", "MM^{2}_{P} e D#rightarrow e'P'#gamma(N_{miss}) (GeV^{2})"},
      {"_mm2_eNg_N", "MM^{2}_{X} e P#rightarrow e'P'#gamma (GeV^{2})"},
      {"_mm2_eNX_N", "MM^{2}_{#gamma} e P#rightarrow e'P'(#gamma_{miss}) (GeV^{2})"},
      {"_strip_Q2", "Q^{2}"},
      {"_strip_Xbj", "x_{B}"},
      {"_t_Nuc", "t_{Nuc}"},
      {"_t_Ph", "t_{Ph}"},
      {"_delta_t", "#Delta t"},
      {"_Phi_Nuc", "#Phi_{Nuc}"},
      {"_Phi_Ph", "#Phi_{Ph}"},
      {"_delta_Phi", "#Delta#Phi"}};

  std::vector<TString> cut_labels = {
      "cut0_theta_gamma_e",
      "cut1_chi2pid",
      "cut2_delta_t",
      "cut3_delta_phi",
      "cut4_mm2_eNg_neutron_expected",
      "cut5_mm2_eNg_N_nothing_expected",
      "cut6_mm2_eNX_N_photon_expected",
      "cut7_mm2_eg_proton_expected"};

  gStyle->SetOptTitle(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  // TCanvas *canvas_base = new TCanvas("canvas_base_all", "MC_Data-Bkg base comparison", 1600, 1200);
  // canvas_base->Divide(3, 4);

  int plots_per_canvas = 4;
  int canvas_index = 0;
  TCanvas *canvas_base = nullptr;

  for (size_t i = 0; i < branch_names.size(); ++i)
  {
    if (i % plots_per_canvas == 0)
    {
      if (i > 0)
      {
        canvas_base->SaveAs(Form("./comparison/base_all_%d.png", canvas_index));
        canvas_base->SaveAs(Form("./comparison/base_all_%d.pdf", canvas_index));
        delete canvas_base;
      }
      canvas_base = new TCanvas(Form("canvas_base_%d", canvas_index), Form("MC-data-Bkg base comparison %d", canvas_index), 1200, 900);
      canvas_base->Divide(2, 2);
      ++canvas_index;
    }

    const auto &var = branch_names[i];

    TString base_MC = Form("h%s_base_MCsignal", var.Data());
    TString base_sig = Form("h%s_base_data", var.Data());
    TString base_bkg = Form("h%s_base_background", var.Data());

    TH1D *h_base_MC = (TH1D *)file_MC->Get(base_MC);
    TH1D *h_base_data = (TH1D *)file_data->Get(base_sig);
    TH1D *h_base_bkg = (TH1D *)file_background->Get(base_bkg);

    if (!h_base_MC || !h_base_data || !h_base_bkg)
    {
      std::cerr << "Skipping " << var << " due to missing base histograms." << std::endl;
      continue;
    }

    AdjustHistogramRange(h_base_MC);
    AdjustHistogramRange(h_base_data);
    AdjustHistogramRange(h_base_bkg);

    double int_base_data = h_base_data->Integral();
    double int_base_MC = h_base_MC->Integral();
    double int_base_bkg = h_base_bkg->Integral();

    h_base_MC->Scale(1.0 / int_base_MC);
    h_base_data->Scale(1.0 / int_base_data);
    h_base_bkg->Scale(1.0 / int_base_bkg);

    h_base_MC->SetLineColor(kBlack);

    h_base_data->SetLineColor(kRed);
    h_base_data->SetFillColor(kRed - 9);
    h_base_data->SetFillStyle(3004);

    h_base_bkg->SetLineColor(kBlue);
    h_base_bkg->SetFillColor(kBlue - 9);
    h_base_bkg->SetFillStyle(3004);

    gStyle->SetOptStat("emr");

    // canvas_base->cd(i + 1);
    canvas_base->cd(i % plots_per_canvas + 1);

    // stats_legend(h_base_MC, h_base_data, h_base_bkg, var, latex_labels, x_axis_ranges);
    stats_legend(h_base_MC, h_base_data, h_base_bkg, var, latex_labels);
  }

  if (canvas_base)
  {
    canvas_base->SaveAs(Form("./comparison/base_all_%d.png", canvas_index));
    canvas_base->SaveAs(Form("./comparison/base_all_%d.pdf", canvas_index));
    delete canvas_base;
  }

  // canvas_base->SaveAs("./comparison/AllTogether/base_all.png");
  // canvas_base->SaveAs("./comparison/AllTogether/base_all.pdf");
  // delete canvas_base;

  for (const auto &label_cut : cut_labels)
  {
    // TCanvas *canvas_cut = new TCanvas(Form("canvas_cut_%s", label_cut.Data()), Form("MC-data-Bkg cut comparison %s", label_cut.Data()), 1600, 1200);
    // canvas_cut->Divide(3, 4);

    int canvas_cut_index = 0;
    int plot_index = 0;
    TCanvas *canvas_cut = nullptr;

    bool canvas_has_been_created = false;

    for (size_t i = 0; i < branch_names.size(); i++)
    {
      const auto &var = branch_names[i];

      bool is_missing_momentum_var = var.Contains("_Pmiss");

      bool is_last_cut = (label_cut.Contains("_mm2_eNX_N_photon_expected") ||
                          label_cut.Contains("_mm2_eg_proton_expected"));

      if ((is_missing_momentum_var && !is_last_cut) || var.Contains("_chi2pid"))
        continue;

      if (plot_index % plots_per_canvas == 0)
      {
        if (canvas_has_been_created)
        {
          canvas_cut->SaveAs(Form("./comparison/%s_%d.png", label_cut.Data(), canvas_cut_index));
          canvas_cut->SaveAs(Form("./comparison/%s_%d.pdf", label_cut.Data(), canvas_cut_index));
          delete canvas_cut;
        }
        canvas_cut = new TCanvas(Form("canvas_cut_%s_%d", label_cut.Data(), canvas_cut_index), Form("MC-data-Bkg cut comparison %s %d", label_cut.Data(), canvas_cut_index), 1200, 900);
        canvas_cut->Divide(2, 2);
        canvas_cut_index++;
        canvas_has_been_created = true;
      }

      auto cut_MC = Form("h%s_%s_MCsignal", var.Data(), label_cut.Data());
      auto cut_data = Form("h%s_%s_data", var.Data(), label_cut.Data());
      auto cut_bkg = Form("h%s_%s_background", var.Data(), label_cut.Data());

      TH1D *h_cut_MC = (TH1D *)file_MC->Get(cut_MC);
      TH1D *h_cut_data = (TH1D *)file_data->Get(cut_data);
      TH1D *h_cut_bkg = (TH1D *)file_background->Get(cut_bkg);

      if (!h_cut_MC || !h_cut_data || !h_cut_bkg)
      {
        std::cerr << "Skipping " << var << " in cut " << label_cut << " due to missing histograms." << std::endl;
        continue;
      }

      double int_cut_data = h_cut_data->Integral();
      double int_cut_MC = h_cut_MC->Integral();
      double int_cut_bkg = h_cut_bkg->Integral();

      h_cut_MC->Scale(int_cut_data / int_cut_MC);
      h_cut_bkg->Scale(int_cut_data / int_cut_bkg);

      h_cut_MC->SetLineColor(kBlack);

      h_cut_data->SetLineColor(kRed);
      h_cut_data->SetFillColor(kRed - 9);
      h_cut_data->SetFillStyle(3004);

      h_cut_bkg->SetLineColor(kBlue);
      h_cut_bkg->SetFillColor(kBlue - 9);
      h_cut_bkg->SetFillStyle(3004);

      // canvas_cut->cd(i + 1);
      canvas_cut->cd((plot_index % plots_per_canvas) + 1);

      // stats_legend(h_cut_MC, h_cut_data, h_cut_bkg, var, latex_labels, x_axis_ranges);
      stats_legend(h_cut_MC, h_cut_data, h_cut_bkg, var, latex_labels);

      ++plot_index;
    }

    if (canvas_cut)
    {
      canvas_cut->SaveAs(Form("./comparison/%s_%d.png", label_cut.Data(), canvas_cut_index));
      canvas_cut->SaveAs(Form("./comparison/%s_%d.pdf", label_cut.Data(), canvas_cut_index));
      delete canvas_cut;
    }
    // canvas_cut->SaveAs(Form("./comparison/AllTogether/%s.png", label_cut.Data()));
    // canvas_cut->SaveAs(Form("./comparison/AllTogether/%s.pdf", label_cut.Data()));
    // delete canvas_cut;
  }
  file_MC->Close();
  file_data->Close();
  file_background->Close();
}