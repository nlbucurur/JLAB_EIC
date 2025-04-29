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
#include <TLorentzVector.h>
#include <utility>

bool should_set_logy(const TString &branch_name)
{
    std::vector<TString> logy_branches = {
        "_t_Nuc",
        "_t_Ph"};

    return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}

std::pair<TString, TString> auto_cut(const TString &var, TH1D *hist, const TString &cut_label)
{
    double mean = hist->GetMean();
    double sigma = hist->GetStdDev();

    double min_cut = mean - 3.0 * sigma;
    double max_cut = mean + 3.0 * sigma;

    TString cut = Form("%s >= %.5f && %s <= %.5f", var.Data(), min_cut, var.Data(), max_cut);

    return {cut_label, cut};
}

std::vector<std::pair<TString, TCut>> generate_cuts(const std::map<TString, TH1D *> &hs_base)
{
    std::vector<std::pair<TString, TString>> cuts_definitions = {
        {"_theta_gamma_e", "_theta_gamma_e > 6"},
        {"_chi2pid", ""},
        {"_delta_t", ""},
        {"_delta_phi", "abs(fmod(_delta_Phi, 180)) <= 1.5"},
        {"_mm2_eNg_neutron_expected", ""},
        {"_mm2_eNg_N_nothing_expected", ""},
        {"_mm2_eNX_N_photon_expected", ""},
        {"_mm2_eg_proton_expected", ""}};

    // for (const auto& cut_def : cuts_definitions) {
    //   std::cout << cut_def.first << ": " << cut_def.second << std::endl;
    // }

    for (auto &[cut_label, cut] : cuts_definitions)
    {
        if (cut == "")
        {
            if (cut_label == "_chi2pid")
            {
                TString var_El = "_strip_El_chi2pid";
                TString var_Nuc = "_strip_Nuc_chi2pid";

                TString cut_El = auto_cut(var_El, hs_base.at(var_El), cut_label).second;
                TString cut_Nuc = auto_cut(var_Nuc, hs_base.at(var_Nuc), cut_label).second;

                cut = cut_El + " && " + cut_Nuc;

                // std::cout << "Auto cut" << cut_label << " : " << cut << std::endl;
            }
            else
            {
                TString var;
                if (cut_label == "_mm2_eNg_neutron_expected")
                    var = "_mm2_eNg";
                else if (cut_label == "_mm2_eNg_N_nothing_expected")
                    var = "_mm2_eNg_N";
                else if (cut_label == "_mm2_eNX_N_photon_expected")
                    var = "_mm2_eNX_N";
                else if (cut_label == "_mm2_eg_proton_expected")
                    var = "_mm2_eg";
                else if (cut_label == "_delta_t")
                    var = "_delta_t";
                else
                    continue;

                cut = auto_cut(var, hs_base.at(var), cut_label).second;

                // std::cout << "Auto cut" << cut_label << " : " << cut << std::endl;
            }
        }
    }

    std::vector<std::pair<TString, TCut>> cuts;

    TCut all_cuts = "";

    for (size_t i = 0; i < cuts_definitions.size(); ++i)
    {
        const auto &[cut_label, cut] = cuts_definitions[i];

        TCut this_cut = TCut(cut.Data());
        all_cuts += this_cut;

        TString cut_name = Form("cut%lu%s", i, cut_label.Data());
        cuts.emplace_back(cut_name, all_cuts);

        std::cout << cut_name << " : " << all_cuts << std::endl;
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

    // printf("h_base_MCsignal %s: %f\n", branch_name.Data(), htemp->GetMinimum());
    // printf("h_cut_MCsignal %s: %f\n", branch_name.Data(), htemp_cut->GetMinimum());

    // double max_total = std::max({htemp->GetMaximum(), htemp_cut->GetMaximum()});

    // if (should_set_logy(branch_name))
    //         {
    //             gPad->SetLogy();
    //             htemp_cut->SetMaximum(100 * max_total);
    //         }
    //         else
    //         {
    //             htemp_cut->SetMaximum(1.4 * max_total);
    //         }

    // gPad->Update();

    TPaveStats *stats1 = (TPaveStats *)htemp->FindObject("stats");
    TPaveStats *stats2 = (TPaveStats *)htemp_cut->FindObject("stats");

    stats1->SetX1NDC(0.15);
    stats1->SetX2NDC(0.33);
    stats1->SetY1NDC(0.78);
    stats1->SetY2NDC(0.88);

    stats2->SetX1NDC(0.15);
    stats2->SetX2NDC(0.33);
    stats2->SetY1NDC(0.66);
    stats2->SetY2NDC(0.76);

    stats1->SetTextColor(kBlack);
    stats2->SetTextColor(kRed);

    TLegend *legend = new TLegend(0.36, 0.78, 0.46, 0.88);

    legend->AddEntry(htemp, "No cuts", "l");
    legend->AddEntry(htemp_cut, "Cuts", "f");
    legend->Draw();

    // gPad->ModifiedUpdate();
}

void analysis_MCsignal()
{
    TFile *file = TFile::Open("./data/1pDVCS_simulation.root");
    TFile *output_file = new TFile("./output_root_hists/analysis_MCsignal.root", "RECREATE");

    TTree *tree = (TTree *)file->Get("pDVCS_stripped");

    //   tree->Print();

    double Pmass = 0.938272;
    double Nmass = 0.9395654;
    double Dmass = 1.8756;
    double Ebeam = 10.6;

    TLorentzVector ElectronBeam;
    TLorentzVector NucTarget_Vec;
    TLorentzVector Target_Vec;
    TLorentzVector Ph_Vec;
    TLorentzVector Nuc_Vec;
    TLorentzVector El_Vec;

    ElectronBeam.SetXYZT(0, 0, Ebeam, Ebeam);
    Target_Vec.SetXYZT(0, 0, 0, Dmass);
    NucTarget_Vec.SetXYZT(0, 0, 0, Pmass);

    double _strip_El_px, _strip_El_py, _strip_El_pz, _strip_El_E;
    double _strip_Nuc_px, _strip_Nuc_py, _strip_Nuc_pz, _strip_Nuc_E;
    double _strip_Ph_px, _strip_Ph_py, _strip_Ph_pz, _strip_Ph_E;

    tree->SetBranchAddress("_strip_El_px", &_strip_El_px);
    tree->SetBranchAddress("_strip_El_py", &_strip_El_py);
    tree->SetBranchAddress("_strip_El_pz", &_strip_El_pz);
    tree->SetBranchAddress("_strip_El_E", &_strip_El_E);
    tree->SetBranchAddress("_strip_Nuc_px", &_strip_Nuc_px);
    tree->SetBranchAddress("_strip_Nuc_py", &_strip_Nuc_py);
    tree->SetBranchAddress("_strip_Nuc_pz", &_strip_Nuc_pz);
    tree->SetBranchAddress("_strip_Nuc_E", &_strip_Nuc_E);
    tree->SetBranchAddress("_strip_Ph_px", &_strip_Ph_px);
    tree->SetBranchAddress("_strip_Ph_py", &_strip_Ph_py);
    tree->SetBranchAddress("_strip_Ph_pz", &_strip_Ph_pz);
    tree->SetBranchAddress("_strip_Ph_E", &_strip_Ph_E);

    std::map<TString, TH1D *> hs_base_MCsignal;

    TLorentzVector Pmiss;
    TLorentzVector Pmiss_Nuc;

    double Pmiss_mag;
    double Pmiss_Nuc_mag;
    double Pmiss_perpendicular;
    double Pmiss_Nuc_perpendicular;

    TTree *newtree = tree->CloneTree(0);
    TBranch *branch_Pmiss_mag = newtree->Branch("_Pmiss_mag", &Pmiss_mag, "Pmiss_mag/D");
    TBranch *branch_Pmiss_Nuc_mag = newtree->Branch("_Pmiss_Nuc_mag", &Pmiss_Nuc_mag, "Pmiss_Nuc_mag/D");
    TBranch *branch_Pmiss_perp = newtree->Branch("_Pmiss_perp", &Pmiss_perpendicular, "Pmiss_perp/D");
    TBranch *branch_Pmiss_Nuc_perp = newtree->Branch("_Pmiss_Nuc_perp", &Pmiss_Nuc_perpendicular, "Pmiss_Nuc_perp/D");

    for (int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);

        Ph_Vec.SetXYZT(_strip_Ph_px, _strip_Ph_py, _strip_Ph_pz, _strip_Ph_E);
        Nuc_Vec.SetXYZT(_strip_Nuc_px, _strip_Nuc_py, _strip_Nuc_pz, _strip_Nuc_E);
        El_Vec.SetXYZT(_strip_El_px, _strip_El_py, _strip_El_pz, _strip_El_E);

        Pmiss = ElectronBeam + Target_Vec - El_Vec - Ph_Vec - Nuc_Vec;
        Pmiss_Nuc = ElectronBeam + NucTarget_Vec - El_Vec - Ph_Vec - Nuc_Vec;

        Pmiss_mag = Pmiss.P();
        Pmiss_Nuc_mag = Pmiss_Nuc.P();

        Pmiss_perpendicular = Pmiss.Perp();
        Pmiss_Nuc_perpendicular = Pmiss_Nuc.Perp();

        newtree->Fill();
    }

    std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
        {"_Pmiss_mag", {-0.1, 1.5}},
        {"_Pmiss_perp", {-0.1, 0.8}},
        {"_Pmiss_Nuc_mag", {-0.1, 1.5}},
        {"_Pmiss_Nuc_perp", {-0.1, 0.8}},
        {"_mm2_eg", {-2, 5.5}},
        {"_mm2_eNg", {-1.5, 5}},
        {"_mm2_eNg_N", {-1, 1}},
        {"_mm2_eNX_N", {-5, 10}},
        {"_strip_Q2", {1, 8}},
        {"_strip_Xbj", {0, 0.7}},
        {"_t_Nuc", {-3, 0}},
        {"_t_Ph", {-2, 1}},
        {"_delta_t", {-2, 2}},
        {"_Phi_Nuc", {0, 360}},
        {"_Phi_Ph", {0, 360}},
        {"_delta_Phi", {-4, 3}},
        {"_strip_El_chi2pid", {-5.5, 5.5}},
        {"_strip_Ph_chi2pid", {-0.2, 10100}},
        {"_strip_Nuc_chi2pid", {-6, 6}}};

    // printf("Number of branches: %lu\n", branch_names.size());
    // for (const auto& branch_name : branch_names) {
    //   std::cout << "Branch name: " << branch_name.first << ", Min: " << branch_name.second.first << ", Max: " << branch_name.second.second << std::endl;
    // }

    std::map<TString, TString> latex_labels = {
        {"_Pmiss_mag", "|P_{miss}| (GeV)"},
        {"_Pmiss_perp", "|P_{miss}^{Perp}| (GeV)"},
        {"_Pmiss_Nuc_mag", "|P_{miss} (Nuc)| (GeV)"},
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
        {"_delta_Phi", "#Delta#Phi"},
        {"_strip_El_chi2pid", "#chi^{2}_{pid}^{e}"},
        {"_strip_Ph_chi2pid", "#chi^{2}_{pid}^{#gamma}"},
        {"_strip_Nuc_chi2pid", "#chi^{2}_{pid}^{N}"}};

    for (const auto &[var, range] : branch_names)
    {
        const auto &[min, max] = range;

        TString base_hist_name_MCsignal = Form("h%s_base_MCsignal", var.Data());
        TH1D *h_base_MCsignal = nullptr;

        h_base_MCsignal = new TH1D(base_hist_name_MCsignal, Form("DVCS%s_MCsignal", var.Data()), 100, min, max);

        if (var == "_Pmiss_mag" || var == "_Pmiss_Nuc_mag" || var == "_Pmiss_perp" || var == "_Pmiss_Nuc_perp")
        {
            newtree->Project(base_hist_name_MCsignal, var, "");
        }
        else
        {
            tree->Project(base_hist_name_MCsignal, var, "");
        }

        h_base_MCsignal->SetMaximum(1.5 * h_base_MCsignal->GetMaximum());
        h_base_MCsignal->SetMinimum(10.0);

        // printf("h_base_MCsignal %s: %f\n", var.Data(), h_base_MCsignal->GetMinimum());

        h_base_MCsignal->SetLineColor(kBlack);
        h_base_MCsignal->SetStats(true);
        hs_base_MCsignal[var] = h_base_MCsignal;
        gStyle->SetOptStat("emr");

        output_file->cd();
        hs_base_MCsignal[var]->Write();
    }

    auto cuts = generate_cuts(hs_base_MCsignal);

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

            bool is_missing_momentum_var = var.Contains("_Pmiss");

            bool is_last_cut = (label_cut.Contains("_mm2_eNX_N_photon_expected") ||
                                label_cut.Contains("_mm2_eg_proton_expected"));

            if ((is_missing_momentum_var && !is_last_cut) || var.Contains("_chi2pid"))
                continue;

            if (plot_index % plots_per_canvas == 0)
            {
                if (canvas_has_been_created)
                {
                    canvas->SaveAs(Form("./cutsMC/optimization_%s_MCsignal_%d.png", label_cut.Data(), canvas_index));
                    canvas->SaveAs(Form("./cutsMC/optimization_%s_MCsignal_%d.pdf", label_cut.Data(), canvas_index));
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

            TString cut_hist_name_MCsignal = Form("h%s_%s_MCsignal", var.Data(), label_cut.Data());
            TH1D *h_cut_MCsignal = nullptr;

            h_cut_MCsignal = new TH1D(cut_hist_name_MCsignal, Form("DVCS%s", var.Data()), 100, min, max);

            if (is_missing_momentum_var)
            {
                newtree->Project(cut_hist_name_MCsignal, var, cut);
            }
            else
            {
                tree->Project(cut_hist_name_MCsignal, var, cut);
            }

            h_cut_MCsignal->SetMinimum(10.0);
            // printf("h_cut_MCsignal %s: %f\n", var.Data(), h_cut_MCsignal->GetMinimum());

            h_cut_MCsignal->SetLineColor(kRed);
            h_cut_MCsignal->SetFillColor(kRed - 9);
            h_cut_MCsignal->SetFillStyle(3004);
            h_cut_MCsignal->SetStats(true);

            gStyle->SetOptStat("emr");

            AdjustHistogramRange(hs_base_MCsignal[var]);
            AdjustHistogramRange(h_cut_MCsignal);

            THStack *stack = new THStack(Form("stack%s_MCsignal", var.Data()), Form("DVCS%s_MCsignal", var.Data()));
            stack->Add(hs_base_MCsignal[var]);
            stack->Add(h_cut_MCsignal);
            stack->Draw("nostack");

            stats_legend(hs_base_MCsignal[var], h_cut_MCsignal, var, latex_labels);
            // gPad->Modified();

            output_file->cd();
            h_cut_MCsignal->Write();
            ++plot_index;
        }

        // // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.png", label_cut.Data());
        // TString filename = Form("./cutsMC/AllTogether/optimization_%s_MCsignal.png", label_cut.Data());
        // TString filename_pdf = Form("./cutsMC/AllTogether/optimization_%s_MCsignal.pdf", label_cut.Data());
        // canvas->SaveAs(filename);
        // canvas->SaveAs(filename_pdf);
        // delete canvas;
        if (canvas_has_been_created)
        {
            canvas->SaveAs(Form("./cutsMC/optimization_%s_MCsignal_%d.png", label_cut.Data(), canvas_index));
            canvas->SaveAs(Form("./cutsMC/optimization_%s_MCsignal_%d.pdf", label_cut.Data(), canvas_index));
            delete canvas;
        }
    }

    file->Close();
    output_file->Close();
    delete output_file;
}
