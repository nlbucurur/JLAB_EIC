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

bool should_set_logy(const TString &branch_name)
{
    std::vector<TString> logy_branches = {
        "_t_Nuc",
        "_t_Ph"};

    return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}

bool should_move_stats(const TString &branch_name)
{
    std::vector<TString> move_stats = {
        "_t_Nuc",
        "_Phi_Nuc",
        "_Phi_Ph",
        "_strip_Ph_chi2pid"};

    return std::find(move_stats.begin(), move_stats.end(), branch_name) != move_stats.end();
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

        // std::cout << cut_name << " : " << all_cuts << std::endl;
    }

    return cuts;
}

void stats_legend(TH1D *htemp, TH1D *htemp_cut, const TString &branch_name)
{

    gPad->cd();

    htemp->Draw("HIST");
    htemp_cut->Draw("HIST SAMES");
    gPad->Update();

    htemp->GetXaxis()->SetTitle(Form("DVCS%s", branch_name.Data()));
    htemp->GetYaxis()->SetTitle("Events");
    htemp->SetMinimum(10.0);

    TPaveStats *stats1 = (TPaveStats *)htemp->FindObject("stats");
    TPaveStats *stats2 = (TPaveStats *)htemp_cut->FindObject("stats");

    bool move_stats = should_move_stats(branch_name);

    if (move_stats)
    {
        stats1->SetX1NDC(0.15);
        stats1->SetX2NDC(0.33);
        stats1->SetY1NDC(0.78);
        stats1->SetY2NDC(0.88);

        stats2->SetX1NDC(0.15);
        stats2->SetX2NDC(0.33);
        stats2->SetY1NDC(0.66);
        stats2->SetY2NDC(0.76);
    }
    else
    {
        stats1->SetX1NDC(0.80);
        stats1->SetX2NDC(0.98);
        stats1->SetY1NDC(0.85);
        stats1->SetY2NDC(0.95);

        stats2->SetX1NDC(0.80);
        stats2->SetX2NDC(0.98);
        stats2->SetY1NDC(0.7);
        stats2->SetY2NDC(0.8);
    }
    stats1->SetTextColor(kBlack);
    stats2->SetTextColor(KRed);

    TLegend *legend;
    if (move_stats)
    {
        // legend = new TLegend(0.15, 0.45, 0.33, 0.55);
        legend = new TLegend(0.36, 0.78, 0.59, 0.88);
    }
    else
    {
        legend = new TLegend(0.8, 0.55, 0.98, 0.65);
    }

    legend->AddEntry(htemp, "No cuts", "l");
    legend->AddEntry(htemp_cut, "Cuts", "f");
    legend->Draw();
}

void analysis_MCsignal()
{

    double P_mass = 0.938272;
    double N_mass = 0.9395654;

    TFile *file = TFile::Open("./data/1pDVCS_simulation.root");
    TFile *output_file = new TFile("./output_files/analysis_MCsignal.root", "RECREATE");

    TTree *tree = (TTree *)file->Get("pDVCS_stripped");

    //   tree->Print();

    std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
        {"_mm2_eg", {-1, 3}},
        {"_mm2_eNg", {-1, 3.3}},
        {"_mm2_eNg_N", {-0.5, 0.5}},
        {"_mm2_eNX_N", {-10, 6}},
        {"_strip_Q2", {0, 8}},
        {"_strip_Xbj", {0, 0.7}},
        {"_t_Nuc", {-2, 0}},
        {"_t_Ph", {-1.6, 0.2}},
        {"_delta_t", {-1, 1}},
        {"_Phi_Nuc", {0, 360}},
        {"_Phi_Ph", {0, 360}},
        {"_delta_Phi", {-3, 3}},
        {"_strip_El_chi2pid", {-5.5, 5.5}},
        {"_strip_Ph_chi2pid", {-0.2, 10100}},
        {"_strip_Nuc_chi2pid", {-0.2, 10100}}};

    // printf("Number of branches: %lu\n", branch_names.size());
    // for (const auto& branch_name : branch_names) {
    //   std::cout << "Branch name: " << branch_name.first << ", Min: " << branch_name.second.first << ", Max: " << branch_name.second.second << std::endl;
    // }

    std::map<TString, TH1D *> hs_base_MCsignal;

    for (const auto &[var, range] : branch_names)
    {
        const auto &[min, max] = range;

        TString base_hist_name_MCsignal = Form("h%s_base_MCsignal", var.Data());
        TH1D *h_base_MCsignal = new TH1D(base_hist_name_MCsignal, Form("DVCS%s_MCsignal", var.Data()), 60, min, max);

        tree->Project(base_hist_name_MCsignal, var, "");
        h_base_MCsignal->SetMaximum(1.2 * h_base_MCsignal->GetMaximum());

        h_base_MCsignal->SetLineColor(kBlack);
        h_base_MCsignal->SetStats(true);
        hs_base_MCsignal[var] = h_base_MCsignal;
        gStyle->SetOptStat("emr");

        output_file->cd();
        hs_base_MCsignal[var]->Write();
    }

    auto cuts = generate_cuts(hs_base_MCsignal);

    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    for (const auto &[label_cut, cut] : cuts)
    {

        TCanvas *canvas = new TCanvas("canvas", label_cut, 1920, 1080);
        canvas->Divide(4, 4);

        for (size_t i = 0; i < branch_names.size(); ++i)
        {
            const auto &[var, range] = branch_names[i];
            const auto &[min, max] = range;

            canvas->cd(i + 1);
            if (should_set_logy(var.Data()))
                gPad->SetLogy();

            TString cut_hist_name_MCsignal = Form("h%s_%s_MCsignal", var.Data(), label_cut.Data());

            TH1D *h_cut_MCsignal = new TH1D(cut_hist_name_MCsignal, Form("DVCS%s", var.Data()), 60, min, max);

            tree->Project(cut_hist_name_MCsignal, var, cut);

            h_cut_MCsignal->SetMinimum(10.0);
            h_cut_MCsignal->SetLineColor(KRed);
            h_cut_MCsignal->SetFillColor(kRed - 9);
            h_cut_MCsignal->SetFillStyle(3004);
            h_cut_MCsignal->SetStats(true);

            gStyle->SetOptStat("emr");

            THStack *stack = new THStack(Form("stack%s_MCsignal", var.Data()), Form("DVCS%s_MCsignal", var.Data()));
            stack->Add(hs_base_MCsignal[var]);
            stack->Add(h_cut_MCsignal);
            stack->Draw("nostack");

            stats_legend(hs_base_MCsignal[var], h_cut_MCsignal, var.Data());
            gPad->Modified();

            output_file->cd();
            h_cut_MCsignal->Write();
        }

        // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.png", label_cut.Data());
        TString filename = Form("./cutsMC/optimization_%s_MCsignal.png", label_cut.Data());
        TString filename_pdf = Form("./cutsMC/optimization_%s_MCsignal.pdf", label_cut.Data());
        canvas->SaveAs(filename);
        canvas->SaveAs(filename_pdf);
        delete canvas;
    }

    file->Close();
    output_file->Close();
    delete output_file;
}
