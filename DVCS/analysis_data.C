// To run: root -l -b -q 'analysis_data.C'

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

#include <TH2D.h>

struct branch_cut_2D
{
    TString name;
    TString x_branch;
    TString y_branch;
    double x_min, x_max;
    double y_min, y_max;
    bool use_newtree = false; // true only if one axis comes from Pmiss_* branches
};

std::map<TString, TString> axis_labels_2D = {
    {"strip_El_Phi", "#phi_{e'}"},
    {"strip_El_P", "P_{e'}"},
    {"strip_El_Theta", "#theta_{e'}"},

    {"strip_Ph_Phi", "#phi_{#gamma}"},
    {"strip_Ph_P", "P_{#gamma}"},
    {"strip_Ph_Theta", "#theta_{#gamma}"},

    {"strip_Nuc_Phi", "#phi_{p'}"},
    {"strip_Nuc_P", "P_{p'}"},
    {"strip_Nuc_Theta", "#theta_{p'}"},

    {"Pmiss_mag", "|P_{miss}| (GeV)"},

    {"strip_Xbj", "x_{B}"},
    {"strip_Q2", "Q^{2}"},
    {"t_Nuc", "t_{Nuc}"},
    {"t_Ph", "t_{Ph}"}};

std::map<TString, TString> plot_titles_2D = {
    {"Electron_P_vs_Phi", "P_{e'} vs #phi_{e'}"},
    {"Electron_Theta_vs_Phi", "#theta_{e'} vs #phi_{e'}"},

    {"Photon_P_vs_Phi", "P_{#gamma} vs #phi_{#gamma}"},
    {"Photon_Theta_vs_Phi", "#theta_{#gamma} vs #phi_{#gamma}"},

    {"Nucleon_P_vs_Phi", "P_{p'} vs #phi_{p'}"},
    {"Nucleon_Theta_vs_Phi", "#theta_{p'} vs #phi_{p'}"},

    {"Pmiss_vs_Bj", "Pmiss vs x_{B}"},

    {"tNuc_vs_xB", "t_{Nuc} vs x_{B}"},
    {"tPh_vs_xB", "t_{Ph} vs x_{B}"},
    {"Q2_vs_xB", "Q^{2} vs x_{B}"}
    };

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
    // htemp_cut->SetMaximum(100 * max_total);
    // }
    // else
    // {
    // htemp_cut->SetMaximum(1.4 * max_total);
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

void analysis_data()
{

    // TFile *file = TFile::Open("/w/hallb-scshelf2102/clas12/nlbucuru/PhD_DVCS/stripped_data_spring2019_pDVCS_1.root");
    // TFile *file = TFile::Open("/work/clas12/nlbucuru/PhD_DVCS/outputs/stripped_data_spring2019_pDVCS_merged.root");
    TFile *file = TFile::Open("/work/clas12/nlbucuru/PhD_DVCS/stripped_fall2019_pDVCS_sidisdvcs_011093.root");
    TFile *output_file = new TFile("/w/hallb-scshelf2102/clas12/nlbucuru/JLAB_EIC/DVCS/output_root_hists/analysis_data_hipo.root", "RECREATE");

    TTree *tree = (TTree *)file->Get("pDVCS");

    //   tree->Print();

    double Pmass = 0.938272;
    double Nmass = 0.9395654;
    double Dmass = 1.8756;
    int RunNumber = 0;
    tree->SetBranchAddress("RunNumber", &RunNumber);

    tree->GetEntry(0);
    double Ebeam = (RunNumber > 10000) ? 10.4 : 10.2;

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

    // std::map<TString, TH1D *> hs_base_MCsignal;

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
        {"strip_Xbj", {0, 1.5}},
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

    std::map<TString, TH1D *> hs_base_data;

    for (const auto &[var, range] : branch_names)
    {
        const auto &[min, max] = range;

        TString base_hist_name_data = Form("h%s_base_data", var.Data());
        TH1D *h_base_data = nullptr;

        if (var == "t_Nuc" || var == "t_Ph")
        {
            h_base_data = new TH1D(base_hist_name_data, Form("DVCS%s_data", var.Data()), 200, min, max);
        }
        else
        {
            h_base_data = new TH1D(base_hist_name_data, Form("DVCS%s_data", var.Data()), 200, min, max);
        }

        if (var == "Pmiss_mag" || var == "Pmiss_Nuc_mag" || var == "Pmiss_perp" || var == "Pmiss_Nuc_perp")
        {
            newtree->Project(base_hist_name_data, var, "");
        }
        else
        {
            tree->Project(base_hist_name_data, var, "");
        }

        h_base_data->SetMaximum(1.5 * h_base_data->GetMaximum());
        h_base_data->SetMinimum(10.0);

        h_base_data->SetLineColor(kBlack);
        h_base_data->SetStats(true);
        hs_base_data[var] = h_base_data;
        gStyle->SetOptStat("emr");

        output_file->cd();
        hs_base_data[var]->Write();
    }

    auto cuts = generate_cuts(hs_base_data);
    bool plot_all_cuts = false; // true = plots after every cumulative cut and flase = plots only the last one, which in your case is

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    int plots_per_canvas = 4;

    size_t first_cut = plot_all_cuts ? 0 : cuts.size() - 1;

    for (size_t icut = first_cut; icut < cuts.size(); ++icut)
    {
        const auto &[label_cut, cut] = cuts[icut];

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
                    canvas->SaveAs(Form("./cuts_data/optimization_%s_data_%d.png", label_cut.Data(), canvas_index));
                    canvas->SaveAs(Form("./cuts_data/optimization_%s_data_%d.pdf", label_cut.Data(), canvas_index));
                    delete canvas;
                }

                canvas = new TCanvas(Form("canvas_%s_%d", label_cut.Data(), canvas_index),
                                     Form("%s - Part %d", label_cut.Data(), canvas_index),
                                     1920, 1080);
                canvas->Divide(2, 2);
                ++canvas_index;
                canvas_has_been_created = true;
            }

            canvas->cd((plot_index % plots_per_canvas) + 1);
            // canvas->cd(i + 1);
            if (should_set_logy(var.Data()))
                gPad->SetLogy();

            TString cut_hist_name_data = Form("h%s_%s_data", var.Data(), label_cut.Data());

            TH1D *h_cut_data = nullptr;

            if (var == "t_Nuc" || var == "t_Ph" || var == "mp_eg" || var == "mp_eNg" || var == "mp_eNg_N" || var == "mp_eNX_N")
            {
                h_cut_data = new TH1D(cut_hist_name_data, Form("DVCS%s", var.Data()), 200, min, max);
            }
            else
            {
                h_cut_data = new TH1D(cut_hist_name_data, Form("DVCS%s", var.Data()), 200, min, max);
            }

            if (is_missing_momentum_var)
            {
                newtree->Project(cut_hist_name_data, var, cut);
            }
            else
            {
                tree->Project(cut_hist_name_data, var, cut);
            }

            h_cut_data->SetMinimum(10.0);
            h_cut_data->SetLineColor(kRed);
            h_cut_data->SetFillColor(kRed - 9);
            h_cut_data->SetFillStyle(3004);
            h_cut_data->SetStats(true);

            gStyle->SetOptStat("emr");

            AdjustHistogramRange(hs_base_data[var]);
            AdjustHistogramRange(h_cut_data);

            THStack *stack = new THStack(Form("stack%s_data", var.Data()), Form("DVCS%s_data", var.Data()));
            stack->Add(hs_base_data[var]);
            stack->Add(h_cut_data);
            stack->Draw("nostack");

            stats_legend(hs_base_data[var], h_cut_data, var, latex_labels);
            // gPad->Modified();

            output_file->cd();
            h_cut_data->Write();
            ++plot_index;
        }

        // // TString filename = Form("./cuts_no_chi2pid/optimization_%s_no_chi2pid.png", label_cut.Data());
        // TString filename = Form("./cuts_data/AllTogether/optimization_%s_data.png", label_cut.Data());
        // TString filename_pdf = Form("./cuts_data/AllTogether/optimization_%s_data.pdf", label_cut.Data());
        // canvas->SaveAs(filename);
        // canvas->SaveAs(filename_pdf);
        // delete canvas;

        if (canvas)
        {
            canvas->SaveAs(Form("./cuts_data/optimization_%s_data_%d.png", label_cut.Data(), canvas_index));
            canvas->SaveAs(Form("./cuts_data/optimization_%s_data_%d.pdf", label_cut.Data(), canvas_index));
            delete canvas;
        }
    }

    // ******** //
    // 2D DATA  //
    // ******** //

    std::vector<branch_cut_2D> branch_names_2D = {
        {"Electron_P_vs_Phi", "strip_El_Phi", "strip_El_P", -180, 180, 0, 9, false},
        {"Electron_Theta_vs_Phi", "strip_El_Phi", "strip_El_Theta", -180, 180, 0, 50, false},

        {"Photon_P_vs_Phi", "strip_Ph_Phi", "strip_Ph_P", -180, 180, 0, 11, false},
        {"Photon_Theta_vs_Phi", "strip_Ph_Phi", "strip_Ph_Theta", -180, 180, 0, 40, false},

        {"Nucleon_P_vs_Phi", "strip_Nuc_Phi", "strip_Nuc_P", -180, 180, 0, 8, false},
        {"Nucleon_Theta_vs_Phi", "strip_Nuc_Phi", "strip_Nuc_Theta", -180, 180, 0, 100, false},

        // Example if later you want a Pmiss one from newtree:
        {"Pmiss_vs_Bj", "strip_Xbj", "Pmiss_mag", 0, 1.5, 0, 1.5, true},

        {"tNuc_vs_xB", "strip_Xbj", "t_Nuc", 0, 1.5, -8, 1, false},
        {"tPh_vs_xB", "strip_Xbj", "t_Ph", 0, 1.5, -8, 1, false},
        {"Q2_vs_xB", "strip_Xbj", "strip_Q2", 0, 1.5, 1, 8, false}};

    gStyle->SetOptStat(0);

    for (const auto &h2def : branch_names_2D)
    {
        TTree *srcTree = h2def.use_newtree ? newtree : tree;

        TString hist_name = Form("h2_%s_data", h2def.name.Data());

        gROOT->cd();
        if (gDirectory->Get(hist_name))
            gDirectory->Delete(hist_name + ";*");

        TString draw_expr = Form("%s:%s>>%s(60,%f,%f,60,%f,%f)",
                                 h2def.y_branch.Data(), h2def.x_branch.Data(),
                                 hist_name.Data(),
                                 h2def.x_min, h2def.x_max,
                                 h2def.y_min, h2def.y_max);

        TString cut_expr = Form("%s>=%f && %s<=%f && %s>=%f && %s<=%f",
                                h2def.x_branch.Data(), h2def.x_min,
                                h2def.x_branch.Data(), h2def.x_max,
                                h2def.y_branch.Data(), h2def.y_min,
                                h2def.y_branch.Data(), h2def.y_max);

        // srcTree->Draw(draw_expr, cut_expr, "goff");
        // h2->Draw("COLZ");

        TCut final_cut = cuts.back().second;
        srcTree->Draw(draw_expr, final_cut && TCut(cut_expr), "goff");

        TH2D *h2 = (TH2D *)gDirectory->Get(hist_name);
        if (!h2)
        {
            std::cerr << "Could not build " << hist_name << std::endl;
            continue;
        }

        h2->SetDirectory(nullptr);
        h2->SetTitle(Form("DVCS %s data", plot_titles_2D[h2def.name].Data()));
        h2->GetXaxis()->SetTitle(axis_labels_2D[h2def.x_branch].Data());
        h2->GetYaxis()->SetTitle(axis_labels_2D[h2def.y_branch].Data());
        // h2->SetTitle(plot_titles_2D[h2def.name].Data());

        // Write to analysis_data.root
        output_file->cd();
        // h2->Write();
        h2->Write("", TObject::kOverwrite);

        gROOT->cd();

        // Save image immediately
        TCanvas *c2 = new TCanvas(Form("c2_%s", h2def.name.Data()),
                                  h2def.name, 1300, 600);
        h2->Draw("COLZ");
        c2->SaveAs(Form("./cuts_data/%s_data.pdf", h2def.name.Data()));
        c2->SaveAs(Form("./cuts_data/%s_data.png", h2def.name.Data()));

        delete c2;
        delete h2;
    }

    file->Close();
    output_file->Close();
    delete output_file;
}
