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
#include <TSystem.h>

#include <TTreeFormula.h>
#include <TEntryList.h>
#include <TObjArray.h>
#include <limits>

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

    {"strip_Nuc_Phi", "#phi_{n'}"},
    {"strip_Nuc_P", "P_{n'}"},
    {"strip_Nuc_Theta", "#theta_{n'}"},

    {"strip_Spec_Phi", "#phi_{p_{s}}"},
    {"strip_Spec_P", "P_{p_{s}}"},
    {"strip_Spec_Theta", "#theta_{p_{s}}"},

    {"miss_mom_eSg", "|P_{n', miss}| (GeV)"},

    {"strip_Xbjd", "x'_{B,d}"},
    {"t_Nuc", "t_{Nuc} (GeV^{2})"},
    {"t_Ph", "t_{Ph} (GeV^{2})"},
    {"strip_Q2", "Q^{2} (GeV^{2})"},

    {"theta_N_q", "#theta_{n'q} (degrees)"}
};

std::map<TString, TString> plot_titles_2D = {
    {"Electron_P_vs_Phi", "P_{e'} vs #phi_{e'}"},
    {"Electron_Theta_vs_Phi", "#theta_{e'} vs #phi_{e'}"},

    {"Photon_P_vs_Phi", "P_{#gamma} vs #phi_{#gamma}"},
    {"Photon_Theta_vs_Phi", "#theta_{#gamma} vs #phi_{#gamma}"},

    {"Nucleon_P_vs_Phi", "P_{n'} vs #phi_{n'}"},
    {"Nucleon_Theta_vs_Phi", "#theta_{n'} vs #phi_{n'}"},

    {"Spectator_P_vs_Phi", "P_{p_{s}} vs #phi_{p_{s}}"},
    {"Spectator_Theta_vs_Phi", "#theta_{p_{s}} vs #phi_{p_{s}}"},
    {"Spectator_P_vs_Theta", "P_{p_{s}} vs #theta_{p_{s}}"},

    {"miss_mom_eSg_vs_xBj", "|P_{n', miss}| vs x_{B}"},

    {"tNuc_vs_xB", "t_{Nuc} vs x_{B}"},
    {"tPh_vs_xB", "t_{Ph} vs x_{B}"},
    {"Q2_vs_xB", "Q^{2} vs x_{B}"},

    {"Nucleon_P_vs_Spectator_P", "P_{n'} vs P_{p_{s}}"},
    {"Nucleon_Theta_vs_Spectator_Theta", "#theta_{n'} vs #theta_{p_{s}}"},
    {"theta_N_q_vs_Nucleon_Theta", "#theta_{n'q} vs #theta_{n'}"}};

bool should_set_logy(const TString &branch_name)
{
    std::vector<TString> logy_branches = {
        "t_Nuc",
        "t_Ph",
        "mm2_eSg"};

    return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}

TCut kinematic_cut()
{
    // Common DVCS phase-space cuts applied to base and final-cut histograms.
    return TCut("strip_W*strip_W >= 4.0 && strip_Q2 >= 1.0 && strip_Ph_P >= 2.0 && strip_El_P >= 1.0");
}

TString find_branch_for_selection(TTree *tree,
                                  const std::vector<TString> &exact_candidates,
                                  const std::vector<TString> &contains_tokens)
{
    if (!tree)
        return "";

    for (const auto &candidate : exact_candidates)
    {
        if (tree->GetBranch(candidate))
            return candidate;
    }

    TObjArray *branches = tree->GetListOfBranches();
    if (!branches)
        return "";

    for (int i = 0; i < branches->GetEntries(); ++i)
    {
        TObject *obj = branches->At(i);
        if (!obj)
            continue;

        TString branch_name = obj->GetName();
        TString lower_name = branch_name;
        lower_name.ToLower();

        bool matches = true;
        for (auto token : contains_tokens)
        {
            token.ToLower();
            if (!lower_name.Contains(token))
            {
                matches = false;
                break;
            }
        }

        if (matches)
            return branch_name;
    }

    return "";
}

TEntryList *build_highest_energy_entry_list(TTree *tree,
                                            const char *list_name,
                                            const char *proton_p_branch = "strip_Spec_P")
{
    if (!tree)
        return nullptr;

    TString run_branch = find_branch_for_selection(
        tree,
        {"RunNumber"},
        {"run"});

    TString event_branch = find_branch_for_selection(
        tree,
        {"EventNumber"},
        {"event"});

    if (event_branch.IsNull())
    {
        std::cerr << "WARNING: could not find an event-number branch. "
                  << "The highest-energy-per-event preselection will not be applied." << std::endl;
        return nullptr;
    }

    const std::vector<TString> needed = {"strip_El_P", "strip_Ph_P", proton_p_branch};
    for (const auto &bname : needed)
    {
        if (!tree->GetBranch(bname))
        {
            std::cerr << "WARNING: branch needed for highest-energy preselection was not found: "
                      << bname << std::endl;
            return nullptr;
        }
    }

    // For electrons and photons, E is effectively P. For the spectator proton,
    // selecting the highest energy is equivalent to selecting the highest momentum.
    TTreeFormula fRun("fRun_highest_energy", run_branch.IsNull() ? "0" : run_branch.Data(), tree);
    TTreeFormula fEvent("fEvent_highest_energy", event_branch.Data(), tree);
    TTreeFormula fElP("fElP_highest_energy", "strip_El_P", tree);
    TTreeFormula fPhP("fPhP_highest_energy", "strip_Ph_P", tree);
    TTreeFormula fPrP("fPrP_highest_energy", proton_p_branch, tree);

    auto eval_formula = [](TTreeFormula &f) -> double
    {
        f.GetNdata();
        return f.EvalInstance(0);
    };

    auto event_key = [&](double run, double event) -> TString
    {
        if (run_branch.IsNull())
            return Form("%.0f", event);
        return Form("%.0f:%.0f", run, event);
    };

    struct BestEntryInfo
    {
        Long64_t entry = -1;
        double score = -std::numeric_limits<double>::infinity();
        double el_e = -std::numeric_limits<double>::infinity();
        double ph_e = -std::numeric_limits<double>::infinity();
        double pr_e = -std::numeric_limits<double>::infinity();
    };

    const double me = 0.000511;  // GeV
    const double mp = 0.938272;  // GeV
    const double eps = 1e-12;

    std::map<TString, BestEntryInfo> best_by_event;
    Long64_t nentries = tree->GetEntries();

    for (Long64_t ientry = 0; ientry < nentries; ++ientry)
    {
        tree->GetEntry(ientry);

        const double run = eval_formula(fRun);
        const double event = eval_formula(fEvent);
        const double el_p = eval_formula(fElP);
        const double ph_p = eval_formula(fPhP);
        const double pr_p = eval_formula(fPrP);

        if (!std::isfinite(el_p) || !std::isfinite(ph_p) || !std::isfinite(pr_p))
            continue;

        const double el_e = std::sqrt(el_p * el_p + me * me);
        const double ph_e = ph_p;
        const double pr_e = std::sqrt(pr_p * pr_p + mp * mp);

        const double score = el_e + ph_e + pr_e;
        TString key = event_key(run, event);

        auto &best = best_by_event[key];

        bool replace = false;
        if (best.entry < 0 || score > best.score + eps)
        {
            replace = true;
        }
        else if (std::fabs(score - best.score) <= eps)
        {
            // Deterministic tie-breaker: electron, then photon, then proton.
            if (el_e > best.el_e + eps ||
                (std::fabs(el_e - best.el_e) <= eps && ph_e > best.ph_e + eps) ||
                (std::fabs(el_e - best.el_e) <= eps && std::fabs(ph_e - best.ph_e) <= eps && pr_e > best.pr_e + eps))
            {
                replace = true;
            }
        }

        if (replace)
        {
            best.entry = ientry;
            best.score = score;
            best.el_e = el_e;
            best.ph_e = ph_e;
            best.pr_e = pr_e;
        }
    }

    TEntryList *selected_entries = new TEntryList(list_name, list_name);

    for (const auto &item : best_by_event)
    {
        if (item.second.entry >= 0)
            selected_entries->Enter(item.second.entry);
    }

    std::cout << "Highest-energy preselection using event branch '" << event_branch << "'";
    if (!run_branch.IsNull())
        std::cout << " and run branch '" << run_branch << "'";
    std::cout << ": kept " << selected_entries->GetN()
              << " candidate entries out of " << nentries
              << " entries, grouped into " << best_by_event.size()
              << " events." << std::endl;
    
    if (selected_entries->GetN() == 0)
    {
        std::cerr << "WARNING: highest-energy preselection selected zero entries. "
                  << "The preselection will not be applied." << std::endl;
        delete selected_entries;
        return nullptr;
    }

    return selected_entries;
}

std::vector<std::pair<TString, TCut>> generate_cuts(const std::map<TString, TH1D *> &hs_base)
{

    std::vector<std::pair<TString, TString>> cuts_definitions = {
        // {"bestCand", "bestCandidateFlag==1"},
        // {"theta_gamma_e", "theta_gamma_e > 6"}//,
        // {"delta_t", "delta_t >= -1.03182 && delta_t <= 1.27318"},
        {"W2", "strip_W*strip_W >= 4.0"},
        {"Q2", "strip_Q2 >= 1.0"},
        {"gammaP", "strip_Ph_P >= 2.0"},
        {"electronP", "strip_El_P >= 1.0"},
        {"strip_Ph_P", "strip_Ph_P >= 2.0"},
        {"best4DChi2Flag", "best4DChi2Flag == 1"},
        {"Exclusive_4DChi2", "Exclusive_4DChi2 <= 3.0"},
        {"theta_gamma_X", "theta_gamma_X <= 10"},
        {"mm2_eSg_neutron_expected", "mm2_eSg >= -5.02876722725 && mm2_eSg <= 4.84389202852"},
        // {"mm2_eNg_S_nothing_expected", "mm2_eNg_S >= 0 && mm2_eNg_S <= 75"},
        // {"mm2_eg_N_spectator_expected", "mm2_eg_N >= -1 && mm2_eg_N <= 2"},
        // {"mm2_eS_N_photon_expected", "mm2_eS_N >= -1 && mm2_eS_N <= 1"}
        {"delta_phi", "abs(fmod(delta_Phi, 180)) <= 1.5"}
    };

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
        hist->GetXaxis()->SetRangeUser(xMin - 0.1 * (xMax - xMin), xMax + 0.1 * (xMax - xMin));
    }
}

void stats_legend(TH1D *htemp, TH1D *htemp_cut, const TString &branch_name, const std::map<TString, TString> &latex_labels)
{

    gPad->cd();

    htemp->Draw("HIST");
    htemp_cut->Draw("HIST SAMES");
    gPad->Update();

    htemp->SetFillStyle(0);

    TString x_label = branch_name;

    auto label_it = latex_labels.find(branch_name);
    if (label_it != latex_labels.end())
    {
        x_label = label_it->second;
    }
    else
    {
        std::cerr << "WARNING: no latex label found for "
                  << branch_name
                  << ". Using branch name as axis label." << std::endl;
    }

    htemp->GetXaxis()->SetTitle(Form("DVCS %s", x_label.Data()));
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

    TFile *file = TFile::Open("/work/clas12/nlbucuru/PhD_DVCS/outputs_from_hipo/merged_spring2019_nDVCS.root");
    // TFile *file = TFile::Open("/work/clas12/nlbucuru/PhD_DVCS/outputs_from_hipo/merged_fall2019_nDVCS.root");

    if (!file || file->IsZombie())
    {
        std::cerr << "Error: could not open input file" << std::endl;
        return;
    }

    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */
    TString plotDir = "./cuts_data_spring2019_nDVCS";
    // TString plotDir = "./cuts_data_fall2019_10_4_nDVCS";
    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */

    if (gSystem->AccessPathName(plotDir))
    {
        if (gSystem->mkdir(plotDir, true) != 0)
        {
            std::cerr << "Error: could not create plot directory " << plotDir << std::endl;
            return;
        }
    }

    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */
    // TString outputDir = "/w/hallb-scshelf2102/clas12/nlbucuru/JLAB_EIC/DVCS/output_root_hists";
    TString outputDir = "/w/hallb-scshelf2102/clas12/nlbucuru/JLAB_EIC/DVCS/output_root_hists_spring2019_nDVCS";
    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */

    if (gSystem->AccessPathName(outputDir))
    {
        gSystem->mkdir(outputDir, true);
    }

    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */
    TString outputName = outputDir + "/analysis_data_spring2019_nDVCS.root";
    // TString outputName = outputDir + "/analysis_data_fall2019_nDVCS.root";
    /****************************************************** */
    /* To modify to change the name of the output directory */
    /****************************************************** */

    TFile *output_file = new TFile(outputName, "RECREATE");

    if (!output_file || output_file->IsZombie())
    {
        std::cerr << "Error: could not create output file " << outputName << std::endl;
        return;
    }

    TTree *tree = (TTree *)file->Get("nDVCS");

    if (!tree)
    {
        std::cerr << "Error: tree nDVCS was not found in the input file" << std::endl;
        file->ls();
        return;
    }

    auto require_branch = [&](const char *bname) -> bool
    {
        if (!tree->GetBranch(bname))
        {
            std::cerr << "Error: required branch not found in nDVCS tree: " << bname << std::endl;
            return false;
        }
        return true;
    };

    if (!require_branch("mm2_eSg") ||
        !require_branch("mm2_eNg_S") ||
        !require_branch("mm2_eg_N") ||
        !require_branch("mm2_eS_N") ||
        !require_branch("mm2_eS_N") ||
        !require_branch("strip_W") ||
        !require_branch("strip_Q2") ||
        !require_branch("strip_Ph_P") ||
        !require_branch("strip_El_P"))
    {
        tree->GetListOfBranches()->Print();
        return;
    }

    //   tree->Print();

    gROOT->cd();
    TEntryList *highest_energy_entries = build_highest_energy_entry_list(tree,
                                                                      Form("highest_energy_entries_%s", tree->GetName()),
                                                                      "strip_Spec_P");
    if (highest_energy_entries)
        tree->SetEntryList(highest_energy_entries);

    std::vector<std::pair<TString, std::pair<double, double>>> branch_names = {
        {"miss_mom_eSg", {-0.1, 10}}, // |P_{n', miss}| (GeV)
        {"strip_Spec_P", {0, 6}},     // |P_{p_{s}}| (GeV)
        // {"strip_Nuc_P", {0, 10}},     // |P_{n'}| (GeV)
        {"strip_El_P", {0, 10}}, // |P_{e'}| (GeV)
        {"strip_Ph_P", {0, 10}}, // |P_{#gamma}| (GeV)
        // {"p_perp", {-0.1, 3}},

        {"mm2_eSg", {-26, 20}},          // (ElectronBeam + Target_Vec - Spec_Vec - El_Vec - Ph_Vec).M2() Expected mass squared of the missing neutron (should peak at neutron mass squared)
        // {"mm2_eg_N", {0.88, 0.89}},      // (ElectronBeam + Target_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2() Expected mass squared of the missing spectator nucleon (should peak at proton mass squared)
        // {"mm2_eNg_S", {-0, 15e-30}},     // (ElectronBeam + Target_Vec - Nuc_Vec - Spec_Vec - El_Vec - Ph_Vec).M2() Expected zero
        {"mm2_eS_N", {-60e-15, 60e-15}}, // (ElectronBeam + Target_Vec - Spec_Vec - El_Vec - Nuc_Vec).M2() Expected mass squared of the missing photon (should peak at zero)
        {"strip_Nuc_Theta", {0, 180}},
        {"strip_Spec_Theta", {0, 180}},

        {"strip_Xbjd", {0, 2.0}},
        {"strip_Xbj", {0, 2.0}},
        {"strip_Q2", {1, 10}},
        {"strip_W", {0, 5.5}},

        {"delta_t", {-1, 1}},
        {"delta_Phi", {-400, 400}},
        {"t_Nuc", {-10, 2}},
        {"t_Ph", {-10, 2}},

        {"Phi_Nuc", {0, 360}},
        {"Phi_Ph", {0, 360}}
        // {"strip_El_chi2pid", {-5.5, 5.5}},
        // {"strip_Ph_chi2pid", {-0.2, 10100}},
        // {"strip_Spec_chi2pid", {-100, 100}}
    };

    // printf("Number of branches: %lu\n", branch_names.size());
    // for (const auto& branch_name : branch_names) {
    //   std::cout << "Branch name: " << branch_name.first << ", Min: " << branch_name.second.first << ", Max: " << branch_name.second.second << std::endl;
    // }

    std::map<TString, TString> latex_labels = {
        {"miss_mom_eSg", "|P_{n', miss}| (GeV)"},
        {"strip_Spec_P", "P_{p_{s}} (GeV)"},
        // {"p_perp", "|P_{n', miss}^{Perp}| (GeV)"},
        // {"strip_Nuc_P", "P_{n'} (GeV)"},
        {"strip_El_P", "P_{e'} (GeV)"},
        {"strip_Ph_P", "P_{#gamma} (GeV)"},

        {"mm2_eSg", "MM^{2}_{n} eD#rightarrow e'#gamma p_{s}(n^{'}_{miss}) (GeV^{2})"},      // (ElectronBeam + Target_Vec - Spec_Vec - El_Vec - Ph_Vec).M2() Expected mass squared of the missing neutron (should peak at neutron mass squared)
        // {"mm2_eg_N", "MM^{2}_{p} eD#rightarrow e'#gamma n' (p_{s,miss}) (GeV^{2})"},         // (ElectronBeam + Target_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2() Expected mass squared of the missing spectator nucleon (should peak at proton mass squared)
        // {"mm2_eNg_S", "MM^{2}_{X} eD#rightarrow e'#gamma p_{s} n' (X_{miss}) (GeV^{2})"},    // (ElectronBeam + Target_Vec - Nuc_Vec - Spec_Vec - El_Vec - Ph_Vec).M2() Expected zero
        {"mm2_eS_N", "MM^{2}_{#gamma} eD#rightarrow e' n' p_{s} (#gamma_{miss}) (GeV^{2})"}, // (ElectronBeam + Target_Vec - Spec_Vec - El_Vec - Nuc_Vec).M2() Expected mass squared of the missing photon (should peak at zero)
        {"strip_Nuc_Theta", "#theta_{n'} (degrees)"},
        {"strip_Spec_Theta", "#theta_{p_{s}} (degrees)"},

        {"strip_Xbjd", "x'_{B,d}"},
        {"strip_Xbj", "x_{B}"},
        {"strip_Q2", "Q^{2}"},
        {"strip_W", "W (GeV)"},

        {"delta_t", "#Delta t (GeV^{2})"},
        {"delta_Phi", "#Delta #Phi (degrees)"},
        {"t_Nuc", "t_{Nuc} (GeV^{2})"},
        {"t_Ph", "t_{Ph} (GeV^{2})"},

        {"Phi_Nuc", "#Phi_{Nuc} (degrees)"},
        {"Phi_Ph", "#Phi_{#gamma} (degrees)"}
    };

    std::map<TString, TH1D *> hs_base_data;

    for (const auto &[var, range] : branch_names)
    {
        const auto &[min, max] = range;

        TString base_hist_name_data = Form("h%s_base_data", var.Data());
        TH1D *h_base_data = nullptr;

        if (var == "t_Nuc" || var == "t_Ph" || var == "mm2_eg_N" || var == "mm2_eNg_S" || var == "delta_t")
        {
            h_base_data = new TH1D(base_hist_name_data, Form("DVCS%s_data", var.Data()), 1000, min, max);
        }
        else
        {
            h_base_data = new TH1D(base_hist_name_data, Form("DVCS%s_data", var.Data()), 200, min, max);
        }

        // tree->Project(base_hist_name_data, var, "best4DChi2Flag == 1"); //  && Exclusive_4DChi2 <= 3.0
        TCut base_cut_data = kinematic_cut() && TCut("best4DChi2Flag == 1");
        tree->Project(base_hist_name_data, var, base_cut_data); // base = kinematic cuts + best candidate

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

            bool is_last_cut = (label_cut.Contains("mm2_eg_N_spectator_expected") ||
                                label_cut.Contains("mm2_eS_N_photon_expected") ||
                                label_cut.Contains("delta_phi"));

            if ((!is_last_cut) || var.Contains("chi2pid"))
                continue;

            if (plot_index % plots_per_canvas == 0)
            {
                if (canvas_has_been_created)
                {
                    canvas->SaveAs(Form("%s/optimization_%s_data_%d.png", plotDir.Data(), label_cut.Data(), canvas_index));
                    canvas->SaveAs(Form("%s/optimization_%s_data_%d.pdf", plotDir.Data(), label_cut.Data(), canvas_index));
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

            if (var == "t_Nuc" || var == "t_Ph" || var == "mp_eg" || var == "mp_eNg" || var == "mp_eNg_N" || var == "mp_eNX_N" || var == "delta_t")
            {
                h_cut_data = new TH1D(cut_hist_name_data, Form("DVCS%s", var.Data()), 1000, min, max);
            }
            else
            {
                h_cut_data = new TH1D(cut_hist_name_data, Form("DVCS%s", var.Data()), 200, min, max);
            }

            tree->Project(cut_hist_name_data, var, kinematic_cut() && cut);

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
            canvas->SaveAs(Form("%s/optimization_%s_data_%d.png", plotDir.Data(), label_cut.Data(), canvas_index));
            canvas->SaveAs(Form("%s/optimization_%s_data_%d.pdf", plotDir.Data(), label_cut.Data(), canvas_index));
            delete canvas;
        }
    }

    // ******** //
    // 2D DATA  //
    // ******** //

    std::vector<branch_cut_2D> branch_names_2D = {
        // {name, x_branch, y_branch, x_min, x_max, y_min, y_max, use_newtree}
        {"Electron_P_vs_Phi", "strip_El_Phi", "strip_El_P", -180, 180, 0, 9, false},
        {"Electron_Theta_vs_Phi", "strip_El_Phi", "strip_El_Theta", -180, 180, 0, 50, false},

        {"Photon_P_vs_Phi", "strip_Ph_Phi", "strip_Ph_P", -180, 180, 0, 11, false},
        {"Photon_Theta_vs_Phi", "strip_Ph_Phi", "strip_Ph_Theta", -180, 180, 0, 40, false},

        {"Nucleon_P_vs_Phi", "strip_Nuc_Phi", "strip_Nuc_P", -180, 180, 0, 30, false},
        {"Nucleon_Theta_vs_Phi", "strip_Nuc_Phi", "strip_Nuc_Theta", -180, 180, 0, 180, false},
        {"Nucleon_P_vs_Theta", "strip_Nuc_Theta", "strip_Nuc_P", 0, 180, 0, 30, false},

        {"Nucleon_P_vs_Spectator_P", "strip_Spec_P", "strip_Nuc_P", 0, 6.0, 0, 10, false},
        {"Nucleon_Theta_vs_Spectator_Theta", "strip_Spec_Theta", "strip_Nuc_Theta", 0, 180, 0, 180, false},
        {"theta_N_q_vs_Nucleon_Theta", "strip_Nuc_Theta", "theta_N_q", 0, 180, 0, 180, false},

        {"Spectator_P_vs_Phi", "strip_Spec_Phi", "strip_Spec_P", -180, 180, 0, 6.0, false},
        {"Spectator_Theta_vs_Phi", "strip_Spec_Phi", "strip_Spec_Theta", -180, 180, 0, 180, false},
        {"Spectator_P_vs_Theta", "strip_Spec_Theta", "strip_Spec_P", 0, 180, 0, 6.0, false},

        // Example if later you want a Pmiss one from newtree:
        {"miss_mom_eSg_vs_xBj", "strip_Xbjd", "miss_mom_eSg", 0, 1.5, 0, 10, false},

        {"tNuc_vs_xB", "strip_Xbjd", "t_Nuc", 0, 1.5, -8, 1, false},
        {"tPh_vs_xB", "strip_Xbjd", "t_Ph", 0, 1.5, -8, 1, false},
        {"Q2_vs_xB", "strip_Xbjd", "strip_Q2", 0, 1.5, 1, 8, false}};

    gStyle->SetOptStat(0);

    for (const auto &h2def : branch_names_2D)
    {
        // TTree *srcTree = h2def.use_newtree ? newtree : tree;
        TTree *srcTree = tree;

        if (!srcTree->GetBranch(h2def.x_branch) || !srcTree->GetBranch(h2def.y_branch))
        {
            std::cerr << "WARNING: skipping 2D plot " << h2def.name
                      << " because branch " << h2def.x_branch
                      << " or " << h2def.y_branch << " is missing." << std::endl;
            continue;
        }

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

        // TCut final_cut = cuts.back().second;
        // srcTree->Draw(draw_expr, final_cut && TCut(cut_expr), "goff");

        // Full 2D selection = phase-space cuts + the last cumulative cut from generate_cuts().
        TCut final_2D_cut = kinematic_cut();
        if (!cuts.empty())
            final_2D_cut = final_2D_cut && cuts.back().second;

        srcTree->Draw(draw_expr, final_2D_cut && TCut(cut_expr), "goff");

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
        c2->SaveAs(Form("%s/%s_data.pdf", plotDir.Data(), h2def.name.Data()));
        c2->SaveAs(Form("%s/%s_data.png", plotDir.Data(), h2def.name.Data()));

        delete c2;
        delete h2;
    }

    tree->SetEntryList(nullptr);
    delete highest_energy_entries;

    file->Close();
    output_file->Close();
    delete output_file;
}
