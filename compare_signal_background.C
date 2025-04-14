#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>

bool should_set_logy (const TString& branch_name) {
    std::vector<TString> logy_branches = {
      "_t_Nuc",
      "_t_Ph"
    };
    
    return std::find(logy_branches.begin(), logy_branches.end(), branch_name) != logy_branches.end();
}


bool should_move_stats (const TString& branch_name) {
    std::vector<TString> move_stats = {
      "_t_Nuc",
      "_Phi_Nuc",
      "_Phi_Ph",
      "_strip_Ph_chi2pid"
    };
    
    return std::find(move_stats.begin(), move_stats.end(), branch_name) != move_stats.end();
}

void stats_legend (TH1D* htemp_sig, TH1D* htemp_bkg, const TString& branch_name) {

    gPad->cd();
  
    htemp_sig->Draw("HIST"); 
    htemp_bkg->Draw("HIST SAMES");
    
    htemp_sig->GetXaxis()->SetTitle(Form("DVCS%s", branch_name.Data()));
    htemp_sig->GetYaxis()->SetTitle("Events");
    htemp_sig->SetMinimum(10.0);
    gPad->Update();
    // htemp_sig->SetMaximum(1.2 * htemp_sig->GetMaximum());

    // gPad->Modified();
  
    TPaveStats* stats1 = (TPaveStats*)htemp_sig->FindObject("stats");
    TPaveStats* stats2 = (TPaveStats*)htemp_bkg->FindObject("stats");
  
    bool move_stats = should_move_stats(branch_name);
    
    if (move_stats) {
      stats1->SetX1NDC(0.15); stats1->SetX2NDC(0.33);
      stats1->SetY1NDC(0.78); stats1->SetY2NDC(0.88);
  
      stats2->SetX1NDC(0.15); stats2->SetX2NDC(0.33);
      stats2->SetY1NDC(0.66); stats2->SetY2NDC(0.76);
    } else {
      stats1->SetX1NDC(0.80); stats1->SetX2NDC(0.98);
      stats1->SetY1NDC(0.85); stats1->SetY2NDC(0.95);
  
      stats2->SetX1NDC(0.80); stats2->SetX2NDC(0.98);
      stats2->SetY1NDC(0.7); stats2->SetY2NDC(0.8);
    }
    stats1->SetTextColor(kBlack);
    stats2->SetTextColor(kBlue);
  
    TLegend* legend;
    if (move_stats) {
      // legend = new TLegend(0.15, 0.45, 0.33, 0.55);
      legend = new TLegend(0.36, 0.78, 0.59, 0.88);
    } else {
      legend = new TLegend(0.8, 0.55, 0.98, 0.65);
    }
  
    legend->AddEntry(htemp_sig, "Signal", "l");
    legend->AddEntry(htemp_bkg, "Bkg", "f");
    legend->Draw();
  }


void compare_signal_background () {
    TFile *file_signal = TFile::Open("./output_files/analysis_signal.root");
    TFile *file_background = TFile::Open("./output_files/analysis_background.root");

    std::vector<TString> branch_names = {
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
        "_delta_Phi",
        "_strip_El_chi2pid",
        "_strip_Ph_chi2pid",
        "_strip_Nuc_chi2pid"
    };

    std::vector<TString> cut_labels = {
        "cut0_theta_gamma_e",
        "cut1_chi2pid",
        "cut2_delta_t",
        "cut3_delta_phi",
        "cut4_mm2_eNg_neutron_expected",
        "cut5_mm2_eNg_N_nothing_expected",
        "cut6_mm2_eNX_N_photon_expected",
        "cut7_mm2_eg_proton_expected"
    };

    TCanvas *canvas_base = new TCanvas("canvas_base_all", "Signal-Backgroun bsae comparison", 1600, 1200);
    canvas_base->Divide(4, 4);

    // gStyle->SetOptStat(0);
    // gStyle->SetPadGridX(true);
    // gStyle->SetPadGridY(true);

    
    for (size_t i = 0; i < branch_names.size(); ++i) {
        const auto& var = branch_names[i];
        
        
        TString base_sig = Form("h%s_base_signal", var.Data());
        TString base_bkg = Form("h%s_base_background", var.Data());
        
        TH1D* h_base_signal = (TH1D*) file_signal->Get(base_sig);
        TH1D* h_base_bkg = (TH1D*) file_background->Get(base_bkg);

        h_base_signal->SetLineColor(kBlack);
        h_base_bkg->SetLineColor(kBlue);
        h_base_bkg->SetFillColor(kBlue - 9);
        h_base_bkg->SetFillStyle(3004);
        
        gStyle->SetOptStat("emr");
        
        canvas_base->cd(i + 1);
        if (should_set_logy(var.Data())) gPad->SetLogy();
        h_base_signal->SetTitle(Form("Signal vs Bkg base for %s", var.Data()));
        
        stats_legend(h_base_signal, h_base_bkg, var);
    }

    canvas_base->SaveAs("./comparison/base_all.png");
    canvas_base->SaveAs("./comparison/base_all.pdf");
    delete canvas_base;

    for (const auto& label_cut : cut_labels) {
        TCanvas* canvas_cut = new TCanvas(Form("canvas_cut_%s", label_cut.Data()), Form("Signal-Bkg cut comparison %s", label_cut.Data()), 1600, 1200);
        canvas_cut->Divide(4, 4);

        for (size_t i = 0; i < branch_names.size(); i++) {
            const auto& var = branch_names[i];

            auto cut_signal = Form("h%s_%s_signal", var.Data(), label_cut.Data());
            auto cut_bkg = Form("h%s_%s_background", var.Data(), label_cut.Data());

            TH1D* h_cut_signal = (TH1D*) file_signal->Get(cut_signal);
            TH1D* h_cut_bkg = (TH1D*) file_background->Get(cut_bkg);
            
            h_cut_signal->SetLineColor(kRed);
            h_cut_bkg->SetLineColor(kBlue);
            h_cut_bkg->SetFillColor(kBlue - 9);
            h_cut_bkg->SetFillStyle(3004);

            canvas_cut->cd(i + 1);
            h_cut_signal->SetTitle(Form("Signal vs Bkg cut %s", var.Data()));
            
            stats_legend(h_cut_signal, h_cut_bkg, var);
        }
        canvas_cut->SaveAs(Form("./comparison/cut_%s.png", label_cut.Data()));
        canvas_cut->SaveAs(Form("./comparison/cut_%s.pdf", label_cut.Data()));
        delete canvas_cut;
    }
    file_signal->Close();
    file_background->Close();
}