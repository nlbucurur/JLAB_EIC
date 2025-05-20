#define BSA_data_fit_cxx
#include "BSA_data_fit.h"

void superimposehist(TH1F *h1, TH1F *h2)
{

    TLatex tex; //
    tex.SetTextAlign(13);
    tex.SetTextColor(17);
    tex.SetTextSize(0.08);
    // tex.SetTextAngle(26.15998);
    tex.SetLineWidth(1);

    if (h1->GetMaximum() > h2->GetMaximum())
    {
        h1->Draw();
        tex.DrawLatex((h1->GetBinLowEdge(0) + h1->GetBinLowEdge(h1->GetNbinsX() / 2)) / 2 - tex.GetXsize(), h1->GetMaximum() / 2, "CLAS12 Preliminary");
        h1->Draw("same");
        h2->Draw("hsame");
    }
    else
    {
        h2->Draw();
        tex.DrawLatex((h2->GetBinLowEdge(0) + h2->GetBinLowEdge(h2->GetNbinsX() / 2)) / 2 - tex.GetXsize(), h2->GetMaximum() / 2, "CLAS12 Preliminary");
        h2->Draw("hsame");
        h1->Draw("same");
    }
}

Double_t funzione(Double_t *x, Double_t *par)
{
    return par[0] * TMath::Sin(x[0] / 180. * TMath::Pi()) / (1 + par[1] * TMath::Cos(x[0] / 180. * TMath::Pi()));
}

//===============================================================
void BSA_data_fit::DoTheJob(int mask = 1, string period = "fall2019", string energy = "", string photonTP = "", bool raw = false, bool sweighted = false, bool Q2binnedBSA = false, bool tbinnedBSA = false, bool xbbinnedBSA = false, bool PmissbinnedBSA = false, bool simultaneous = false, bool MC = false)
{

    using namespace RooFit;
    using namespace RooStats;

    bool do_CD_FD = false;

    bool binned = true;
    bool superimpose = false;
    if (raw == false)
    {
        superimpose = false;
    }

    bool MVACut = true;

    bool output = true; // produce output
    bool filter = true; // use filtered tuples : MANDATORY FOR UNBINNED FIT - WARNING  SELECTION-DEPENDENT
    bool minos = false;
    bool logY = false; // log-scale
    int nCPU = 2;
    bool w2Error = true; // MUST BE TRUE - EXCEPT WHEN SCALING ERRORS TO DATA STATs

    // dataset & fitting range
    double min = 0.;
    double max = 360.;
    int nbin = 12;

    double average_beam_polarization = 0.849759; // weighted average of three data sets (s19,f19,s20)

    if (period == "spring2019")
    {
        average_beam_polarization = 0.847658;
    }
    if (period == "fall2019")
    {
        average_beam_polarization = 0.84983;
    }
    if (period == "spring2020")
    {
        average_beam_polarization = 0.854082;
    }

    gROOT->Reset();

    TLatex tex; //
    tex.SetTextAlign(13);
    tex.SetTextColor(17);
    tex.SetTextSize(0.09);
    // tex.SetTextAngle(26.15998);
    tex.SetLineWidth(1);

    std::string purity = (raw) ? "raw_" : "BKG_free_";
    double HFlip = 1;

    /* if (period == "spring2019")
    {
        HFlip = -1;
    }*/
    if (period == "fall2019" || period == "spring2020")
    {
        energy = "10p4";
    }
    if (period == "all" && energy == "")
    {
        energy = "";
    }
    double tmin = -100, tmax = 100;
    double Q2min = -100, Q2max = 100;
    double xbmin = -100, xbmax = 100;
    double Pmissmin = -100, Pmissmax = 100;

    double sinPhi_ampl = 0;
    double sinPhi_ampl_err = 0;

    int tnbins = (mask & 1) ? 3 : 3;
    int Q2nbins = (mask & 1) ? 5 : 3;
    int xbnbins = (mask & 1) ? 4 : 3;

    double min_t[3] = {-100, -0.5, -0.3};
    double max_t[3] = {-0.5, -0.3, 0};

    double min_Q2[3] = {1, 1.9, 2.9};
    double max_Q2[3] = {1.9, 2.9, 100};

    double min_xb[3] = {0.05, 0.14, 0.2};
    double max_xb[3] = {0.14, 0.2, 100};

    // double min_t[2] = {-100, -0.29};
    // double max_t[2] = {-0.29, 0.};

    /*double min_t[3] = {-12.2452, -0.536325, -0.28278};
    double max_t[3] = {-0.536325, -0.28278, 0.126281};

    double min_Q2[5] = {1, 1.6, 2.0, 2.9, 3.8};
    double max_Q2[5] = {1.6, 2.0, 2.9, 3.8, 100};

    double min_xb[4] = {0.04, 0.13, 0.16, 0.2};
    double max_xb[4] = {0.13, 0.16, 0.2, 100};*/

    double min_Pmiss[7] = {0, 0.05, 0.1, 0.17, 0.3, 0.6, 0.85};
    double max_Pmiss[7] = {0.05, 0.1, 0.17, 0.3, 0.6, 0.85, 1};

    int t_count, Q2_count, xb_count, Pmiss_count;
    int binIndenter = 0;
    int plotIndenter = 1;
    t_count = (tbinnedBSA) ? sizeof(min_t) / sizeof(min_t[0]) : 1;
    Q2_count = (Q2binnedBSA) ? sizeof(min_Q2) / sizeof(min_Q2[0]) : 1;
    xb_count = (xbbinnedBSA) ? sizeof(min_xb) / sizeof(min_xb[0]) : 1;
    Pmiss_count = (PmissbinnedBSA) ? sizeof(min_Pmiss) / sizeof(min_Pmiss[0]) : 1;

    //----- TFile & TTree

    TChain *chain1 = (mask & 1) ? linked_pDVCS("pDVCS_stripped", period) : (TChain *)NULL; // RooTreeDataStore_dataS1_pDVCS_stripped", period) : (TChain *)NULL;
    TChain *chain2 = (mask & 2) ? linked_nDVCS("nDVCS_stripped", period) : (TChain *)NULL; // RooTreeDataStore_dataS1_pDVCS_stripped", period) : (TChain *)NULL;
    TChain *chain3 = (mask & 1) ? linked_DATA_eppi0("eppi0_stripped", period) : linked_DATA_enpi0("enpi0_stripped", period);
    TChain *chain4 = (mask & 1) ? linked_MC_eppi0("eppi0_stripped", energy, period) : linked_MC_enpi0("enpi0_stripped", energy);
    TChain *chain5 = (mask & 1) ? linked_MC_eppi01g("pDVCS_stripped", energy, period) : linked_MC_enpi01g("nDVCS_stripped", energy);

    tree1 = (mask & 1) ? (TTree *)chain1 : (TTree *)NULL; // pDVCS
    tree2 = (mask & 2) ? (TTree *)chain2 : (TTree *)NULL; // nDVCS
    tree3 = (TTree *)chain3;
    tree4 = (TTree *)chain4;
    tree5 = (TTree *)chain5;

    TCanvas *c11 = new TCanvas("c11", "c11");
    TCanvas *cbkgR = new TCanvas("cbkgR", "cbkgR");
    TCanvas *c11_cd = new TCanvas("c11_cd", "c11_cd");
    TCanvas *c11_fd = new TCanvas("c11_fd", "c11_fd");

    string binbin = "";

    int cxaxis, cyaxis;
    /*if (tbinnedBSA)
    {
        cxaxis = 3;
        cyaxis = 1;
    }
    if (Q2binnedBSA)
    {
        cxaxis = 3;
        cyaxis = 1;
    }
    if (xbbinnedBSA)
    {
        cxaxis = 3;
        cyaxis = 1;
    }

    if (Q2binnedBSA && xbbinnedBSA)
    {
        cxaxis = 3;
        cyaxis = 3;
    }*/

    if (tbinnedBSA)
    {
        cxaxis = 3;
        cyaxis = 1;
        binbin += "_t_";
    }
    if (Q2binnedBSA)
    {
        cxaxis = 3; // 5;
        cyaxis = 1;
        binbin += "_Q2_";
    }
    if (xbbinnedBSA)
    {
        cxaxis = 3; // 4;
        cyaxis = 1;
        binbin += "_Xb_";
    }

    if (Q2binnedBSA && xbbinnedBSA)
    {
        cxaxis = 4;
        cyaxis = 5;
    }

    if (PmissbinnedBSA)
    {
        cxaxis = 7;
        cyaxis = 1;
    }

    c11->Divide(cxaxis, cyaxis);
    cbkgR->Divide(cxaxis, cyaxis);
    c11_cd->Divide(cxaxis, cyaxis);
    c11_fd->Divide(cxaxis, cyaxis);

    std::string pref = (sweighted) ? "_sweighted_" : ""; // sweights?
    string channel_ALL = (mask & 1) ? "pDVCS" : "nDVCS";
    string channel_CD = (mask & 1) ? "pDVCS_cd" : "nDVCS_cd";
    string channel_FD = (mask & 1) ? "pDVCS_fd" : "nDVCS_fd";

    ofstream latex(std::string(channel_ALL + energy + purity + binbin + "_latex_plots_pi0comparisonDATAMCandBKGratioinData.dat").c_str());

    std::string Q2bins;
    std::string tbins;
    std::string xbbins;
    std::string Pmissbins;
    std::string bins_;

    if (Q2binnedBSA == true && tbinnedBSA == false && xbbinnedBSA == false)
    {
        bins_ = "Q2_bins";
    }
    else if (Q2binnedBSA == false && tbinnedBSA == true && xbbinnedBSA == false)
    {
        bins_ = "t_bins";
    }
    else if (Q2binnedBSA == false && tbinnedBSA == false && xbbinnedBSA == true)
    {
        bins_ = "Xbj_bins";
    }
    else
    {
        bins_ = "";
    }

    if (PmissbinnedBSA == true)
    {
        bins_ = "Pmiss_bins";
    }

    std::string D4bins = "";

    ofstream flux_sinPhi_ampl(std::string("sinPhi_ampl/" + period + "_" + purity + "BSA_" + channel_ALL + "_all_" + photonTP + "_binned_" + bins_ + pref + energy + ".dat").c_str());

    plotIndenter = 1;
    for (int i = 0; i < t_count; i++)
    {
        if (Q2binnedBSA && xbbinnedBSA)
            plotIndenter = 1;

        tmin = min_t[i];
        tmax = max_t[i];

        if (tbinnedBSA == false)
        {
            tmin = -100;
            tmax = 100.0; //-0.29
        }
        binIndenter = 0;
        for (int j = 0; j < Q2_count; j++)
        {
            Q2min = min_Q2[j];
            Q2max = max_Q2[j];

            if (Q2binnedBSA == false)
            {
                Q2min = -100;
                Q2max = 100;
            }
            cout << Q2min << "and" << Q2max << endl;

            if (j > 0 && j < 3 && xbbinnedBSA == true)
            {
                binIndenter++;
            }
            for (int k = 0; k < xb_count; k++)
            {

                xbmin = min_xb[k];
                xbmax = max_xb[k];

                if (xbbinnedBSA == false)
                {
                    xbmin = -100;
                    xbmax = 100;
                }

                for (int l = 0; l < Pmiss_count; l++)
                {

                    Pmissmin = min_Pmiss[l];
                    Pmissmax = max_Pmiss[l]; // Tag

                    if (PmissbinnedBSA == false)
                    {
                        Pmissmin = -100;
                        Pmissmax = 100;
                    }

                    Q2bins = (Q2binnedBSA) ? "_Q2bin_" + to_string(Q2min) + "_" + to_string(Q2max) : "";
                    tbins = (tbinnedBSA) ? "_tbin_" + to_string(tmin) + "_" + to_string(tmax) : "";
                    xbbins = (xbbinnedBSA) ? "_xbbin_" + to_string(xbmin) + "_" + to_string(xbmax) : "";
                    Pmissbins = (PmissbinnedBSA) ? "_Pmissbin_" + to_string(Pmissmin) + "_" + to_string(Pmissmax) : "";

                    D4bins = "";
                    D4bins += Q2bins;
                    D4bins += tbins;
                    D4bins += xbbins;
                    D4bins += Pmissbins;
                    D4bins += "_";
                    D4bins += energy;

                    string _Q2bins;
                    string _tbins;
                    string _xbbins;
                    string _D4bins;

                    _Q2bins = (Q2binnedBSA) ? " Q2bin " + to_string(Q2min) + " " + to_string(Q2max) : "";
                    _tbins = (tbinnedBSA) ? " tbin " + to_string(tmin) + " " + to_string(tmax) : "";
                    _xbbins = (xbbinnedBSA) ? " xbbin " + to_string(xbmin) + " " + to_string(xbmax) : "";

                    _D4bins = "";
                    _D4bins += _Q2bins;
                    _D4bins += _tbins;
                    _D4bins += _xbbins;

                    std::string suff = ""; // output suffix

                    //=================================
                    suff += "";
                    //================================

                    std::cout << "=============== Hb -> (hh)gamma fit =========== " << std::endl;
                    std::cout << " Usage :  timeFitMC([binned] ,[mask] , [trigger], [period] ,[simult])" << std::endl;
                    std::cout << " - mask    (int)   : define the set of channels to fit - see flags below (default : mask   = 1)" << std::endl;
                    std::cout << " - binned  (bool)  : perform binned  or unbinned fit (default : binned = true)" << std::endl;
                    std::cout << " - period  (string): select data sample '2011','2012' or 'any' (default : period ='all')" << std::endl;
                    std::cout << " - sim     (string): SIM version 'sim6', 'sim8', 'all' (default : sim='all')" << std::endl;
                    std::cout << " - simult  (bool)  : perform simultaneous or independant fits (default : simult = true)" << std::endl;
                    std::cout << " Mask  flags  : " << std::endl;

                    std::cout << "========================================= " << std::endl;
                    std::cout << " Fit mask = " << mask << std::endl;
                    if (!output)
                        std::cout << "No output will be stored" << std::endl;

                    if (0 == mask)
                    {
                        std::cout << "mask must be > 0" << std::endl;
                        simultaneous = false;
                    }

                    //====== Sanity checks
                    if (!binned && !filter)
                    {
                        std::cout << " === WARNING : filtered tuples must be used to perform unbinned fits - change 'filter' to true " << endl;
                        filter = true;
                    }

                    if ((mask == 1 || mask == 2 || mask == 4 || mask == 8) && simultaneous)
                    {
                        std::cout << " No need for simultaneous flag for simple fit - change 'simultaneous' to false" << endl;
                        simultaneous = false;
                    }

                    //==== Set fit flag
                    std::string flag = (binned) ? "_Binned_" : "_Unbinned_";
                    if (simultaneous)
                    {
                        std::ostringstream type("");
                        type << mask;
                        flag += "_simult" + type.str();
                    }

                    //==== Get data
                    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

                    if (minos)
                        std::cout << " ===== MINOS fit is requested =======" << std::endl;
                    else
                        std::cout << " ===== MINOS fit is disabled =======" << std::endl;

                    //====== Define Minos set
                    RooCmdArg myMinos;
                    if (minos)
                        // myMinos =  RooCmdArg(Minos(parList));
                        myMinos = RooCmdArg(Minos(true));
                    else
                        myMinos = RooCmdArg(Minos(false));

                    //----- define TCanvas
                    TCanvas *c1 = NULL;
                    TCanvas *c1_cd = NULL;
                    TCanvas *c1_fd = NULL;

                    c1 = new TCanvas((channel_ALL).c_str(), "");
                    prepareCanvas(c1, 1);
                    c1_cd = new TCanvas((channel_CD).c_str(), "");
                    prepareCanvas(c1_cd, 1);
                    c1_fd = new TCanvas((channel_FD).c_str(), "");
                    prepareCanvas(c1_fd, 1);

                    //----- Channel selection

                    TCut weight_mm2 = (MC) ? TCut("weighter(_mm2_eNg)") : TCut("1"); //_weight_mm2_epg;

                    TCut sWeights = (sweighted) ? TCut("sig_yield_sw") : TCut("1");
                    // TCut sWeights = (sweighted) ? TCut("bkg_yield_insignal_sw") : TCut("1");
                    // TCut sWeights = (sweighted) ? TCut("bkg_yield_comb_sw") : TCut("1");

                    TCut cut1_all = (mask & 1) ? processSelection(mask, "ALL", "") : TCut("0");
                    TCut cut1_cd = (mask & 1) ? processSelection(mask, "CD", "") : TCut("0");
                    TCut cut1_fd = (mask & 1) ? processSelection(mask, "FD", "") : TCut("0");

                    TCut cut2_all = (mask & 2) ? processSelection(mask, "ALL", "loose") : TCut("0");
                    TCut cut2_cd = (mask & 2) ? processSelection(mask, "CD", "tight") : TCut("0");
                    TCut cut2_fd = (mask & 2) ? processSelection(mask, "FD", "loose") : TCut("0");

                    TCut cut3_all = (mask & 1) ? processSelection(4, "ALL", "loose") : processSelection(8, "ALL", "loose");
                    TCut cut3_cd = (mask & 1) ? processSelection(4, "CD", "tight") : processSelection(8, "CD", "tight");
                    TCut cut3_fd = (mask & 1) ? processSelection(4, "FD", "loose") : processSelection(8, "FD", "loose");

                    string Q2bin = "_strip_Q2>" + to_string(Q2min) + "&&" + "_strip_Q2<" + to_string(Q2max);
                    string tbin = "_t_Nuc>" + to_string(tmin) + "&&" + "_t_Nuc<" + to_string(tmax);
                    string xbbin = "_strip_Xbj>" + to_string(xbmin) + "&&" + "_strip_Xbj<" + to_string(xbmax);
                    string Pmissbin = "TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal)>" + to_string(Pmissmin) + "&&" + "TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal)<" + to_string(Pmissmax);

                    TCut Q2bincut = TCut((Q2bin).c_str());
                    TCut tbincut = TCut((tbin).c_str());
                    TCut xbbincut = TCut((xbbin).c_str());
                    TCut Pmissbincut = TCut((Pmissbin).c_str());

                    if (Q2binnedBSA)
                    {
                        cut1_all += Q2bincut;
                        cut1_cd += Q2bincut;
                        cut1_fd += Q2bincut;
                        cut2_all += Q2bincut;
                        cut2_cd += Q2bincut;
                        cut2_fd += Q2bincut;
                        cut3_all += Q2bincut;
                        cut3_cd += Q2bincut;
                        cut3_fd += Q2bincut;
                    }
                    if (tbinnedBSA)
                    {
                        cut1_all += tbincut;
                        cut1_cd += tbincut;
                        cut1_fd += tbincut;
                        cut2_all += tbincut;
                        cut2_cd += tbincut;
                        cut2_fd += tbincut;
                        cut3_all += tbincut;
                        cut3_cd += tbincut;
                        cut3_fd += tbincut;
                    }
                    if (xbbinnedBSA)
                    {
                        cut1_all += xbbincut;
                        cut1_cd += xbbincut;
                        cut1_fd += xbbincut;
                        cut2_all += xbbincut;
                        cut2_cd += xbbincut;
                        cut2_fd += xbbincut;
                        cut3_all += xbbincut;
                        cut3_cd += xbbincut;
                        cut3_fd += xbbincut;
                    }

                    if (PmissbinnedBSA)
                    {
                        cut1_all += Pmissbincut;
                        cut1_cd += Pmissbincut;
                        cut1_fd += Pmissbincut;
                        cut2_all += Pmissbincut;
                        cut2_cd += Pmissbincut;
                        cut2_fd += Pmissbincut;
                        cut3_all += Pmissbincut;
                        cut3_cd += Pmissbincut;
                        cut3_fd += Pmissbincut;
                    }

                    TCut PhotonFT = TCut("_strip_Ph_status>=1000  && _strip_Ph_status<2000");
                    TCut PhotonEC = TCut("_strip_Ph_status>=2000");

                    if (photonTP == "FT")
                    {
                        cut1_all += PhotonFT;
                        cut1_cd += PhotonFT;
                        cut1_fd += PhotonFT;
                        cut2_all += PhotonFT;
                        cut2_cd += PhotonFT;
                        cut2_fd += PhotonFT;
                        cut3_all += PhotonFT;
                        cut3_cd += PhotonFT;
                        cut3_fd += PhotonFT;
                    }

                    if (photonTP == "EC")
                    {
                        cut1_all += PhotonEC;
                        cut1_cd += PhotonEC;
                        cut1_fd += PhotonEC;
                        cut2_all += PhotonEC;
                        cut2_cd += PhotonEC;
                        cut2_fd += PhotonEC;
                        cut3_all += PhotonEC;
                        cut3_cd += PhotonEC;
                        cut3_fd += PhotonEC;
                    }

                    TCut energy_period = TCut("");

                    if (energy == "10p6")
                    {
                        energy_period += TCut("_RunNumber<6420");
                    }
                    if (energy == "10p2")
                    {
                        energy_period += TCut("_RunNumber>=6420 && _RunNumber<10000");
                    }
                    if (energy == "10p4")
                    {
                        energy_period += TCut("_RunNumber>10000");
                    }

                    TCut proton_contamination_removal = TCut("_strip_Nuc_BDT>0.05");

                    if (MVACut == true)
                    {

                        cut2_all += proton_contamination_removal;
                        cut2_cd += proton_contamination_removal;
                        cut2_fd += proton_contamination_removal;
                        // cut3_all += proton_contamination_removal;
                        // cut3_cd += proton_contamination_removal;
                        // cut3_fd += proton_contamination_removal;
                    }

                    // to account for helicity flip
                    TCut pHEL = TCut("(_RunNumber <6700 && _RunNumber != 6378 && _Helicity<0) || (_RunNumber==6378 && _Helicity>0) || (_RunNumber >6700 && _Helicity>0)");
                    TCut nHEL = TCut("(_RunNumber <6700 && _RunNumber != 6378 && _Helicity>0) || (_RunNumber==6378 && _Helicity<0) || (_RunNumber >6700 && _Helicity<0)");

                    TCut polarization = TCut("weighterpolar(_RunNumber)"); //(_RunNumber <6700 && (1/0.847658)) || (_RunNumber>=6700 && _RunNumber<11300 && (1/0.84983)) || (_RunNumber >=11300 && (1/0.854082))");

                    //================= BUILD  binned dataset

                    DefineHistos(mask);

                    RooDataHist *dataH1 = NULL;
                    RooDataHist *dataH1_cd = NULL;
                    RooDataHist *dataH1_fd = NULL;

                    RooDataHist *dataH2 = NULL;

                    RooDataHist *dataH1_temp = NULL;
                    RooDataHist *dataH2_temp = NULL;

                    TH1F *ratioscale;
                    TH1F *ratiofinal, *ratiofinal_cd, *ratiofinal_fd;

                    RooRealVar x("x", "#Phi (degrees)", 0, 360);

                    //---------------------------------------------------------------------
                    int bins;
                    double syste = 0;

                    if (mask & 1)
                    {

                        bins = 12;

                        tree1->Project("all_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_all + pHEL + energy_period));
                        tree1->Project("all_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_all + nHEL + energy_period));

                        cout << "total number of events " << all_phi_cuts_hp->GetEntries() + all_phi_cuts_hn->GetEntries() << endl;
                        if (all_phi_cuts_hp->GetEntries() + all_phi_cuts_hn->GetEntries() == 0)
                        {
                            plotIndenter++;
                            continue;
                        }

                        tree5->Project("MC_Pi01g_all_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut1_all));
                        tree5->Project("MC_Pi01g_all_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut1_all));

                        tree4->Project("MC_Pi0_all_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_all));
                        tree4->Project("MC_Pi0_all_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_all));

                        tree3->Project("DATA_Pi0_all_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_all + pHEL + energy_period));
                        tree3->Project("DATA_Pi0_all_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_all + nHEL + energy_period));

                        MC_Pi01g_Pi0_all_phi_cuts_hp_ratio->Divide(MC_Pi01g_all_phi_cuts_hp, MC_Pi0_all_phi_cuts_hp, 1, 1);
                        MC_Pi01g_Pi0_all_phi_cuts_hn_ratio->Divide(MC_Pi01g_all_phi_cuts_hn, MC_Pi0_all_phi_cuts_hn, 1, 1);

                        // histos with no polarisation weight #######################################################################

                        tree1->Project("all_phi_cuts_hp_noP", "_Phi_Ph", sWeights * weight_mm2 * (cut1_all + pHEL + energy_period));
                        tree1->Project("all_phi_cuts_hn_noP", "_Phi_Ph", sWeights * weight_mm2 * (cut1_all + nHEL + energy_period));

                        tree3->Project("DATA_Pi0_all_phi_cuts_hp_noP", "_Phi_Pi0", weight_mm2 * (cut3_all + pHEL + energy_period));
                        tree3->Project("DATA_Pi0_all_phi_cuts_hn_noP", "_Phi_Pi0", weight_mm2 * (cut3_all + nHEL + energy_period));

                        DATA_Pi01g_all_phi_cuts_hp_noP->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hp_ratio, DATA_Pi0_all_phi_cuts_hp_noP, 1, 1);
                        DATA_Pi01g_all_phi_cuts_hn_noP->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hn_ratio, DATA_Pi0_all_phi_cuts_hn_noP, 1, 1);

                        all_phi_cuts_hp_BKGFree_noP->Add(all_phi_cuts_hp_noP, DATA_Pi01g_all_phi_cuts_hp_noP, 1, -1);
                        all_phi_cuts_hn_BKGFree_noP->Add(all_phi_cuts_hn_noP, DATA_Pi01g_all_phi_cuts_hn_noP, 1, -1);

                        // ##########################################################################################################

                        DATA_Pi01g_all_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hp_ratio, DATA_Pi0_all_phi_cuts_hp, 1, 1);
                        DATA_Pi01g_all_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hn_ratio, DATA_Pi0_all_phi_cuts_hn, 1, 1);

                        all_phi_cuts_hp_BKGFree->Add(all_phi_cuts_hp, DATA_Pi01g_all_phi_cuts_hp, 1, -1);
                        all_phi_cuts_hn_BKGFree->Add(all_phi_cuts_hn, DATA_Pi01g_all_phi_cuts_hn, 1, -1);

                        all_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_all_phi_cuts_hp, all_phi_cuts_hp, 1, 1);
                        all_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_all_phi_cuts_hn, all_phi_cuts_hn, 1, 1);

                        if (raw == true)
                        {
                            sum->Add(all_phi_cuts_hp_noP, all_phi_cuts_hn_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp, all_phi_cuts_hn, 1, -1);
                            ratio->Divide(diff, sum, 1, 1);
                        }
                        else
                        {
                            sum->Add(all_phi_cuts_hp_BKGFree_noP, all_phi_cuts_hn_BKGFree_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp_BKGFree, all_phi_cuts_hn_BKGFree, 1, -1);
                            ratio->Divide(diff, sum, 1, 1);
                        }

                        if (superimpose)
                        {
                            sum->Add(all_phi_cuts_hp_BKGFree_noP, all_phi_cuts_hn_BKGFree_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp_BKGFree, all_phi_cuts_hn_BKGFree, 1, -1);
                            ratio_bkgf->Divide(diff, sum, 1, 1);
                        }

                        if (do_CD_FD == true)
                        {

                            tree1->Project("cd_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_cd + pHEL + energy_period));
                            tree1->Project("cd_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_cd + nHEL + energy_period));

                            cout << "proton in CD, number of events " << cd_phi_cuts_hp->GetEntries() + cd_phi_cuts_hn->GetEntries() << endl;

                            tree5->Project("MC_Pi01g_cd_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut1_cd));
                            tree5->Project("MC_Pi01g_cd_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut1_cd));

                            tree4->Project("MC_Pi0_cd_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_cd));
                            tree4->Project("MC_Pi0_cd_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_cd));

                            tree3->Project("DATA_Pi0_cd_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_cd + pHEL + energy_period));
                            tree3->Project("DATA_Pi0_cd_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_cd + nHEL + energy_period));

                            MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio->Divide(MC_Pi01g_cd_phi_cuts_hp, MC_Pi0_cd_phi_cuts_hp, 1, 1);
                            MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio->Divide(MC_Pi01g_cd_phi_cuts_hn, MC_Pi0_cd_phi_cuts_hn, 1, 1);

                            DATA_Pi01g_cd_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio, DATA_Pi0_cd_phi_cuts_hp, 1, 1);
                            DATA_Pi01g_cd_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio, DATA_Pi0_cd_phi_cuts_hn, 1, 1);

                            cd_phi_cuts_hp_BKGFree->Add(cd_phi_cuts_hp, DATA_Pi01g_cd_phi_cuts_hp, 1, -1);
                            cd_phi_cuts_hn_BKGFree->Add(cd_phi_cuts_hn, DATA_Pi01g_cd_phi_cuts_hn, 1, -1);

                            cd_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_cd_phi_cuts_hp, cd_phi_cuts_hp, 1, 1);
                            cd_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_cd_phi_cuts_hn, cd_phi_cuts_hn, 1, 1);

                            if (raw == true)
                            {
                                cd_sum->Add(cd_phi_cuts_hp, cd_phi_cuts_hn, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp, cd_phi_cuts_hn, 1, -1);
                                cd_ratio->Divide(cd_diff, cd_sum, 1, 1);
                            }
                            else
                            {
                                cd_sum->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, 1, -1);
                                cd_ratio->Divide(cd_diff, cd_sum, 1, 1);
                            }

                            if (superimpose)
                            {
                                cd_sum->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, 1, -1);
                                cd_ratio_bkgf->Divide(cd_diff, cd_sum, 1, 1);
                            }

                            tree1->Project("fd_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_fd + pHEL + energy_period));
                            tree1->Project("fd_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut1_fd + nHEL + energy_period));

                            cout << "proton in FD, number of events " << fd_phi_cuts_hp->GetEntries() + fd_phi_cuts_hn->GetEntries() << endl;

                            tree5->Project("MC_Pi01g_fd_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut1_fd));
                            tree5->Project("MC_Pi01g_fd_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut1_fd));

                            tree4->Project("MC_Pi0_fd_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_fd));
                            tree4->Project("MC_Pi0_fd_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_fd));

                            tree3->Project("DATA_Pi0_fd_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_fd + pHEL + energy_period));
                            tree3->Project("DATA_Pi0_fd_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_fd + nHEL + energy_period));

                            MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio->Divide(MC_Pi01g_fd_phi_cuts_hp, MC_Pi0_fd_phi_cuts_hp, 1, 1);
                            MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio->Divide(MC_Pi01g_fd_phi_cuts_hn, MC_Pi0_fd_phi_cuts_hn, 1, 1);

                            DATA_Pi01g_fd_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio, DATA_Pi0_fd_phi_cuts_hp, 1, 1);
                            DATA_Pi01g_fd_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio, DATA_Pi0_fd_phi_cuts_hn, 1, 1);

                            fd_phi_cuts_hp_BKGFree->Add(fd_phi_cuts_hp, DATA_Pi01g_fd_phi_cuts_hp, 1, -1);
                            fd_phi_cuts_hn_BKGFree->Add(fd_phi_cuts_hn, DATA_Pi01g_fd_phi_cuts_hn, 1, -1);

                            fd_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_fd_phi_cuts_hp, fd_phi_cuts_hp, 1, 1);
                            fd_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_fd_phi_cuts_hn, fd_phi_cuts_hn, 1, 1);

                            if (raw == true)
                            {
                                fd_sum->Add(fd_phi_cuts_hp, fd_phi_cuts_hn, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp, fd_phi_cuts_hn, 1, -1);
                                fd_ratio->Divide(fd_diff, fd_sum, 1, 1);
                            }
                            else
                            {
                                fd_sum->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, 1, -1);
                                fd_ratio->Divide(fd_diff, fd_sum, 1, 1);
                            }

                            if (superimpose)
                            {
                                fd_sum->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, 1, -1);
                                fd_ratio_bkgf->Divide(fd_diff, fd_sum, 1, 1);
                            }
                        }
                    }
                    if (mask & 2)
                    {

                        bins = 7;

                        tree2->Project("all_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_all + pHEL + energy_period));
                        tree2->Project("all_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_all + nHEL + energy_period));

                        cout << "total number of events " << all_phi_cuts_hp->GetEntries() + all_phi_cuts_hn->GetEntries() << endl;
                        if (all_phi_cuts_hp->GetEntries() + all_phi_cuts_hn->GetEntries() == 0)
                        {
                            plotIndenter++;
                            continue;
                        }

                        tree5->Project("MC_Pi01g_all_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut2_all));
                        tree5->Project("MC_Pi01g_all_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut2_all));

                        tree4->Project("MC_Pi0_all_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_all));
                        tree4->Project("MC_Pi0_all_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_all));

                        tree3->Project("DATA_Pi0_all_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_all + pHEL + energy_period));
                        tree3->Project("DATA_Pi0_all_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_all + nHEL + energy_period));

                        MC_Pi01g_Pi0_all_phi_cuts_hp_ratio->Divide(MC_Pi01g_all_phi_cuts_hp, MC_Pi0_all_phi_cuts_hp, 1, 1);
                        MC_Pi01g_Pi0_all_phi_cuts_hn_ratio->Divide(MC_Pi01g_all_phi_cuts_hn, MC_Pi0_all_phi_cuts_hn, 1, 1);

                        // histos with no polarisation weight #######################################################################

                        tree2->Project("all_phi_cuts_hp_noP", "_Phi_Ph", sWeights * weight_mm2 * (cut2_all + pHEL + energy_period));
                        tree2->Project("all_phi_cuts_hn_noP", "_Phi_Ph", sWeights * weight_mm2 * (cut2_all + nHEL + energy_period));

                        tree3->Project("DATA_Pi0_all_phi_cuts_hp_noP", "_Phi_Pi0", weight_mm2 * (cut3_all + pHEL + energy_period));
                        tree3->Project("DATA_Pi0_all_phi_cuts_hn_noP", "_Phi_Pi0", weight_mm2 * (cut3_all + nHEL + energy_period));

                        DATA_Pi01g_all_phi_cuts_hp_noP->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hp_ratio, DATA_Pi0_all_phi_cuts_hp_noP, 1, 1);
                        DATA_Pi01g_all_phi_cuts_hn_noP->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hn_ratio, DATA_Pi0_all_phi_cuts_hn_noP, 1, 1);

                        all_phi_cuts_hp_BKGFree_noP->Add(all_phi_cuts_hp_noP, DATA_Pi01g_all_phi_cuts_hp_noP, 1, -1);
                        all_phi_cuts_hn_BKGFree_noP->Add(all_phi_cuts_hn_noP, DATA_Pi01g_all_phi_cuts_hn_noP, 1, -1);

                        // ##########################################################################################################

                        DATA_Pi01g_all_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hp_ratio, DATA_Pi0_all_phi_cuts_hp, 1, 1);
                        DATA_Pi01g_all_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_all_phi_cuts_hn_ratio, DATA_Pi0_all_phi_cuts_hn, 1, 1);

                        all_phi_cuts_hp_BKGFree->Add(all_phi_cuts_hp, DATA_Pi01g_all_phi_cuts_hp, 1, -1);
                        all_phi_cuts_hn_BKGFree->Add(all_phi_cuts_hn, DATA_Pi01g_all_phi_cuts_hn, 1, -1);

                        all_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_all_phi_cuts_hp, all_phi_cuts_hp, 1, 1);
                        all_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_all_phi_cuts_hn, all_phi_cuts_hn, 1, 1);

                        if (raw == true)
                        {

                            sum->Add(all_phi_cuts_hp_noP, all_phi_cuts_hn_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp, all_phi_cuts_hn, 1, -1);
                            ratio->Divide(diff, sum, 1, 1);
                        }
                        else
                        {
                            sum->Add(all_phi_cuts_hp_BKGFree_noP, all_phi_cuts_hn_BKGFree_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp_BKGFree, all_phi_cuts_hn_BKGFree, 1, -1);
                            ratio->Divide(diff, sum, 1, 1);
                        }

                        if (superimpose)
                        {
                            sum->Add(all_phi_cuts_hp_BKGFree_noP, all_phi_cuts_hn_BKGFree_noP, HFlip * 1, HFlip * 1);
                            diff->Add(all_phi_cuts_hp_BKGFree, all_phi_cuts_hn_BKGFree, 1, -1);
                            ratio_bkgf->Divide(diff, sum, 1, 1);
                        }
                        if (do_CD_FD == true)
                        {
                            tree2->Project("cd_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_cd + pHEL + energy_period));
                            tree2->Project("cd_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_cd + nHEL + energy_period));

                            cout << "neutron in CD, number of events " << cd_phi_cuts_hp->GetEntries() + cd_phi_cuts_hn->GetEntries() << endl;

                            tree5->Project("MC_Pi01g_cd_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut2_cd));
                            tree5->Project("MC_Pi01g_cd_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut2_cd));

                            tree4->Project("MC_Pi0_cd_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_cd));
                            tree4->Project("MC_Pi0_cd_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_cd));

                            tree3->Project("DATA_Pi0_cd_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_cd + pHEL + energy_period));
                            tree3->Project("DATA_Pi0_cd_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_cd + nHEL + energy_period));

                            MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio->Divide(MC_Pi01g_cd_phi_cuts_hp, MC_Pi0_cd_phi_cuts_hp, 1, 1);
                            MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio->Divide(MC_Pi01g_cd_phi_cuts_hn, MC_Pi0_cd_phi_cuts_hn, 1, 1);

                            DATA_Pi01g_cd_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio, DATA_Pi0_cd_phi_cuts_hp, 1, 1);
                            DATA_Pi01g_cd_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio, DATA_Pi0_cd_phi_cuts_hn, 1, 1);

                            cd_phi_cuts_hp_BKGFree->Add(cd_phi_cuts_hp, DATA_Pi01g_cd_phi_cuts_hp, 1, -1);
                            cd_phi_cuts_hn_BKGFree->Add(cd_phi_cuts_hn, DATA_Pi01g_cd_phi_cuts_hn, 1, -1);

                            cd_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_cd_phi_cuts_hp, cd_phi_cuts_hp, 1, 1);
                            cd_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_cd_phi_cuts_hn, cd_phi_cuts_hn, 1, 1);

                            if (raw == true)
                            {
                                cd_sum->Add(cd_phi_cuts_hp, cd_phi_cuts_hn, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp, cd_phi_cuts_hn, 1, -1);
                                cd_ratio->Divide(cd_diff, cd_sum, 1, 1);
                            }
                            else
                            {
                                cd_sum->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, 1, -1);
                                cd_ratio->Divide(cd_diff, cd_sum, 1, 1);
                            }

                            if (superimpose)
                            {
                                cd_sum->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                cd_diff->Add(cd_phi_cuts_hp_BKGFree, cd_phi_cuts_hn_BKGFree, 1, -1);
                                cd_ratio_bkgf->Divide(cd_diff, cd_sum, 1, 1);
                            }

                            tree2->Project("fd_phi_cuts_hp", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_fd + pHEL + energy_period));
                            tree2->Project("fd_phi_cuts_hn", "_Phi_Ph", polarization * sWeights * weight_mm2 * (cut2_fd + nHEL + energy_period));

                            cout << "neutron in FD, number of events " << fd_phi_cuts_hp->GetEntries() + fd_phi_cuts_hn->GetEntries() << endl;

                            tree5->Project("MC_Pi01g_fd_phi_cuts_hp", "_Phi_Ph", weight_mm2 * (cut2_fd));
                            tree5->Project("MC_Pi01g_fd_phi_cuts_hn", "_Phi_Ph", weight_mm2 * (cut2_fd));

                            tree4->Project("MC_Pi0_fd_phi_cuts_hp", "_Phi_Pi0", weight_mm2 * (cut3_fd));
                            tree4->Project("MC_Pi0_fd_phi_cuts_hn", "_Phi_Pi0", weight_mm2 * (cut3_fd));

                            tree3->Project("DATA_Pi0_fd_phi_cuts_hp", "_Phi_Pi0", polarization * weight_mm2 * (cut3_fd + pHEL + energy_period));
                            tree3->Project("DATA_Pi0_fd_phi_cuts_hn", "_Phi_Pi0", polarization * weight_mm2 * (cut3_fd + nHEL + energy_period));

                            MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio->Divide(MC_Pi01g_fd_phi_cuts_hp, MC_Pi0_fd_phi_cuts_hp, 1, 1);
                            MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio->Divide(MC_Pi01g_fd_phi_cuts_hn, MC_Pi0_fd_phi_cuts_hn, 1, 1);

                            DATA_Pi01g_fd_phi_cuts_hp->Multiply(MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio, DATA_Pi0_fd_phi_cuts_hp, 1, 1);
                            DATA_Pi01g_fd_phi_cuts_hn->Multiply(MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio, DATA_Pi0_fd_phi_cuts_hn, 1, 1);

                            fd_phi_cuts_hp_BKGFree->Add(fd_phi_cuts_hp, DATA_Pi01g_fd_phi_cuts_hp, 1, -1);
                            fd_phi_cuts_hn_BKGFree->Add(fd_phi_cuts_hn, DATA_Pi01g_fd_phi_cuts_hn, 1, -1);

                            fd_phi_cuts_hp_BKGRatio->Divide(DATA_Pi01g_fd_phi_cuts_hp, fd_phi_cuts_hp, 1, 1);
                            fd_phi_cuts_hn_BKGRatio->Divide(DATA_Pi01g_fd_phi_cuts_hn, fd_phi_cuts_hn, 1, 1);

                            if (raw == true)
                            {
                                fd_sum->Add(fd_phi_cuts_hp, fd_phi_cuts_hn, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp, fd_phi_cuts_hn, 1, -1);
                                fd_ratio->Divide(fd_diff, fd_sum, 1, 1);
                            }
                            else
                            {
                                fd_sum->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, 1, -1);
                                fd_ratio->Divide(fd_diff, fd_sum, 1, 1);
                            }

                            if (superimpose)
                            {
                                fd_sum->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, HFlip * 1, HFlip * 1);
                                fd_diff->Add(fd_phi_cuts_hp_BKGFree, fd_phi_cuts_hn_BKGFree, 1, -1);
                                fd_ratio_bkgf->Divide(fd_diff, fd_sum, 1, 1);
                            }
                        }
                    }
                    ofstream flux_nevents(std::string("bincontents/" + period + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + energy + "nevents.dat").c_str());
                    for (int kl = 1; kl <= ratio->GetNbinsX(); kl++)
                    {
                        flux_nevents << sum->GetBinContent(kl) << " ";
                    }
                    flux_nevents.close();

                    ofstream Flux_merging(std::string("datasets_merging/" + period + purity + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + energy + ".dat").c_str());

                    for (int kl = 1; kl <= ratio->GetNbinsX(); kl++)
                    {
                        Flux_merging << ratio->GetBinContent(kl) << " " << ratio->GetBinError(kl) << " ";
                    }
                    Flux_merging.close();

                    ofstream Flux(std::string("systematic_selection/" + period + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".dat").c_str());

                    for (int kl = 1; kl <= ratio->GetNbinsX(); kl++)
                    {
                        Flux << ratio->GetBinContent(kl) << " ";
                    }
                    Flux.close();
                    ifstream Fluxin(std::string("estimate_bkgSubtraction_systematics_BSA/systematics_selection/" + period + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".dat").c_str());

                    for (int kl = 0; kl <= bins; kl++)
                    {
                        Fluxin >> syste;
                        ratio_sys->SetBinContent(kl + 1, syste * ratio->GetBinContent(kl + 1));
                    }

                    Fluxin.close();
                    syste = 0;
                    ifstream Fluxin_sw(std::string("estimate_bkgSubtraction_systematics_BSA_sw/systematics_sw/" + period + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".dat").c_str());

                    for (int kl = 0; kl <= bins; kl++)
                    {
                        Fluxin_sw >> syste;
                        ratio_sys_sw->SetBinContent(kl + 1, syste * ratio->GetBinContent(kl + 1));
                    }

                    Fluxin_sw.close();
                    syste = 0;
                    ifstream Fluxin_BDT(std::string("estimate_systematics_BSA_BDT/systematics/" + period + "BSA_bins_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".dat").c_str());

                    for (int kl = 0; kl <= bins; kl++)
                    {
                        Fluxin_BDT >> syste;
                        ratio_sys_BDT->SetBinContent(kl + 1, syste * ratio->GetBinContent(kl + 1));
                    }

                    Fluxin_BDT.close();

                    gStyle->SetOptStat(0);
                    gStyle->SetOptFit(1);

                    c1->cd();
                    c1->SetLogy(true);
                    // MC_Pi0_all_phi_cuts_hp->Draw();
                    // DATA_Pi0_all_phi_cuts_hp->Draw("same");
                    superimposehist(MC_Pi0_all_phi_cuts_hp, DATA_Pi0_all_phi_cuts_hp);
                    c1->Print(("BKG_ratio_plots/singleplots/" + period + "_Pi0_MCDATA_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1->SetLogy(false);
                    c1->cd();
                    all_phi_cuts_hp_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/singleplots/" + period + "_hp_BKGRatioInData_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());

                    // produce latex format of plots to include in analysis note

                    latex << "\\begin{figure}[htb]" << endl;
                    latex << "\\begin{center}" << endl;
                    latex << "\\includegraphics[width=0.6\\textwidth]{new_plots/BSA_related_plots/" << channel_ALL << "/BKG_ratios/" << period << "_Pi0_MCDATA_" << channel_ALL << "_all_" << photonTP << D4bins << pref << ".pdf}" << endl;
                    latex << "\\includegraphics[width=0.6\\textwidth]{new_plots/BSA_related_plots/" << channel_ALL << "/BKG_ratios/" << period << "_hp_BKGRatioInData_" << channel_ALL << "_all_" << photonTP << D4bins << pref << ".pdf}" << endl;
                    latex << "\\caption {Channel: " << channel_ALL << " Bin: " << _D4bins << ". Top: data (blue) versus MC (red) comparison between the number of $\\pi^0$ events in each kinematical bin. Bottom: the fraction of $\\pi^0$ contamination in data per each kinematics bin.}" << endl;
                    latex << "\\label{" << D4bins << "}" << endl;
                    latex << "\\end{center}" << endl;
                    latex << "\\end{figure}" << endl;

                    c1->cd();
                    MC_Pi01g_all_phi_cuts_hp->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_1g_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    MC_Pi01g_cd_phi_cuts_hp->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_1g_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    MC_Pi01g_fd_phi_cuts_hp->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_1g_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    MC_Pi0_all_phi_cuts_hp->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_Pi0_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    MC_Pi0_cd_phi_cuts_hp->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_Pi0_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    MC_Pi0_fd_phi_cuts_hp->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_Pi0_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    MC_Pi01g_Pi0_all_phi_cuts_hp_ratio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_R_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_R_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_R_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    DATA_Pi0_all_phi_cuts_hp->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_dataPi0HP_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    DATA_Pi0_cd_phi_cuts_hp->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_dataPi0HP_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    DATA_Pi0_fd_phi_cuts_hp->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_dataPi0HP_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    DATA_Pi0_all_phi_cuts_hn->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_dataPi0HN_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    DATA_Pi0_cd_phi_cuts_hn->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_dataPi0HN_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    DATA_Pi0_fd_phi_cuts_hn->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_dataPi0HN_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    DATA_Pi01g_all_phi_cuts_hp->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hp_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    DATA_Pi01g_cd_phi_cuts_hp->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hp_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    DATA_Pi01g_fd_phi_cuts_hp->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hp_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    DATA_Pi01g_all_phi_cuts_hn->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hn_BKG_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    DATA_Pi01g_cd_phi_cuts_hn->Draw();
                    c1_cd->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hn_BKG_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    DATA_Pi01g_fd_phi_cuts_hn->Draw();
                    c1_fd->Print(("BKG_ratio_plots/" + period + "_dataPi01g_hn_BKG_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    all_phi_cuts_hp_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hp_BKGRatio_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    cd_phi_cuts_hp_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hp_BKGRatio_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    fd_phi_cuts_hp_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hp_BKGRatio_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    c1->cd();
                    all_phi_cuts_hn_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hn_BKGRatio_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    cd_phi_cuts_hn_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hn_BKGRatio_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    fd_phi_cuts_hn_BKGRatio->Draw();
                    c1->Print(("BKG_ratio_plots/" + period + "_hn_BKGRatio_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    ratio->SetMinimum(-0.4);
                    ratio->SetMaximum(0.4);
                    cd_ratio->SetMinimum(-0.4);
                    cd_ratio->SetMaximum(0.4);
                    fd_ratio->SetMinimum(-0.4);
                    fd_ratio->SetMaximum(0.4);

                    if (mask & 2)
                    {
                        ratio->SetMinimum(-0.15);
                        ratio->SetMaximum(0.15);
                        cd_ratio->SetMinimum(-0.15);
                        cd_ratio->SetMaximum(0.15);
                        fd_ratio->SetMinimum(-0.15);
                        fd_ratio->SetMaximum(0.15);
                    }

                    ratio->SetTitle("");
                    cd_ratio->SetTitle("");
                    fd_ratio->SetTitle("");

                    ratio->GetXaxis()->SetTitle("#phi (deg)");
                    cd_ratio->GetXaxis()->SetTitle("#phi (deg)");
                    fd_ratio->GetXaxis()->SetTitle("#phi (deg)");

                    ratio->GetYaxis()->SetTitle("A_{LU}");
                    cd_ratio->GetYaxis()->SetTitle("A_{LU}");
                    fd_ratio->GetYaxis()->SetTitle("A_{LU}");

                    ratio->GetXaxis()->SetTitleSize(0.055);
                    ratio->GetXaxis()->SetTitleOffset(0.9);
                    ratio->GetXaxis()->SetLabelSize(0.055);
                    ratio->GetYaxis()->SetTitleSize(0.055);
                    ratio->GetYaxis()->SetTitleOffset(0.8);
                    ratio->GetYaxis()->SetLabelSize(0.05);
                    ratio->GetYaxis()->CenterTitle(true);

                    TF1 *fitFcn = new TF1("fitFcn", funzione, 0, 360, 2);

                    fitFcn->SetParNames("a", "b");
                    fitFcn->SetParameter(0, 0.0);
                    fitFcn->SetParameter(1, 0.0);
                    fitFcn->SetLineColor(kBlue);

                    /* if (period == "fall2019" && xbmin == 0.17)
                    {
                        fitFcn->FixParameter(1, 0);
                    }*/

                    if (mask & 2)
                    {
                        fitFcn->FixParameter(1, 0);
                    }

                    if (Q2binnedBSA == true || tbinnedBSA == true || xbbinnedBSA == true || PmissbinnedBSA == true)
                    {

                        c11->cd(plotIndenter);

                        ratio->Fit("fitFcn");
                        ratio->Draw("");
                        ratio_sys->Draw("b SAME");
                        // ratio_sys_sw->Draw("b SAME");
                        ratio_sys_BDT->Draw("b SAME");
                        ratio->Draw("same");
                        if (superimpose)
                            ratio_bkgf->Draw("same");
                        tex.DrawLatex(5, -0.1, "CLAS12 Preliminary");
                        c11_cd->cd(plotIndenter);
                        // cd_ratio->Fit("fitFcn", "e");
                        cd_ratio->Draw();
                        ratio_sys->Draw("b SAME");
                        if (superimpose)
                            cd_ratio_bkgf->Draw("same");
                        c11_fd->cd(plotIndenter);
                        // fd_ratio->Fit("fitFcn", "e");
                        fd_ratio->Draw();
                        if (superimpose)
                            fd_ratio_bkgf->Draw("same");

                        cbkgR->cd(plotIndenter);
                        all_phi_cuts_hp_BKGRatio->Draw();

                        plotIndenter++;
                        // if (plotIndenter == 4)
                        //{
                        //     plotIndenter++;
                        // }
                        // if (plotIndenter == 7)
                        //{
                        //     plotIndenter += 2;
                        // }
                    }

                    c1->cd();
                    ratio->Fit("fitFcn", "e");
                    sinPhi_ampl = fitFcn->GetParameter(0);
                    sinPhi_ampl_err = fitFcn->GetParError(0);
                    ratio->Draw("");
                    ratio_sys->Draw("b SAME");
                    // ratio_sys_sw->Draw("b SAME");
                    ratio_sys_BDT->Draw("b SAME");
                    ratio->Draw("same");
                    // fitFcn->Draw("same");
                    if (superimpose)
                    {
                        // fitFcn->SetLineColor(kRed);
                        // ratio_bkgf->Fit("fitFcn");
                        ratio_bkgf->Draw("same");
                    }
                    tex.DrawLatex(5, -0.1, "CLAS12 Preliminary");
                    c1->Print((period + "_" + purity + "BSA_" + channel_ALL + "_all_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_cd->cd();
                    cd_ratio->Fit("fitFcn", "e");
                    cd_ratio->Draw();
                    ratio_sys->Draw("b SAME");
                    if (superimpose)
                    {
                        fitFcn->SetLineColor(kRed);
                        ratio_bkgf->Fit("fitFcn");
                        cd_ratio_bkgf->Draw("same");
                    }
                    c1_cd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_cd_" + photonTP + D4bins + pref + ".pdf").c_str());
                    c1_fd->cd();
                    fd_ratio->Fit("fitFcn", "e");
                    fd_ratio->Draw();
                    if (superimpose)
                    {
                        fitFcn->SetLineColor(kRed);
                        ratio_bkgf->Fit("fitFcn");
                        fd_ratio_bkgf->Draw("same");
                    }
                    c1_fd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_fd_" + photonTP + D4bins + pref + ".pdf").c_str());

                    flux_sinPhi_ampl << "t min " << tmin << " t max " << tmax << " Q2 min " << Q2min << " Q2 max " << Q2max << " xb min " << xbmin << " xb max " << xbmax << " sin phi ampl " << sinPhi_ampl << " " << sinPhi_ampl_err << endl;

                    /*
        ratioscale = new TH1F("ratioscale", "ratioscale", bins, 0, 360);
        ratioscale->Sumw2();
        for (int i = 0; i < bins; i++)
        {
            ratioscale->SetBinContent(i + 1, 3);
            ratioscale->SetBinError(i + 1, 0);
        }

        ratiofinal = new TH1F("ratiofinal", "ratiofinal", bins, 0, 360);
        ratiofinal->Sumw2();

        ratiofinal_cd = new TH1F("ratiofinal_cd", "ratiofinal_cd", bins, 0, 360);
        ratiofinal_cd->Sumw2();

        ratiofinal_fd = new TH1F("ratiofinal_fd", "ratiofinal_fd", bins, 0, 360);
        ratiofinal_fd->Sumw2();

        RooRealVar a("a^{blind}", "a", 0.15, -1., 1.);
        RooRealVar b("b^{blind}", "b", 0.00, -1., 1.);
        //blinding variables
        RooUnblindUniform *a_blind = new RooUnblindUniform("a_blind", "a_blind", "pdvcsrocks", 1, a);
        RooUnblindUniform *b_blind = new RooUnblindUniform("b_blind", "b_blind", "pdvcsrocks", 1, b);
        //blinding variables
        //RooMyPdf bsa("bsa","bsa",x,*a_blind,*b_blind) ;
        RooMyPdf bsa("bsa", "bsa", x, a, b);
        RooRealVar extend("extend", "extend", 10., 0., 1000.); //
        RooAddPdf finalbsa("finalbsa", "finalbsa", bsa, extend);

        RooRealVar a_cd("a_{cd}^{blind}", "a_cd", 0.15, -1., 1.);
        RooRealVar b_cd("b_{cd}^{blind}", "b_cd", 0.00, -1., 1.);
        //blinding variables
        RooUnblindUniform *a_cd_blind = new RooUnblindUniform("a_cd_blind", "a_cd_blind", "pdvcsrocks", 1, a_cd);
        RooUnblindUniform *b_cd_blind = new RooUnblindUniform("b_cd_blind", "b_cd_blind", "pdvcsrocks", 1, b_cd);
        //blinding variables
        //RooMyPdf bsa_cd("bsa_cd","bsa_cd",x,*a_cd_blind,*b_cd_blind) ;
        RooMyPdf bsa_cd("bsa_cd", "bsa_cd", x, a_cd, b_cd);
        RooRealVar extend_cd("extend_cd", "extend_cd", 10., 0., 1000.); //
        RooAddPdf finalbsa_cd("finalbsa_cd", "finalbsa_cd", bsa_cd, extend_cd);

        RooRealVar a_fd("a_{fd}^{blind}", "a_fd", 0.05, -1., 1.);
        RooRealVar b_fd("b_{fd}^{blind}", "b_fd", 0.00); //,-1.,1.);
        //blinding variables
        RooUnblindUniform *a_fd_blind = new RooUnblindUniform("a_fd_blind", "a_fd_blind", "pdvcsrocks", 1, a_fd);
        RooUnblindUniform *b_fd_blind = new RooUnblindUniform("b_fd_blind", "b_fd_blind", "pdvcsrocks", 1, b_fd);
        //blinding variables
        //RooMyPdf bsa_fd("bsa_fd","bsa_fd",x,*a_fd_blind,*b_fd_blind) ;
        RooMyPdf bsa_fd("bsa_fd", "bsa_fd", x, a_fd, b_fd);
        RooRealVar extend_fd("extend_fd", "extend_fd", 10., 0., 1000.); //
        RooAddPdf finalbsa_fd("finalbsa_fd", "finalbsa_fd", bsa_fd, extend_fd);

        //RooPlot* frame1 = x.frame(Title(ratio->GetTitle())) ;
        //RooPlot* frame1_cd = x.frame(Title(ratio->GetTitle())) ;
        //RooPlot* frame1_fd = x.frame(Title(ratio->GetTitle())) ;

        RooPlot *frame1 = x.frame(Title(" "));
        RooPlot *frame1_cd = x.frame(Title(" "));
        RooPlot *frame1_fd = x.frame(Title(" "));

        // all
        ratiofinal->Add(ratio, ratioscale, 1, 1);

        dataH1 = new RooDataHist("dataH1", "dataH1", x, ratiofinal);
        //dataH1_temp = new RooDataHist("dataH1_temp","dataH1_temp", x, ratio);
        //RooMyPdftemp bsa_temp("bsa_temp","bsa_temp",x,a,b) ;
        //RooAddPdf finalbsa_temp("finalbsa_temp", "finalbsa_temp", bsa_temp, extend);

        std::cout << *dataH1 << std::endl;

        RooFitResult *result1 = finalbsa.fitTo(*dataH1, Extended(kTRUE), Save(), myMinos, NumCPU(nCPU), SumW2Error(w2Error));

        rooDisplayFit(x, dataH1, frame1, c1, finalbsa, result1, (period+"_"+purity + "BSA_" + channel_ALL + D4bins).c_str(), RooArgList(x, a, b), false, "SumW2", RooCmdArg(), RooFit::Layout(0.55, 0.8, 0.9));

        //cd

        ratiofinal_cd->Add(cd_ratio, ratioscale, 1, 1);

        dataH1_cd = new RooDataHist("dataH1_cd", "dataH1_cd", x, ratiofinal_cd);
        std::cout << *dataH1_cd << std::endl;

        RooFitResult *result1_cd = finalbsa_cd.fitTo(*dataH1_cd, Extended(kTRUE), Save(), myMinos, NumCPU(nCPU), SumW2Error(w2Error));

        rooDisplayFit(x, dataH1_cd, frame1_cd, c1_cd, finalbsa_cd, result1_cd, (period+"_"+purity + "BSA_" + channel_ALL + "_cd" + D4bins).c_str(), RooArgList(x, a_cd, b_cd), false, "SumW2", RooCmdArg(), RooFit::Layout(0.55, 0.8, 0.9));

        //fd

        ratiofinal_fd->Add(fd_ratio, ratioscale, 1, 1);

        dataH1_fd = new RooDataHist("dataH1_fd", "dataH1_fd", x, ratiofinal_fd);
        std::cout << *dataH1 << std::endl;

        RooFitResult *result1_fd = finalbsa_fd.fitTo(*dataH1_fd, Extended(kTRUE), Save(), myMinos, NumCPU(nCPU), SumW2Error(w2Error));

        cout << result1_fd << endl;

        rooDisplayFit(x, dataH1_fd, frame1_fd, c1_fd, finalbsa_fd, result1_fd, (period+"_"+purity + "BSA_" + channel_ALL + "_fd" + D4bins).c_str(), RooArgList(x, a_fd, b_fd), false, "SumW2", RooCmdArg(), RooFit::Layout(0.55, 0.8, 0.9));
        */
                }
            }
        }
        if (Q2binnedBSA == true && tbinnedBSA == true && xbbinnedBSA == true)
        {

            c11->Print((period + "_" + purity + "BSA_" + channel_ALL + "_all_" + photonTP + "_binned_" + bins_ + tbins + pref + energy + ".pdf").c_str());
            cbkgR->Print((period + "_" + purity + "R_BKG_" + channel_ALL + "_all_" + photonTP + "_binned_" + bins_ + tbins + pref + energy + ".pdf").c_str());

            c11_cd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_cd_" + photonTP + "_binned_" + bins_ + tbins + pref + energy + ".pdf").c_str());
            c11_fd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_fd_" + photonTP + "_binned_" + bins_ + tbins + pref + energy + ".pdf").c_str());
        }
    }
    if (Q2binnedBSA == true || tbinnedBSA == true || xbbinnedBSA == true || PmissbinnedBSA == true)
    {

        c11->Print((period + "_" + purity + "BSA_" + channel_ALL + "_all_" + photonTP + "_binned_" + bins_ + pref + energy + ".pdf").c_str());
        cbkgR->Print((period + "_" + purity + "R_BKG_" + channel_ALL + "_all_" + photonTP + "_binned_" + bins_ + pref + energy + ".pdf").c_str());
        c11_cd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_cd_" + photonTP + "_binned_" + bins_ + pref + energy + ".pdf").c_str());
        c11_fd->Print((period + "_" + purity + "BSA_" + channel_ALL + "_fd_" + photonTP + "_binned_" + bins_ + pref + energy + ".pdf").c_str());
    }
    return;
}
