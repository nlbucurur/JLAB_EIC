#ifndef BSA_data_fit_h
#define BSA_data_fit_h

// #include "rooLib.C"
// #include "getTreeMC.C"
#include <TMath.h>
#include "RooNumConvPdf.h"
#include "RooMultiVarGaussian.h"
#include "RooUnblindUniform.h"
#include "RooAddition.h"
#include "RooMinuit.h"

#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <cmath>
// ROOT includes
#include <TTree.h>
#include <TObject.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TCut.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TPaveStats.h>
#include <TRandom.h>
// RooFit includes
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooAbsArg.h>
#include <RooDataSet.h>
#include <RooFormulaVar.h>
#include <RooCategory.h>
#include <RooGaussian.h>
#include <RooVoigtian.h>
#include <RooExtendPdf.h>
#include <RooCBShape.h>
#include <RooEffProd.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooAbsArg.h>
#include <RooAbsPdf.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooSimultaneous.h>
#include <RooArgusBG.h>
#include "RooStats/SPlot.h"
#include <RooCategory.h>
#include <RooStringVar.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooMCStudy.h>
#include <TPaveText.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include "RooKeysPdf.h"
#include "RooTFnBinding.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include <TCanvas.h>
#include <TGraphErrors.h>
// ROOSTATS
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "chain_pDVCS.C"
// #include "chain_nDVCS.C"
// #include "chain_eppi0_data.C"
// #include "chain_eppi0_MC.C"
#include "chain_eppi01gamma_MC.C"
// #include "chain_enpi0_data.C"
// #include "chain_enpi0_MC.C"
// #include "chain_enpi01gamma_MC.C"

class BSA_data_fit
{
public:
    TTree *tree1;
    TTree *tree2;
    TTree *tree3;
    TTree *tree4;
    TTree *tree5;
    TTree *tree8;
    TTree *tree;

    TFile *rootFileHistos;
    TFile *fin;

    TH1F *all_phi_cuts_hp_noP, *all_phi_cuts_hn_noP, *DATA_Pi0_all_phi_cuts_hp_noP, *DATA_Pi0_all_phi_cuts_hn_noP, *all_phi_cuts_hp_BKGFree_noP, *all_phi_cuts_hn_BKGFree_noP, *DATA_Pi01g_all_phi_cuts_hp_noP, *DATA_Pi01g_all_phi_cuts_hn_noP;

    TH1F *all_phi_cuts_hp, *all_phi_cuts_hn, *sum, *diff, *ratio, *ratio_sys, *ratio_sys_sw, *ratio_sys_BDT;
    TH1F *cd_phi_cuts_hp, *cd_phi_cuts_hn, *cd_sum, *cd_diff, *cd_ratio;
    TH1F *fd_phi_cuts_hp, *fd_phi_cuts_hn, *fd_sum, *fd_diff, *fd_ratio;

    TH1F *ratio_bkgf;
    TH1F *cd_ratio_bkgf;
    TH1F *fd_ratio_bkgf;

    TH1F *DATA_Pi0_all_phi_cuts_hp, *DATA_Pi0_all_phi_cuts_hn;
    TH1F *DATA_Pi0_cd_phi_cuts_hp, *DATA_Pi0_cd_phi_cuts_hn;
    TH1F *DATA_Pi0_fd_phi_cuts_hp, *DATA_Pi0_fd_phi_cuts_hn;

    TH1F *MC_Pi0_all_phi_cuts_hp, *MC_Pi0_all_phi_cuts_hn;
    TH1F *MC_Pi0_cd_phi_cuts_hp, *MC_Pi0_cd_phi_cuts_hn;
    TH1F *MC_Pi0_fd_phi_cuts_hp, *MC_Pi0_fd_phi_cuts_hn;

    TH1F *MC_Pi01g_all_phi_cuts_hp, *MC_Pi01g_all_phi_cuts_hn;
    TH1F *MC_Pi01g_cd_phi_cuts_hp, *MC_Pi01g_cd_phi_cuts_hn;
    TH1F *MC_Pi01g_fd_phi_cuts_hp, *MC_Pi01g_fd_phi_cuts_hn;

    TH1F *MC_Pi01g_Pi0_all_phi_cuts_hp_ratio, *MC_Pi01g_Pi0_all_phi_cuts_hn_ratio;
    TH1F *MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio, *MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio;
    TH1F *MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio, *MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio;

    TH1F *DATA_Pi01g_all_phi_cuts_hp, *DATA_Pi01g_all_phi_cuts_hn;
    TH1F *DATA_Pi01g_cd_phi_cuts_hp, *DATA_Pi01g_cd_phi_cuts_hn;
    TH1F *DATA_Pi01g_fd_phi_cuts_hp, *DATA_Pi01g_fd_phi_cuts_hn;

    TH1F *all_phi_cuts_hp_BKGFree, *all_phi_cuts_hn_BKGFree;
    TH1F *cd_phi_cuts_hp_BKGFree, *cd_phi_cuts_hn_BKGFree;
    TH1F *fd_phi_cuts_hp_BKGFree, *fd_phi_cuts_hn_BKGFree;

    TH1F *all_phi_cuts_hp_BKGRatio, *all_phi_cuts_hn_BKGRatio;
    TH1F *cd_phi_cuts_hp_BKGRatio, *cd_phi_cuts_hn_BKGRatio;
    TH1F *fd_phi_cuts_hp_BKGRatio, *fd_phi_cuts_hn_BKGRatio;

    BSA_data_fit();

    virtual ~BSA_data_fit();
    virtual RooHist *myResidHist(const RooHist &data, const RooCurve &curve, bool normalize, int method, RooAddPdf pdf, RooRealVar M);
    virtual int prepareCanvas(TCanvas *canvas, int addResiduals);
    virtual void plotResiduals(RooRealVar M, RooPlot &frame, TString flag, RooAddPdf pdf);
    virtual void rooDisplay(RooRealVar x, RooDataHist *dataH, RooPlot *frame, TCanvas *c, std::string out, bool logY, std::string dataError, const RooCmdArg &range, const RooCmdArg &cmd1);
    virtual void rooDisplayFit(RooRealVar x, RooDataHist *dataH, RooPlot *frame, TCanvas *c, RooAddPdf pdf, RooFitResult *result, std::string out, RooArgList pars, bool logY, std::string dataError, const RooCmdArg &range, const RooCmdArg &cmd1);

    virtual void DefineHistos(int selflag);
    virtual TCut ALL_cuts(int selflag);
    virtual TCut processSelection(int selFlag, std::string DetectorRegion, std::string selCriteria);

    virtual void DoTheJob(int mask, string period, string energy, string photonTP, bool raw, bool sw, bool Q2binnedBSA, bool tbinnedBSA, bool xbbinnedBSA, bool PmissbinnedBSA, bool simultaneous, bool MC);

    // virtual TCut weight();
};

BSA_data_fit::BSA_data_fit()
{
    return;
}

BSA_data_fit::~BSA_data_fit()
{
    return;
}

TWeighter1D tweighterpDVCS_mm2_eNg("weight_mm2_eNg_1D.root", "c_mm2_eNg_cuts_weights");

double weighter(double val)
{
    double ww = 0.;
    ww = tweighterpDVCS_mm2_eNg(val);
    // cout<<ww<<endl;
    return ww;
}

double weighterpolar(double val)
{
    double ww = 0.;

    if (val < 5600)
    {
        ww = 1 / 0.86;
    }

    else if (val >= 6142 && val <= 6149)
    {
        ww = 1 / 0.81132;
    }

    else if (val >= 6150 && val <= 6188)
    {
        ww = 1 / 0.82137;
    }

    else if (val >= 6189 && val <= 6260)
    {
        ww = 1 / 0.83598;
    }

    else if (val >= 6261 && val <= 6339)
    {
        ww = 1 / 0.80770;
    }

    else if (val >= 6340 && val <= 6342)
    {
        ww = 1 / 0.85536;
    }

    else if (val >= 6343 && val <= 6399)
    {
        ww = 1 / 0.87038;
    }

    else if (val >= 6400 && val <= 6476)
    {
        ww = 1 / 0.88214;
    }

    else if (val >= 6477 && val <= 6532)
    {
        ww = 1 / 0.86580;
    }

    else if (val >= 6533 && val <= 11000)
    {
        ww = 1 / 0.87887;
    }

    else if (val >= 11010 && val <= 11310)
    {
        ww = 1 / 0.84983;
    }

    else if (val >= 11323 && val <= 11334)
    {
        ww = 1 / 0.87135;
    }

    else if (val >= 11335 && val <= 11387)
    {
        ww = 1 / 0.85048;
    }

    else if (val >= 11388 && val <= 11541)
    {
        ww = 1 / 0.84262;
    }

    else if (val >= 11542)
    {
        ww = 1 / 0.85188;
    }

    else
    {
        ww = 1 / 0.849759;
    }

    // cout << ww << endl;
    return ww;
}

int BSA_data_fit::prepareCanvas(TCanvas *canvas, int addResiduals)
{
    if (addResiduals == 1)
    {
        // Prepare canvas
        canvas->Divide(1, 2);
        TPad *padHisto = (TPad *)canvas->GetListOfPrimitives()->At(0);
        TPad *padResid = (TPad *)canvas->GetListOfPrimitives()->At(1);
        double small = 0.05;
        double r = 0.;
        if (padResid)
            r = .2;
        padHisto->SetPad(0., r, 1., 1.);
        padHisto->SetBottomMargin(small);
        if (padResid)
        {
            padHisto->SetBottomMargin(1.0);
            padResid->SetPad(0., 0., 1., r);
            padResid->SetBottomMargin(0.25);
            padResid->SetTopMargin(small /*+0.02*/);
        }
    }
    return 0;
}

RooHist *BSA_data_fit::myResidHist(const RooHist &data, const RooCurve &curve, bool normalize, int method, RooAddPdf pdf, RooRealVar M)
{

    // Create and return RooHist containing  residuals w.r.t to given curve.
    // If normalize is true, the residuals are normalized by the histogram
    // errors creating a RooHist with pull values

    // Copy all non-content properties from hist1
    RooHist *hist = new RooHist(data.getNominalBinWidth());
    if (normalize)
    {
        hist->SetName(Form("pull_%s_%s", data.GetName(), curve.GetName()));
        hist->SetTitle(Form("Pull of %s and %s", data.GetTitle(), curve.GetTitle()));
    }
    else
    {
        hist->SetName(Form("resid_%s_%s", data.GetName(), curve.GetName()));
        hist->SetTitle(Form("Residual of %s and %s", data.GetTitle(), curve.GetTitle()));
    }

    // Determine range of curve
    Double_t xstart, xstop, y;
    curve.GetPoint(0, xstart, y);
    curve.GetPoint(curve.GetN() - 1, xstop, y);

    // Add histograms, calculate Poisson confidence interval on sum value
    double m = 0;
    double m2 = 0;
    double n = 0;

    //
    double sum = 0;
    if (method == 3)
    {
        const RooArgList &coefs = pdf.coefList();
        RooLinkedListIter it = coefs.iterator();
        RooRealVar *coef;
        while ((coef = (RooRealVar *)it.Next()))
        {
            sum += coef->getVal();
        }
    }

    //
    for (Int_t i = 0; i < data.GetN(); i++)
    {
        Double_t x, point;
        data.GetPoint(i, x, point);

        // Only calculate pull for bins inside curve range
        if (x < xstart || x > xstop)
            continue;

        Double_t yy;
        if (method == 1)
        {
            Double_t exl = 0.5 * data.getNominalBinWidth();
            Double_t exh = 0.5 * data.getNominalBinWidth();
            yy = point - curve.average(x - exl, x + exh);
        }
        else if (method == 2)
        {
            yy = point - curve.interpolate(x);
            // std::cout << " " << i << " " << x << " " << point << " " << curve.interpolate(x) << std::endl;
        }
        else if (method == 3)
        {
            Double_t exl = 0.5 * data.getNominalBinWidth();
            Double_t exh = 0.5 * data.getNominalBinWidth();
            M.setRange("bin", x - exl, x + exh);
            RooAbsReal *integ = pdf.createIntegral(M, RooFit::NormSet(M), RooFit::Range("bin"));
            // std::cout << " Integral " << integ->getVal() << " * " << sum << std::endl;
            yy = point - integ->getVal() * sum;
        }
        else
        {
            std::cout << "unknown method " << std::endl;
            return NULL;
        }

        Double_t dyl = data.GetErrorYlow(i);
        Double_t dyh = data.GetErrorYhigh(i);

        // std::cout << " --> " << x << " " << point << " " << yy  << " " << dyl << " " << dyh << std::endl;

        if (normalize)
        {
            Double_t norm = (yy > 0 ? dyl : dyh);
            if (norm == 0.)
            {
                //          std::cout << "RooHist::makeResisHist(" << data.GetName() << ") WARNING: point " << i << " has zero error, setting residual to zero" << std::endl ;
                yy = 0;
                dyh = 0;
                dyl = 0;
            }
            else
            {
                yy /= norm;
                dyh /= norm;
                dyl /= norm;
            }

            m += yy;
            m2 += yy * yy;
            n += 1.;
        }
        hist->addBinWithError(x, yy, dyl, dyh);
    }
    if (n > 0)
    {
        double mean = m / n;
        double rms = sqrt(m2 / n - m * m / n / n);
        double emean = mean / sqrt(n);
        std::cout << " PULL <mu>= " << mean << " ± " << emean << " - rms = " << rms << " -   chi2 = " << m2 << std::endl;
    }

    return hist;
}

void BSA_data_fit::plotResiduals(RooRealVar M, RooPlot &frame, TString flag, RooAddPdf pdf)
{
    using namespace RooFit;
    using namespace RooStats;
    TString sdata = flag + "Data";
    TString scurve = flag + "Curve";

    const RooCurve *curve = frame.getCurve(scurve);
    const RooHist *data = frame.getHist(sdata);
    RooHist *residuals = myResidHist(*data, *curve, true, 2, pdf, M);
    double xMin = frame.GetXaxis()->GetXmin();
    double xMax = frame.GetXaxis()->GetXmax();
    RooPlot *fpull = M.frame(Title(" "), Range(xMin, xMax));
    fpull->addPlotable(residuals, "E3");
    fpull->SetMinimum(-5.);
    fpull->SetMaximum(5.);
    residuals->SetFillColor(13);
    // residuals->SetMarkerColor( kBlue);
    // residuals->SetMarkerStyle( 20    );
    // residuals->SetMarkerSize ( .8    );
    TAxis *xAxis = fpull->GetXaxis();
    xAxis->SetTickLength(5 * xAxis->GetTickLength());
    xAxis->SetLabelSize(5 * 0.03 /*xAxis->GetLabelSize()*/);
    xAxis->SetTitleSize(5 * 0.04 /*xAxis->GetTitleSize()*/);
    xAxis->SetLabelOffset(5 * xAxis->GetLabelOffset());
    xAxis->SetTitle("");
    TLine *midLine = new TLine(xMin, 0., xMax, 0.);
    TLine *uppLine = new TLine(xMin, 2., xMax, 2.);
    TLine *lowLine = new TLine(xMin, -2., xMax, -2.);
    uppLine->SetLineColor(kRed);
    lowLine->SetLineColor(kRed);
    TAxis *yAxis = fpull->GetYaxis();
    yAxis->SetNdivisions(504);
    yAxis->SetLabelSize(5 * 0.03 /*yAxis->GetLabelSize()*/);
    fpull->Draw();
    uppLine->Draw("same");
    midLine->Draw("same");
    lowLine->Draw("same");

    TAxis *axis = frame.GetXaxis();
    frame.GetYaxis()->SetTitleSize(0.05);
    frame.GetYaxis()->SetLabelSize(0.04);
    TString myTitle = axis->GetTitle();
    axis->SetTitleOffset(0.5);
    axis->SetLabelSize(0);
    axis->SetTitleSize(0.05);
    axis->CenterTitle();
}

void BSA_data_fit::rooDisplayFit(RooRealVar x, RooDataHist *dataH, RooPlot *frame, TCanvas *c, RooAddPdf pdf, RooFitResult *result, std::string out = "",
                                 RooArgList pars = RooArgList(),
                                 bool logY = false,
                                 std::string dataError = "SumW2",
                                 const RooCmdArg &range = RooCmdArg(),
                                 const RooCmdArg &cmd1 = RooFit::Layout(0.55, 0.8, 0.9))
{
    using namespace RooFit;
    TPad *padHisto = NULL;
    TPad *padResid = NULL;
    if (c->GetListOfPrimitives()->LastIndex() > 0)
    {
        padResid = (TPad *)c->GetListOfPrimitives()->At(1);
        padHisto = (TPad *)c->GetListOfPrimitives()->At(0);
        std::cout << "Prepare for residuals" << endl;
        padHisto->cd();
    }

    const static int ncol = 10;
    const static Color_t cols[ncol] = {kRed, kGreen, kMagenta, kCyan, kBlue, kYellow, kPink, kViolet, kTeal, kSpring};
    const static int nsty = 10;
    const static Style_t stys[nsty] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    if (dataError == "Poisson")
        dataH->plotOn(frame, Name("Data"), DataError(RooAbsData::Poisson), range);
    else
    {
        if (dataError != "SumW2")
            std::cout << "Unknown data error type - assume SumW2" << endl;
        dataH->plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2), range);
    }

    // display total pdf
    pdf.plotOn(frame, VisualizeError(*result), FillColor(kBlue), range); // band
    pdf.plotOn(frame, LineColor(kBlue + 2), Name("Curve"), range);       // curve

    // display pdf component(s)
    printf("Display components \n");
    const RooArgList &pdfs = pdf.pdfList();
    const RooArgList &coefs = pdf.coefList();
    RooLinkedListIter it = pdfs.iterator();
    RooAbsPdf *component;
    int k = 0;
    const RooAbsPdf &sigPdf = (RooAbsPdf &)pdfs[0]; // Warning : assume signal PDF is the 1st entry in RooAddPdf
    const RooRealVar &coef = (RooRealVar &)coefs[0];
    while ((component = (RooAbsPdf *)it.Next()))
    {
        Color_t col = (k < ncol) ? cols[k] : cols[ncol - 1];
        Style_t sty = (k < nsty) ? stys[k] : stys[nsty - 1];
        pdf.plotOn(frame, Components(*component), VisualizeError(*result), FillColor(col), range);
        pdf.plotOn(frame, Components(*component), LineStyle(sty), LineColor(col + 2), range);

        k++;
    }
    if (dataError == "Poisson")
    {
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
        dataH->plotOn(frame, DataError(RooAbsData::Poisson));
    }
    else
        dataH->plotOn(frame, DataError(RooAbsData::SumW2));

    if (logY)
    {
        frame->SetMinimum(0.2);
        if (padHisto)
            padHisto->SetLogy();
        else
            c->SetLogy();
    }
    else
        pdf.paramOn(frame, RooFit::Parameters(pars), cmd1);

    printf("Draw main frame \n");
    frame->GetYaxis()->SetTitleOffset(0.72);
    frame->GetYaxis()->SetTitle("A_{LU}(#Phi)");
    frame->Draw();
    c->Update();

    if (padResid)
    {
        padResid->cd();
        printf("Draw residuals \n");
        plotResiduals(x, *frame, "", pdf);
    }
    if (out != "")
    {
        c->Print(std::string(out + ".C").c_str());
        c->Print(std::string(out + ".pdf").c_str());
    }

    // S/B value
    printf("Compute S/B \n");
    if (pars.getSize() >= 2)
    {
        const RooRealVar &mu = (RooRealVar &)pars[0];  // Assume mean is the 1st parameters
        const RooRealVar &sig = (RooRealVar &)pars[1]; // Assume sigma is the 2nd parameter
        double xmin = mu.getVal() - 3. * sig.getVal();
        double xmax = mu.getVal() + 3. * sig.getVal();
        x.setRange("xRange", xmin, xmax);
        RooAbsReal *isig = sigPdf.createIntegral(x, NormSet(x), RooFit::Range("xRange"));
        RooAbsReal *iall = pdf.createIntegral(x, NormSet(x), RooFit::Range("xRange"));
        double S = isig->getVal() * coef.getVal();
        const RooArgSet *nset = new RooArgSet(x);
        double B = iall->getVal() * pdf.expectedEvents(nset) - S;
        double SoB = (B != 0) ? S / B : 0;
        std::cout << " ====== S/B[mean+/-3*sigma] = S/B[" << xmin << "," << xmax << "]  = " << S << "/" << B << " = " << SoB << std::endl;
    }
    // approximate chi2
    int npar = result->floatParsFinal().getSize();
    int nbin = dataH->numEntries();
    double chi2 = frame->chiSquare("Curve", "Data", npar) * (nbin - npar);
    // return chi2;

    return;
}

void BSA_data_fit::rooDisplay(RooRealVar x, RooDataHist *dataH, RooPlot *frame, TCanvas *c, std::string out = "",
                              bool logY = false,
                              std::string dataError = "SumW2",
                              const RooCmdArg &range = RooCmdArg(),
                              const RooCmdArg &cmd1 = RooFit::Layout(0.15, 0.5, 0.9))
{
    using namespace RooFit;
    TPad *padHisto = NULL;
    TPad *padResid = NULL;
    if (c->GetListOfPrimitives()->LastIndex() > 0)
    {
        padResid = (TPad *)c->GetListOfPrimitives()->At(1);
        padHisto = (TPad *)c->GetListOfPrimitives()->At(0);
        std::cout << "Prepare for residuals" << endl;
        padHisto->cd();
    }

    const static int ncol = 10;
    const static Color_t cols[ncol] = {kRed, kGreen, kMagenta, kCyan, kBlue, kYellow, kPink, kViolet, kTeal, kSpring};
    const static int nsty = 10;
    const static Style_t stys[nsty] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    if (dataError == "Poisson")
        dataH->plotOn(frame, Name("Data"), DataError(RooAbsData::Poisson), range);
    else
    {
        if (dataError != "SumW2")
            std::cout << "Unknown data error type - assume SumW2" << endl;
        dataH->plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2), range);
    }

    if (dataError == "Poisson")
    {
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
        dataH->plotOn(frame, DataError(RooAbsData::Poisson));
    }
    else
        dataH->plotOn(frame, DataError(RooAbsData::SumW2));

    if (logY)
    {
        frame->SetMinimum(0.2);
        if (padHisto)
            padHisto->SetLogy();
        else
            c->SetLogy();
    }

    printf("Draw main frame \n");
    frame->Draw();
    c->Update();

    if (out != "")
    {
        c->Print(std::string(out + ".C").c_str());
        c->Print(std::string(out + ".pdf").c_str());
    }

    return;
}

void BSA_data_fit::DefineHistos(int selflag)
{

    int bins = 12;

    if (selflag == 1)
        bins = 12;
    if (selflag == 2)
        bins = 7;

    all_phi_cuts_hp = new TH1F("all_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hn = new TH1F("all_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);
    sum = new TH1F("sum", ";#Phi (deg)", bins, 0, 360);
    diff = new TH1F("diff", ";#Phi (deg)", bins, 0, 360);
    ratio = new TH1F("ratio", "A_{LU}(#Phi)", bins, 0, 360);
    ratio_sys = new TH1F("ratio_sys", "A_{LU}(#Phi)", bins, 0, 360);
    ratio_sys_sw = new TH1F("ratio_sys_sw", "A_{LU}(#Phi)", bins, 0, 360);
    ratio_sys_BDT = new TH1F("ratio_sys_BDT", "A_{LU}(#Phi)", bins, 0, 360);

    ratio_sys->SetFillColor(kRed);
    ratio_sys_sw->SetFillColor(kYellow);
    ratio_sys_sw->SetFillStyle(3001);
    ratio_sys_BDT->SetFillColor(kCyan);
    ratio_sys_BDT->SetFillStyle(3001);
    ratio_sys->Sumw2();
    ratio_sys_sw->Sumw2();
    ratio_sys_BDT->Sumw2();

    all_phi_cuts_hp->Sumw2();
    all_phi_cuts_hn->Sumw2();
    sum->Sumw2();
    diff->Sumw2();
    ratio->Sumw2();

    cd_phi_cuts_hp = new TH1F("cd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    cd_phi_cuts_hn = new TH1F("cd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);
    cd_sum = new TH1F("cd_sum", ";#Phi (deg)", bins, 0, 360);
    cd_diff = new TH1F("cd_diff", ";#Phi (deg)", bins, 0, 360);
    cd_ratio = new TH1F("cd_ratio", "A_{LU}(#Phi)", bins, 0, 360);

    cd_phi_cuts_hp->Sumw2();
    cd_phi_cuts_hn->Sumw2();
    cd_sum->Sumw2();
    cd_diff->Sumw2();
    cd_ratio->Sumw2();

    fd_phi_cuts_hp = new TH1F("fd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    fd_phi_cuts_hn = new TH1F("fd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);
    fd_sum = new TH1F("fd_sum", ";#Phi (deg)", bins, 0, 360);
    fd_diff = new TH1F("fd_diff", ";#Phi (deg)", bins, 0, 360);
    fd_ratio = new TH1F("fd_ratio", "A_{LU}(#Phi)", bins, 0, 360);

    fd_phi_cuts_hp->Sumw2();
    fd_phi_cuts_hn->Sumw2();
    fd_sum->Sumw2();
    fd_diff->Sumw2();
    fd_ratio->Sumw2();

    ratio_bkgf = new TH1F("ratio_bkgf", "A_{LU}(#Phi)", bins, 0, 360);
    cd_ratio_bkgf = new TH1F("cd_ratio_bkgf", "A_{LU}(#Phi)", bins, 0, 360);
    fd_ratio_bkgf = new TH1F("fd_ratio_bkgf", "A_{LU}(#Phi)", bins, 0, 360);

    ratio_bkgf->Sumw2();
    cd_ratio_bkgf->Sumw2();
    fd_ratio_bkgf->Sumw2();

    ratio_bkgf->SetLineColor(kRed);
    cd_ratio_bkgf->SetLineColor(kRed);
    fd_ratio_bkgf->SetLineColor(kRed);

    ratio_bkgf->SetMarkerColor(kRed);
    cd_ratio_bkgf->SetMarkerColor(kRed);
    fd_ratio_bkgf->SetMarkerColor(kRed);
    // #############################

    DATA_Pi0_all_phi_cuts_hp = new TH1F("DATA_Pi0_all_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi0_all_phi_cuts_hn = new TH1F("DATA_Pi0_all_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi0_all_phi_cuts_hp->Sumw2();
    DATA_Pi0_all_phi_cuts_hn->Sumw2();

    DATA_Pi0_cd_phi_cuts_hp = new TH1F("DATA_Pi0_cd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi0_cd_phi_cuts_hn = new TH1F("DATA_Pi0_cd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi0_cd_phi_cuts_hp->Sumw2();
    DATA_Pi0_cd_phi_cuts_hn->Sumw2();

    DATA_Pi0_fd_phi_cuts_hp = new TH1F("DATA_Pi0_fd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi0_fd_phi_cuts_hn = new TH1F("DATA_Pi0_fd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi0_fd_phi_cuts_hp->Sumw2();
    DATA_Pi0_fd_phi_cuts_hn->Sumw2();

    // #############################

    MC_Pi0_all_phi_cuts_hp = new TH1F("MC_Pi0_all_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi0_all_phi_cuts_hn = new TH1F("MC_Pi0_all_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi0_all_phi_cuts_hp->SetLineColor(kRed);

    MC_Pi0_all_phi_cuts_hp->Sumw2();
    MC_Pi0_all_phi_cuts_hn->Sumw2();

    MC_Pi0_cd_phi_cuts_hp = new TH1F("MC_Pi0_cd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi0_cd_phi_cuts_hn = new TH1F("MC_Pi0_cd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi0_cd_phi_cuts_hp->Sumw2();
    MC_Pi0_cd_phi_cuts_hn->Sumw2();

    MC_Pi0_fd_phi_cuts_hp = new TH1F("MC_Pi0_fd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi0_fd_phi_cuts_hn = new TH1F("MC_Pi0_fd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi0_fd_phi_cuts_hp->Sumw2();
    MC_Pi0_fd_phi_cuts_hn->Sumw2();

    // #############################

    MC_Pi01g_all_phi_cuts_hp = new TH1F("MC_Pi01g_all_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_all_phi_cuts_hn = new TH1F("MC_Pi01g_all_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_all_phi_cuts_hp->Sumw2();
    MC_Pi01g_all_phi_cuts_hn->Sumw2();

    MC_Pi01g_cd_phi_cuts_hp = new TH1F("MC_Pi01g_cd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_cd_phi_cuts_hn = new TH1F("MC_Pi01g_cd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_cd_phi_cuts_hp->Sumw2();
    MC_Pi01g_cd_phi_cuts_hn->Sumw2();

    MC_Pi01g_fd_phi_cuts_hp = new TH1F("MC_Pi01g_fd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_fd_phi_cuts_hn = new TH1F("MC_Pi01g_fd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_fd_phi_cuts_hp->Sumw2();
    MC_Pi01g_fd_phi_cuts_hn->Sumw2();

    // #############################

    MC_Pi01g_Pi0_all_phi_cuts_hp_ratio = new TH1F("MC_Pi01g_Pi0_all_phi_cuts_hp_ratio", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_Pi0_all_phi_cuts_hn_ratio = new TH1F("MC_Pi01g_Pi0_all_phi_cuts_hn_ratio", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_Pi0_all_phi_cuts_hp_ratio->Sumw2();
    MC_Pi01g_Pi0_all_phi_cuts_hn_ratio->Sumw2();

    MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio = new TH1F("MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio = new TH1F("MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_Pi0_cd_phi_cuts_hp_ratio->Sumw2();
    MC_Pi01g_Pi0_cd_phi_cuts_hn_ratio->Sumw2();

    MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio = new TH1F("MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio", ";#Phi (deg)", bins, 0, 360);
    MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio = new TH1F("MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio", ";#Phi (deg)", bins, 0, 360);

    MC_Pi01g_Pi0_fd_phi_cuts_hp_ratio->Sumw2();
    MC_Pi01g_Pi0_fd_phi_cuts_hn_ratio->Sumw2();

    // #############################

    DATA_Pi01g_all_phi_cuts_hp = new TH1F("DATA_Pi01g_all_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi01g_all_phi_cuts_hn = new TH1F("DATA_Pi01g_all_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi01g_all_phi_cuts_hp->Sumw2();
    DATA_Pi01g_all_phi_cuts_hn->Sumw2();

    DATA_Pi01g_cd_phi_cuts_hp = new TH1F("DATA_Pi01g_cd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi01g_cd_phi_cuts_hn = new TH1F("DATA_Pi01g_cd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi01g_cd_phi_cuts_hp->Sumw2();
    DATA_Pi01g_cd_phi_cuts_hn->Sumw2();

    DATA_Pi01g_fd_phi_cuts_hp = new TH1F("DATA_Pi01g_fd_phi_cuts_hp", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi01g_fd_phi_cuts_hn = new TH1F("DATA_Pi01g_fd_phi_cuts_hn", ";#Phi (deg)", bins, 0, 360);

    DATA_Pi01g_fd_phi_cuts_hp->Sumw2();
    DATA_Pi01g_fd_phi_cuts_hn->Sumw2();

    // #############################

    all_phi_cuts_hp_BKGFree = new TH1F("all_phi_cuts_hp_BKGFree", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hn_BKGFree = new TH1F("all_phi_cuts_hn_BKGFree", ";#Phi (deg)", bins, 0, 360);

    all_phi_cuts_hp_BKGFree->Sumw2();
    all_phi_cuts_hn_BKGFree->Sumw2();

    cd_phi_cuts_hp_BKGFree = new TH1F("cd_phi_cuts_hp_BKGFree", ";#Phi (deg)", bins, 0, 360);
    cd_phi_cuts_hn_BKGFree = new TH1F("cd_phi_cuts_hn_BKGFree", ";#Phi (deg)", bins, 0, 360);

    cd_phi_cuts_hp_BKGFree->Sumw2();
    cd_phi_cuts_hn_BKGFree->Sumw2();

    fd_phi_cuts_hp_BKGFree = new TH1F("fd_phi_cuts_hp_BKGFree", ";#Phi (deg)", bins, 0, 360);
    fd_phi_cuts_hn_BKGFree = new TH1F("fd_phi_cuts_hn_BKGFree", ";#Phi (deg)", bins, 0, 360);

    fd_phi_cuts_hp_BKGFree->Sumw2();
    fd_phi_cuts_hn_BKGFree->Sumw2();

    // #############################
    all_phi_cuts_hp_BKGRatio = new TH1F("all_phi_cuts_hp_BKGRatio", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hn_BKGRatio = new TH1F("all_phi_cuts_hn_BKGRatio", ";#Phi (deg)", bins, 0, 360);

    all_phi_cuts_hp_BKGRatio->Sumw2();
    all_phi_cuts_hn_BKGRatio->Sumw2();

    cd_phi_cuts_hp_BKGRatio = new TH1F("cd_phi_cuts_hp_BKGRatio", ";#Phi (deg)", bins, 0, 360);
    cd_phi_cuts_hn_BKGRatio = new TH1F("cd_phi_cuts_hn_BKGRatio", ";#Phi (deg)", bins, 0, 360);

    cd_phi_cuts_hp_BKGRatio->Sumw2();
    cd_phi_cuts_hn_BKGRatio->Sumw2();

    fd_phi_cuts_hp_BKGRatio = new TH1F("fd_phi_cuts_hp_BKGRatio", ";#Phi (deg)", bins, 0, 360);
    fd_phi_cuts_hn_BKGRatio = new TH1F("fd_phi_cuts_hn_BKGRatio", ";#Phi (deg)", bins, 0, 360);

    fd_phi_cuts_hp_BKGRatio->Sumw2();
    fd_phi_cuts_hn_BKGRatio->Sumw2();

    // ############################

    all_phi_cuts_hp_noP = new TH1F("all_phi_cuts_hp_noP", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hn_noP = new TH1F("all_phi_cuts_hn_noP", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi0_all_phi_cuts_hp_noP = new TH1F("DATA_Pi0_all_phi_cuts_hp_noP", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi0_all_phi_cuts_hn_noP = new TH1F("DATA_Pi0_all_phi_cuts_hn_noP", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hp_BKGFree_noP = new TH1F("all_phi_cuts_hp_BKGFree_noP", ";#Phi (deg)", bins, 0, 360);
    all_phi_cuts_hn_BKGFree_noP = new TH1F("all_phi_cuts_hn_BKGFree_noP", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi01g_all_phi_cuts_hp_noP = new TH1F("DATA_Pi01g_all_phi_cuts_hp_noP", ";#Phi (deg)", bins, 0, 360);
    DATA_Pi01g_all_phi_cuts_hn_noP = new TH1F("DATA_Pi01g_all_phi_cuts_hn_noP", ";#Phi (deg)", bins, 0, 360);

    all_phi_cuts_hp_noP->Sumw2();
    all_phi_cuts_hn_noP->Sumw2();
    DATA_Pi0_all_phi_cuts_hp_noP->Sumw2();
    DATA_Pi0_all_phi_cuts_hn_noP->Sumw2();
    all_phi_cuts_hp_BKGFree_noP->Sumw2();
    all_phi_cuts_hn_BKGFree_noP->Sumw2();
    DATA_Pi01g_all_phi_cuts_hp_noP->Sumw2();
    DATA_Pi01g_all_phi_cuts_hn_noP->Sumw2();

    return;
}

TCut BSA_data_fit::ALL_cuts(int selFlag)
{

    std::string c;

    if (selFlag == 1)
        // c = " _strip_Ph_P > 2 && _strip_Q2 > 1.0 && _strip_W > 2 && _strip_Nuc_P > 0.3 && _strip_El_P > 1.0 && TMath::Abs(_Phi_Nuc - _Phi_Ph) < 2 &&  TMath::Abs(_mm2_eNg_N) < 1 && _theta_gamma_X < 3  && _strip_El_vz < 10 && _strip_El_vz > -12 && TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal) < 0.75";;

        c = "_Exclusive_4DChi2 < 8 && _strip_Ph_P > 2 && _strip_Q2 > 1.0 && _strip_W > 2 && _strip_Nuc_P > 0.3 && _strip_El_P > 1.0 && TMath::Abs(_Phi_Nuc - _Phi_Ph) < 2 && TMath::Abs(_t_Nuc - _t_Ph) < 0.25 && TMath::Abs(_mm2_eNg_N) < 1 && _theta_gamma_X < 3  && _strip_El_vz < 2 && _strip_El_vz > -8 && TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal) < 0.75 && _t_Nuc > -1.9 && (_strip_El_vz-_strip_Nuc_vz) < 4.445 && (_strip_El_vz-_strip_Nuc_vz) > -4.055";

    //_Exclusive_4DChi2 < 8 && _strip_Ph_P > 2 && _strip_Q2 > 1.0 && _strip_W > 2 && _strip_Nuc_P > 0.3 && _strip_El_P > 1.0 && TMath::Abs(_Phi_Nuc - _Phi_Ph) < 2 && TMath::Abs(_t_Nuc - _t_Ph) < 0.25 && TMath::Abs(_mm2_eNg_N) < 1 && _theta_gamma_X < 3 && _theta_gamma_e > 5 && _strip_El_vz < 10 && _strip_El_vz > -12 && TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal) < 0.75)

    if (selFlag == 2)
        c = "(_Exclusive_4DChi2<8 && _theta_N_e> 5 && _theta_gamma_e>5 &&_strip_Q2 > 1.0 &&_strip_W > 2 &&_strip_Ph_P > 2 &&_strip_Nuc_P > 0.35  &&_strip_El_P > 1.0 &&_strip_Nuc_Theta < 150 &&_strip_El_Theta > 5.5 &&_strip_El_vz < 2 &&_strip_El_vz > -8 && TMath::Abs(_Phi_Nuc - _Phi_Ph) < 2 && TMath::Abs(_t_Nuc - _t_Ph) < 0.5 &&TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal) < 0.75 &&_mm2_eNg < 2.5 &&_mm2_eNg > 0. &&_mm2_eNg_N > -1 &&_mm2_eNg_N < 1 &&_theta_gamma_X < 3 )  && _t_Nuc > -1.9  && (_strip_El_vz-_strip_Nuc_vz) < 4.445 && (_strip_El_vz-_strip_Nuc_vz) > -4.055";

    //(_Exclusive_4DChi2<8 &&_strip_Q2 > 1.0 &&_strip_W > 2 &&_strip_Ph_P > 2 &&_strip_Nuc_P > 0.35 &&_strip_El_P > 1.0 &&_strip_El_Theta > 5.5 &&_strip_El_vz < 2 &&_strip_El_vz > -8 && TMath::Abs(_Phi_Nuc - _Phi_Ph) < 0.8 && TMath::Abs(_t_Nuc - _t_Ph) < 0.5 &&TMath::Sqrt(_Xbal * _Xbal + _Ybal * _Ybal + _Zbal * _Zbal) < 0.75 &&_mm2_eNg < 2.5 &&_mm2_eNg > 0. &&_mm2_eNg_N > -0.2 &&_mm2_eNg_N < 0.2 &&_theta_gamma_X < 3 )

    // (Exclusive_4DChi2<5 && theta_N_e> 5 && theta_gamma_e>5 && bestCandidateFlag &&strip_Q2 > 1.0 &&strip_W > 2 &&strip_Ph_P > 2 &&strip_Nuc_P > 0.35 &&strip_El_P > 1.0 &&strip_Nuc_Theta < 150 &&strip_El_Theta > 5.5 &&strip_El_vz < 10 &&strip_El_vz > -12 && TMath::Abs(Phi_Nuc - Phi_Ph) < 0.8 && TMath::Abs(t_Nuc - t_Ph) < 0.5 &&TMath::Sqrt(Xbal * Xbal + Ybal * Ybal + Zbal * Zbal) < 0.75 &&mm2_eNg < 2.5 &&mm2_eNg > 0. &&mm2_eNg_N > -0.2 &&mm2_eNg_N < 0.2 &&theta_gamma_X < 3 )

    //(Exclusive_4DChi2<2.5 && theta_N_e> 5 && theta_gamma_e>5 && bestCandidateFlag &&strip_Q2 > 1.0 &&strip_W > 2 &&strip_Ph_P > 2 &&strip_Nuc_P > 0.35 &&strip_El_P > 1.0 &&strip_Nuc_Theta < 150 &&strip_El_Theta > 5.5 &&strip_El_vz < 10 &&strip_El_vz > -12 && TMath::Abs(Phi_Nuc - Phi_Ph) < 0.7 && TMath::Abs(t_Nuc - t_Ph) < 0.3 &&TMath::Sqrt(Xbal * Xbal + Ybal * Ybal + Zbal * Zbal) < 0.75 && mm2_eNg < 2.5 && mm2_eNg > 0. && mm2_eNg_N > -0.2 &&mm2_eNg_N < 0.2 && strip_Nuc_P < 1.1 && theta_gamma_X < 3 && ((strip_Nuc_status >=2000 && strip_Nuc_status < 4000 && strip_Nuc_Theta<35) || (strip_Nuc_status >=4000 && strip_Nuc_Theta>35 && !((strip_Nuc_Phi<-70 && strip_Nuc_Phi>-90)||(strip_Nuc_Phi>25 && strip_Nuc_Phi<50)||(strip_Nuc_Phi>150 && strip_Nuc_Phi<170))))) && ((strip_Nuc_status >=4000 && strip_Ph_status>=1000  && strip_Ph_status<2000) || (strip_Nuc_status >=4000 && strip_Ph_status>=2000) || (strip_Nuc_status >=2000 && strip_Nuc_status < 4000 && strip_Ph_status>=1000  && strip_Ph_status<2000) || (strip_Nuc_status >=2000 && strip_Nuc_status < 4000 && strip_Ph_status>=2000))

    if (selFlag == 4)
        c = "_best4DChi2Flag  && abs(_mm2_Pi0_X)<0.7 && _strip_Pi0_IM2>0 && _strip_Pi0_IM2<0.035 && _strip_Pi0_2DChi2<0.12 && _Exclusive_4DChi2<8  && _strip_Q2>1.0 && _strip_W>2 && _strip_Ph1_P>0.35 && _strip_Ph2_P>0.35 && _strip_Nuc_P>0.3 && _strip_El_P>1.0 && TMath::Abs(_Phi_Nuc-_Phi_Pi0)<2 && TMath::Abs(_t_Nuc-_t_Pi0)<2 && TMath::Abs(_mm2_eNgg_N)<1 && _theta_Pi0_X<5 && _theta_Pi0_e>6 &&_strip_El_vz < 2 &&_strip_El_vz > -8 && TMath::Sqrt(_Xbal*_Xbal+_Ybal*_Ybal+_Zbal*_Zbal)<0.75 && _t_Nuc > -1.9 && (_strip_El_vz-_strip_Nuc_vz) < 4.445 && (_strip_El_vz-_strip_Nuc_vz) > -4.055";

    // c = "best4DChi2Flag  && abs(mm2_Pi0_X)<0.7 && strip_Pi0_IM2>0 && strip_Pi0_IM2<0.035 && strip_Pi0_2DChi2<0.12 && Exclusive_4DChi2<8  && strip_Q2>1.0 && strip_W>2 && strip_Ph1_P>0.35 && strip_Ph2_P>0.35 && strip_Nuc_P>0.3 && strip_El_P>1.0 && TMath::Abs(Phi_Nuc-Phi_Pi0)<2 && TMath::Abs(t_Nuc-t_Pi0)<2 && TMath::Abs(mm2_eNgg_N)<1 && theta_Pi0_X<5 && theta_Pi0_e>6 &&strip_El_vz < 10 &&strip_El_vz > -12 && TMath::Sqrt(Xbal*Xbal+Ybal*Ybal+Zbal*Zbal)<0.75";

    if (selFlag == 8)
        c = "_best4DChi2Flag  && abs(_mm2_Pi0_X)<0.7 && _strip_Pi0_IM2>0 && _strip_Pi0_IM2<0.035 && _strip_Pi0_2DChi2<0.12 && _Exclusive_4DChi2<2.5 && _strip_Q2>1.0 && _strip_W>2 && _strip_Ph1_P>0.35 && _strip_Ph2_P>0.35 && _strip_Nuc_P>0.35 && _strip_El_P>1.0 && TMath::Abs(_Phi_Nuc-_Phi_Pi0)<2 && TMath::Abs(_t_Nuc-_t_Pi0)<2 && TMath::Abs(_mm2_eNgg_N)<1 && _theta_N_e> 5 && _theta_Pi0_X<5 && _theta_Pi0_e>6 &&_strip_El_vz < 2 &&_strip_El_vz > -8 && TMath::Sqrt(_Xbal*_Xbal+_Ybal*_Ybal+_Zbal*_Zbal)<0.75 && _t_Nuc > -1.9 && (_strip_El_vz-_strip_Nuc_vz) < 4.445 && (_strip_El_vz-_strip_Nuc_vz) > -4.055";
    // && ((_strip_Nuc_status >=2000 && _strip_Nuc_status < 4000 && _strip_Nuc_Theta<35) || (_strip_Nuc_status >=4000 && _strip_Nuc_Theta>35 && !((_strip_Nuc_Phi<-70 && _strip_Nuc_Phi>-90)||(_strip_Nuc_Phi>25 && _strip_Nuc_Phi<50)||(_strip_Nuc_Phi>150 && _strip_Nuc_Phi<170))))";

    // c = "best4DChi2Flag  && abs(mm2_Pi0_X)<0.7 && strip_Pi0_IM2>0 && strip_Pi0_IM2<0.035 && strip_Pi0_2DChi2<0.12 && Exclusive_4DChi2<8  && strip_Q2>1.0 && strip_W>2 && strip_Ph1_P>0.35 && strip_Ph2_P>0.35 && strip_Nuc_P>0.3 && strip_El_P>1.0 && TMath::Abs(Phi_Nuc-Phi_Pi0)<2 && TMath::Abs(t_Nuc-t_Pi0)<2 && TMath::Abs(mm2_eNgg_N)<1 && theta_N_e> 5 && theta_Pi0_X<5 && theta_Pi0_e>6 &&strip_El_vz < 10 &&strip_El_vz > -12 && TMath::Sqrt(Xbal*Xbal+Ybal*Ybal+Zbal*Zbal)<0.75";

    return TCut(c.c_str());
}

//------- Main selection ------------------------------------------------
TCut BSA_data_fit::processSelection(int selFlag, std::string DetectorRegion = "", std::string selCriteria = "")
{

    // should add error message in case of unknown selection criteria

    TCut TheCut = TCut("");

    if (DetectorRegion == "ALL")
    {

        TheCut += ALL_cuts(selFlag);
        // TheCut += TCut("_RunNumber<11286");
    }
    if (DetectorRegion == "FD")
    {
        TheCut += ALL_cuts(selFlag);
        TheCut += TCut("_strip_Nuc_status >=2000 && _strip_Nuc_status < 4000 ");
        // TheCut += TCut("_RunNumber<11286");
    }
    if (DetectorRegion == "CD")
    {
        TheCut += ALL_cuts(selFlag);
        TheCut += TCut("_strip_Nuc_status >=4000 ");
        // TheCut += TCut("_RunNumber<11286");
    }

    return TheCut;
}

#endif
