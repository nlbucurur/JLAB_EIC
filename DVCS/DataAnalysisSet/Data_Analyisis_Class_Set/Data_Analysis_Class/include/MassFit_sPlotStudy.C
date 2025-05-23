#define MassFit_sPlotStudy_cxx
#include "MassFit_sPlotStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassFit_sPlotStudy::MassFit(TString S_cut_bin)//, bool Pmissbinned = false, double Pmin = 0, double Pmax = 10)
{
  int mask = 1;
  using namespace RooFit;
  using namespace RooStats;

  bool binned = true;

  bool output = true; // produce output
  bool minos = false;
  bool logY = false; // log-scale
  int nCPU = 4;
  bool w2Error = true; // MUST BE TRUE - EXCEPT WHEN SCALING ERRORS TO DATA STATs

  gROOT->Reset();

  std::string pref = ""; // output prefix
  std::string suff = ""; // output suffix

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

  //----- TFile & TTree
  TString name0="~/Quality2f_Data_P.root";
  TFile *f1 = new TFile(name0,"READ");
  TTree *tree0 = (mask & 1) ? (TTree *)f1->Get("pDVCS") : (TTree *)NULL;      // expected

  //Limits on the RooRealVar must be applied as cuts in the data
  TCut cut0 = TCut("bestCandidateFlag==1 && strip_Q2 > 1.0 && strip_W > 2 && strip_Nuc_P > 0.35 && strip_El_P > 1.0 && strip_Ph_P>2  && strip_El_vz < 10 && strip_El_vz > -12 && TMath::Abs(Phi_Nuc - Phi_Ph) < 2 && TMath::Abs(t_Nuc - t_Ph) < 2 && TMath::Sqrt(Xbal * Xbal + Ybal*Ybal + Zbal*Zbal) <1 && mm2_eg > 0 && mm2_eg < 2.5 && abs(mm2_eNg)<0.1");
  TCut cut_bin =TCut(S_cut_bin);
  TCut cut = cut0 + cut_bin; 
  TFile ofile("sw_Data.root","RECREATE");
  TTree *tree1 = tree0->CopyTree(cut); //Not working for some reason
    
  //----- define TCanvas
  TCanvas *c1 = NULL;

  string channel_ALL = "";

  if (mask & 1)
    {
      channel_ALL += "NDE";
      c1 = new TCanvas("c1", "");
      prepareCanvas(c1, 1);
    }

  //If Trees are not aligned, these branches will be written together with the sweights
  RooRealVar mm2_eg("mm2_eg", "mm2_eg", 0.0, 2.5);
  RooRealVar Phi_Ph("Phi_Ph", "Phi_Ph", 0.0, 360);
  RooRealVar delta_Phi("delta_Phi", "delta_Phi", -1000., 1000.);
  RooRealVar delta_t("delta_t", "delta_t", -1000., 1000.);

  RooArgSet parList_(mm2_eg);
  parList_.add(delta_Phi);
  parList_.add(delta_t);
  parList_.add(Phi_Ph);

    
  RooDataSet *dataS1 = NULL;
  RooDataHist *dataH1 = NULL;

  if (mask & 1)
    {
      dataS1 = new RooDataSet("dataS1", "pDVCS", tree1, parList_);
      dataH1 = new RooDataHist("dataH1", "dataH1", RooArgSet(mm2_eg), *dataS1);
    }

  // Construct signal pdf
  RooRealVar Syield("Syield", "Syield", 500000, 0, 2000000);
  RooRealVar B1yield("B1yield", "B1yield", 500000, 0, 2000000);
  RooRealVar B2yield("B2yield", "B2yield", 500000, 0, 2000000);
  RooRealVar B3yield("B3yield", "B3yield", 500000, 0, 2000000);

    
  //Double sided Crystal ball DVCS
  RooRealVar mean("mean", "mean", 0.88);    
  RooRealVar sigma("sigma", "sigma", 0.1,0,2.0); 
  RooRealVar alpha("alpha", "alpha", 1.3307, 0.0 ,100.0);
  RooRealVar n("n", "n", 2.,1.,100.);
  RooRealVar alpha1("alpha1", "alpha1", -1.73,-100.,-1.);
  RooRealVar n1("n1", "n1", 1.487,1.,100.);
  RooCrystalBall ball("ball","DoubleSidedCB",mm2_eg,mean,sigma,alpha,n,alpha1,n1);

  //Double sided Crystal ball Pi0 
  RooRealVar mean_0("mean_0", "mean_0", 1.524,1.0,2.0);    
  RooRealVar sigma_0("sigma_0", "sigma_0", 0.4, 0.,1.0); 
  RooRealVar alpha_0("alpha_0", "alpha_0", 1.103, 0.0 ,200.0);
  RooRealVar n_0("n_0", "n_0", 100, 1., 200.);
  RooRealVar alpha1_0("alpha1_0", "alpha1_0", -2.6523,-200.0 ,0.0);
  RooRealVar n1_0("n1_0", "n1_0", 97, 1., 200.);
  RooCrystalBall ball0("ball0","DoubleSidedCB",mm2_eg,mean_0,sigma_0,alpha_0,n_0,alpha1_0,n1_0);

  // gaussian Pi0 
  RooRealVar gaus_mean0("gaus_mean0", "gaus_mean0", 1.34, 0., 1.4);
  RooRealVar gaus_sigma0("gaus_sigma0", "gaus_sigma0", 0.1549,0,1.0);
  RooGaussian gaus0("gaus0", "gaus0", mm2_eg, gaus_mean0, gaus_sigma0);

  // gaussian Pi0 - fit
  RooRealVar gaus_mean("gaus_mean", "gaus_mean", 1.52);//, 0.8, 2.0);
  RooRealVar gaus_sigma("gaus_sigma", "gaus_sigma", 0.4549);//,0,1.0);
  RooGaussian gaus("gaus", "gaus", mm2_eg, gaus_mean, gaus_sigma);

  //Double sided Crystal ball Pi0 - fit Pi0_Data training
  RooRealVar mean_1("mean_1", "mean_1", 1.7524);    
  RooRealVar sigma_1("sigma_1", "sigma_1", 0.38); 
  RooRealVar alpha_1("alpha_1", "alpha_1", 1.103);
  RooRealVar n_1("n_1", "n_1", 100);
  RooRealVar alpha1_1("alpha1_1", "alpha1_1", -2.6523);
  RooRealVar n1_1("n1_1", "n1_1", 97);
  RooCrystalBall ball1("ball1","DoubleSidedCB",mm2_eg,mean_1,sigma_1,alpha_1,n_1,alpha1_1,n1_1);

  //Double sided Crystal ball Pi0 - fit Pi0_Data training
  RooRealVar mean_3("mean_3", "mean_3", 1.21815);    
  RooRealVar sigma_3("sigma_3", "sigma_3", 0.306); 
  RooRealVar alpha_3("alpha_3", "alpha_3", 2.022);
  RooRealVar n_3("n_3", "n_3", 43.9);
  RooRealVar alpha1_3("alpha1_3", "alpha1_3", -100);
  RooRealVar n1_3("n1_3", "n1_3", 97);
  RooCrystalBall ball3("ball3","DoubleSidedCB",mm2_eg,mean_3,sigma_3,alpha_3,n_3,alpha1_3,n1_3);


  //Double sided Crystal ball Pi0 - fit Pi0_Data training
  RooRealVar mean_4("mean_4", "mean_4", 1.8413);    
  RooRealVar sigma_4("sigma_4", "sigma_4", 0.3329); 
  RooRealVar alpha_4("alpha_4", "alpha_4", 21);
  RooRealVar n_4("n_4", "n_4", 69);
  RooRealVar alpha1_4("alpha1_4", "alpha1_4", -3.6835);
  RooRealVar n1_4("n1_4", "n1_4", 69);
  RooCrystalBall ball4("ball4","DoubleSidedCB",mm2_eg,mean_4,sigma_4,alpha_4,n_4,alpha1_4,n1_4);

  //Pi0 as sum of 3 gaussians
  // gaussian Pi0
  RooRealVar gaus_mean1("gaus_mean1", "gaus_mean1", 1.0677);
  RooRealVar gaus_sigma1("gaus_sigma1", "gaus_sigma1", 0.2072);
  RooGaussian gaus1("gaus1", "gaus1", mm2_eg, gaus_mean1, gaus_sigma1);

  // gaussian Pi0
  RooRealVar gaus_mean2("gaus_mean2", "gaus_mean2", 1.5086);
  RooRealVar gaus_sigma2("gaus_sigma2", "gaus_sigma2", 0.41975);
  RooGaussian gaus2("gaus2", "gaus2", mm2_eg, gaus_mean2, gaus_sigma2);

  // gaussian Pi0
  RooRealVar gaus_mean3("gaus_mean3", "gaus_mean3", 2.1077);
  RooRealVar gaus_sigma3("gaus_sigma3", "gaus_sigma3", 0.2229);
  RooGaussian gaus3("gaus3", "gaus3", mm2_eg, gaus_mean3, gaus_sigma3);

  // to fit complete model
  RooAddPdf pdf_gausswithgauss("pdf_gausswithgauss", "pdf_gausswithgauss", RooArgList(ball, ball0, gaus0), RooArgList(Syield, B1yield, B2yield));

  RooArgList FitparList(mm2_eg);
  FitparList.add(mean_0);
  FitparList.add(sigma_0);
  FitparList.add(alpha_0);
  FitparList.add(n_0);
  FitparList.add(alpha1_0);
  FitparList.add(n1_0);
  FitparList.add(mean);
  FitparList.add(sigma);
  FitparList.add(alpha);
  FitparList.add(n);
  FitparList.add(alpha1);
  FitparList.add(n1);
  FitparList.add(mean_1);
  FitparList.add(sigma_1);
  FitparList.add(alpha_1);
  FitparList.add(n_1);
  FitparList.add(alpha1_1);
  FitparList.add(n1_1);
  FitparList.add(gaus_mean);
  FitparList.add(gaus_sigma);
  FitparList.add(gaus_mean0);
  FitparList.add(gaus_sigma0);
  FitparList.add(gaus_mean1);
  FitparList.add(gaus_sigma1);
  FitparList.add(gaus_mean2);
  FitparList.add(gaus_sigma2);
  FitparList.add(gaus_mean3);
  FitparList.add(gaus_sigma3);
  FitparList.add(Syield);
  FitparList.add(B1yield);
  FitparList.add(B2yield);
  FitparList.add(B3yield);
  RooFitResult *result1;

  if (mask & 1)
    {

      RooPlot *frame1 = mm2_eg.frame(Title(" "));

      result1 = pdf_gausswithgauss.fitTo(*dataH1, Extended(kTRUE), Save(), myMinos, NumCPU(nCPU), SumW2Error(w2Error));

      rooDisplayFit(mm2_eg, dataH1, frame1, c1, pdf_gausswithgauss, result1, "MassFit", FitparList, false, "SumW2", RooCmdArg(), RooFit::Layout(0.55, 0.8, 0.9));
      int sploting = sPlotter(dataS1, tree1, pdf_gausswithgauss, RooArgList(Syield, B1yield, B2yield), 0);
      ofile.Close();
    }


}
