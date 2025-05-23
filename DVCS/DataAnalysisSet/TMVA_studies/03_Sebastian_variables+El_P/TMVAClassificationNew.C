/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set of classifiers is used.
/// The output file "TMVAC.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// Launch the GUI via the command:
///
///     root -l ./TMVAGui.C
///
/// You can also compile and run the example with the following commands
///
///     make
///     ./TMVAClassification <Methods>
///
/// where: `<Methods> = "method1 method2"` are the TMVA classifier names
/// example:
///
///     ./TMVAClassification Fisher LikelihoodPCA BDT
///
/// If no method given, a default set is of classifiers is used
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAClassification
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void CalculateMissingMomentum(TTree *tree, const TString &output_filename)
{
   // Set up 4-vectors
   TLorentzVector ElectronBeam, NucTarget_Vec, Target_Vec;
   TLorentzVector Ph_Vec, Nuc_Vec, El_Vec, Pmiss, Pmiss_Nuc;

   // Beam and target setup
   Int_t RunNumber;
   double Ebeam = 10.2;
   tree->SetBranchAddress("_RunNumber", &RunNumber);

   if (RunNumber >= 6420)
      Ebeam = 10.2;

   if (RunNumber > 10000)
      Ebeam = 10.4;
   double Pmass = 0.938272;
   double Nmass = 0.9395654;
   double Dmass = 1.8756;

   ElectronBeam.SetXYZT(0, 0, Ebeam, Ebeam);
   Target_Vec.SetXYZT(0, 0, 0, Dmass);
   NucTarget_Vec.SetXYZT(0, 0, 0, Pmass);

   // Variables for branch addresses
   static double _strip_El_px, _strip_El_py, _strip_El_pz, _strip_El_E;
   static double _strip_Nuc_px, _strip_Nuc_py, _strip_Nuc_pz, _strip_Nuc_E;
   static double _strip_Ph_px, _strip_Ph_py, _strip_Ph_pz, _strip_Ph_E;
   static double Pmiss_mag, Pmiss_Nuc_mag, Pmiss_perp, Pmiss_Nuc_perp;

   // Set branch addresses
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

   // Create new branches
   TFile *output_file = new TFile(output_filename, "RECREATE");
   TTree *new_tree = tree->CloneTree(-1);

   TBranch *branch_Pmiss_mag = new_tree->Branch("_Pmiss_mag", &Pmiss_mag, "_Pmiss_mag/D");
   TBranch *branch_Pmiss_Nuc_mag = new_tree->Branch("_Pmiss_Nuc_mag", &Pmiss_Nuc_mag, "_Pmiss_Nuc_mag/D");
   TBranch *branch_Pmiss_perp = new_tree->Branch("_Pmiss_perp", &Pmiss_perp, "_Pmiss_perp/D");
   TBranch *branch_Pmiss_Nuc_perp = new_tree->Branch("_Pmiss_Nuc_perp", &Pmiss_Nuc_perp, "_Pmiss_Nuc_perp/D");

   // Process each event
   for (int i = 0; i < tree->GetEntries(); i++)
   {
      tree->GetEntry(i);

      // Set particle 4-vectors
      Ph_Vec.SetPxPyPzE(_strip_Ph_px, _strip_Ph_py, _strip_Ph_pz, _strip_Ph_E);
      Nuc_Vec.SetPxPyPzE(_strip_Nuc_px, _strip_Nuc_py, _strip_Nuc_pz, _strip_Nuc_E);
      El_Vec.SetPxPyPzE(_strip_El_px, _strip_El_py, _strip_El_pz, _strip_El_E);

      // Calculate missing momentum
      Pmiss = ElectronBeam + Target_Vec - El_Vec - Ph_Vec - Nuc_Vec;
      Pmiss_Nuc = ElectronBeam + NucTarget_Vec - El_Vec - Ph_Vec - Nuc_Vec;

      // Magnitude and transverse components
      Pmiss_mag = Pmiss.P();
      Pmiss_Nuc_mag = Pmiss_Nuc.P();
      Pmiss_perp = Pmiss.Perp();
      Pmiss_Nuc_perp = Pmiss_Nuc.Perp();

      // Fill the branches
      branch_Pmiss_mag->Fill();
      branch_Pmiss_Nuc_mag->Fill();
      branch_Pmiss_perp->Fill();
      branch_Pmiss_Nuc_perp->Fill();
   }

   // Save the updated tree
   new_tree->Write();
   output_file->Close();
   delete output_file;
}

int TMVAClassificationNew(TString myMethodList = "")
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // Methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string, int> Use;

   // Cut optimisation
   Use["Cuts"] = 1;
   Use["CutsD"] = 1;
   Use["CutsPCA"] = 0;
   Use["CutsGA"] = 0;
   Use["CutsSA"] = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"] = 1;
   Use["LikelihoodD"] = 0;   // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"] = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"] = 0;
   Use["LikelihoodMIX"] = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"] = 1;
   Use["PDERSD"] = 0;
   Use["PDERSPCA"] = 0;
   // Use["PDEFoam"] = 1;
   // Use["PDEFoamBoost"] = 0; // uses generalised MVA method boosting
   Use["KNN"] = 1; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"] = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"] = 0;
   Use["FisherG"] = 0;
   Use["BoostedFisher"] = 0; // uses generalised MVA method boosting
   Use["HMatrix"] = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"] = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"] = 0;
   Use["FDA_MC"] = 0;
   Use["FDA_MT"] = 0;
   Use["FDA_GAMT"] = 0;
   Use["FDA_MCMT"] = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"] = 0;      // Recommended ANN
   Use["MLPBFGS"] = 0;  // Recommended ANN with optional training method
   Use["MLPBNN"] = 1;   // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"] = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"] = 0;  // ROOT's own ANN
#ifdef R__HAS_TMVAGPU
   Use["DNN_GPU"] = 1; // CUDA-accelerated DNN training.
#else
   Use["DNN_GPU"] = 0;
#endif

#ifdef R__HAS_TMVACPU
   Use["DNN_CPU"] = 1; // Multi-core accelerated DNN.
#else
   Use["DNN_CPU"] = 0;
#endif
   //
   // Support Vector Machine
   Use["SVM"] = 1;
   //
   // Boosted Decision Trees
   Use["BDT"] = 1;  // uses Adaptive Boost
   Use["BDTG"] = 0; // uses Gradient Boost
   Use["BDTB"] = 0; // uses Bagging
   Use["BDTD"] = 0; // decorrelation + Adaptive Boost
   Use["BDTF"] = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"] = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "")
   {
      for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++)
         it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');
      for (UInt_t i = 0; i < mlist.size(); i++)
      {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end())
         {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++)
               std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   /**
    * ## Signal Data
    * This is refered as the signal data
    */
   TFile *input_s(0);
   TFile *input_b(0);
   // TString fname = "./tmva_class_example.root";
   // TString fname_s = "../data/1pDVCS_simulation.root";
   // TString fname_b = "../data/0pDVCS_Pi0dataAsDVCS_10p2.root";

   TString fname_s = "../../data/1pDVCS_simulation_withPmiss.root";
   TString fname_b = "../../data/0pDVCS_Pi0dataAsDVCS_10p2_withPmiss.root";
   if (!gSystem->AccessPathName(fname_s) && !gSystem->AccessPathName(fname_b))
   {
      input_s = TFile::Open(fname_s); // check if file in local directory exists
      input_b = TFile::Open(fname_b); // check if file in local directory exists
   }
   else
   {
      TFile::SetCacheFileDir(".");
      input_s = TFile::Open("https://root.cern/files/tmva_class_example.root", "CACHEREAD");
   }
   if (!input_s || !input_b)
   {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassification       : Using input file: " << input_s->GetName() << std::endl;

   // Register the training and test trees

   TTree *signalTree = (TTree *)input_s->Get("pDVCS_stripped");
   TTree *background = (TTree *)input_b->Get("pDVCS_stripped");

   // Calculate missing momentum and save to new files
   // CalculateMissingMomentum(signalTree, "../data/1pDVCS_simulation_withPmiss.root");
   // CalculateMissingMomentum(background, "../data/0pDVCS_Pi0dataAsDVCS_10p2_withPmiss.root");

   // Now reopen the updated trees
   // input_s->Close();
   // input_b->Close();

   // input_s = TFile::Open("../data/1pDVCS_simulation_withPmiss.root");
   // input_b = TFile::Open("../data/0pDVCS_Pi0dataAsDVCS_10p2_withPmiss.root");

   // signalTree->ResetBranchAddresses();
   // background->ResetBranchAddresses();

   // Lorena's Clock tick duration 2.09e-10 seconds

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName("./output_files/TMVAC.root");
   // TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TFile *outputFile = new TFile(outfileName, "RECREATE");

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   //    dataloader->AddVariable( "myvar1 := var1+var2", 'F' );
   //    dataloader->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
   //    dataloader->AddVariable( "var3",                "Variable 3", "units", 'F' );
   //    dataloader->AddVariable( "var4",                "Variable 4", "units", 'F' );

   dataloader->AddVariable("_mm2_eNg", "mm2_eNg missing Neutron", "GeV", 'D');
   // dataloader->AddVariable("_mm2_eNg_N", "mm2_eNg_N All detected", "GeV", 'D');
   // dataloader->AddVariable("_mm2_eNX_N", "mm2_eNX_N missing photon", "GeV", 'D');
   dataloader->AddVariable("_mm2_eg", "mm2_eg missing proton", "GeV", 'D');
   dataloader->AddVariable("_delta_t", "Delta t", "GeV", 'D');
   dataloader->AddVariable("_delta_Phi", "Delta Phi", "°", 'D');
   dataloader->AddVariable("_theta_gamma_X", "Angle missing photon - photon out", "°", 'D');
   // dataloader->AddVariable("_Pmiss_mag", "Magnitude of Missing Momentum", "GeV", 'D');
   // dataloader->AddVariable("_Pmiss_Nuc_mag", "Magnitude of Missing Nucleon Momentum", "GeV", 'D');
   // dataloader->AddVariable("_Pmiss_perp", "Transverse Missing Momentum", "GeV", 'D');
   // dataloader->AddVariable("_Pmiss_Nuc_perp", "Transverse Missing Nucleon Momentum", "GeV", 'D');
    dataloader->AddVariable("_strip_Ph_P", "photon Momentum", "GeV", 'D');
    dataloader->AddVariable("_strip_El_P", "electron Momentum", "GeV", 'D');

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   //    dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //    dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree(signalTree, signalWeight);
   dataloader->AddBackgroundTree(background, backgroundWeight);

   // Apply additional cuts on the signal and background samples (can be different)
   // TCut mycuts = "abs(_delta_Phi)<5 && _theta_gamma_e > 6 && _delta_t >= -0.46292 && _delta_t <= 0.47175 && _strip_El_chi2pid >= -4.56920 && _strip_El_chi2pid <= 3.61976 && _strip_Nuc_chi2pid >= -195.04711 && _strip_Nuc_chi2pid <= 201.30658 && _mm2_eNg >= -0.37894 && _mm2_eNg <= 2.42267 && _mm2_eNg_N >= -0.19478 && _mm2_eNg_N <= 0.15635 && _mm2_eNX_N >= -3.95236 && _mm2_eNX_N <= 3.74568 && _mm2_eg >= -0.12854 && _mm2_eg <= 2.21362 && _Pmiss_mag < 1.0 && _Pmiss_mag > 0.0 && _Pmiss_Nuc_mag < 1.0 && _Pmiss_Nuc_mag > 0.0 && _Pmiss_perp < 1.0 && _Pmiss_perp > 0.0 && _Pmiss_Nuc_perp < 1.0 && _Pmiss_Nuc_perp > 0.0";  // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   // TCut mycutb = "abs(_delta_Phi)<5  && _theta_gamma_e > 6 && _delta_t >= -0.46292 && _delta_t <= 0.47175 && _strip_El_chi2pid >= -4.56920 && _strip_El_chi2pid <= 3.61976 && _strip_Nuc_chi2pid >= -195.04711 && _strip_Nuc_chi2pid <= 201.30658 && _mm2_eNg >= -0.37894 && _mm2_eNg <= 2.42267 && _mm2_eNg_N >= -0.19478 && _mm2_eNg_N <= 0.15635 && _mm2_eNX_N >= -3.95236 && _mm2_eNX_N <= 3.74568 && _mm2_eg >= -0.12854 && _mm2_eg <= 2.21362 && _Pmiss_mag < 1.0 && _Pmiss_mag > 0.0 && _Pmiss_Nuc_mag < 1.0 && _Pmiss_Nuc_mag > 0.0 && _Pmiss_perp < 1.0 && _Pmiss_perp > 0.0 && _Pmiss_Nuc_perp < 1.0 && _Pmiss_Nuc_perp > 0.0"; // for example: TCut mycutb = "abs(var1)<0.5";
   TCut mycuts = "abs(_delta_Phi)<5";  // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = "abs(_delta_Phi)<5"; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the dataloader how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   //
   // To also specify the number of testing events, use:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut,
   //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   dataloader->PrepareTrainingAndTestTree(mycuts, mycutb,
                                          "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V");

   // ### Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/old_site/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

   if (Use["CutsD"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsD",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

   if (Use["CutsPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsPCA",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

   if (Use["CutsGA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsGA",
                          "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

   if (Use["CutsSA"])
      factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsSA",
                          "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "Likelihood",
                          "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                          "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                          "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA");

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
                          "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50");

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
                          "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50");

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERS",
                          "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");

   if (Use["PDERSD"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERSD",
                          "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");

   if (Use["PDERSPCA"])
      factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERSPCA",
                          "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   // if (Use["PDEFoam"])
   //    factory->BookMethod(dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
   //   "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

   // if (Use["PDEFoamBoost"])
   //    factory->BookMethod(dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
   //                        "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod(dataloader, TMVA::Types::kKNN, "KNN",
                          "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod(dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None");

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod(dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss");

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod(dataloader, TMVA::Types::kFisher, "BoostedFisher",
                          "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MC",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_GA",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1");

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_SA",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

   if (Use["FDA_MT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");

   if (Use["FDA_GAMT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_GAMT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");

   if (Use["FDA_MCMT"])
      factory->BookMethod(dataloader, TMVA::Types::kFDA, "FDA_MCMT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");

   if (Use["MLPBFGS"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");

   if (Use["MLPBNN"])
      factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator"); // BFGS training with bayesian regulators

   // Multi-architecture DNN implementation.
   if (Use["DNN_CPU"] or Use["DNN_GPU"])
   {
      // General layout.
      TString layoutString("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

      // Define Training strategy. One could define multiple strategy string separated by the "|" delimiter

      TString trainingStrategyString = ("TrainingStrategy=LearningRate=1e-2,Momentum=0.9,"
                                        "ConvergenceSteps=20,BatchSize=100,TestRepetitions=1,"
                                        "WeightDecay=1e-4,Regularization=None,"
                                        "DropConfig=0.0+0.5+0.5+0.5");

      // General Options.
      TString dnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                         "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append(":");
      dnnOptions.Append(layoutString);
      dnnOptions.Append(":");
      dnnOptions.Append(trainingStrategyString);

      // Cuda implementation.
      if (Use["DNN_GPU"])
      {
         TString gpuOptions = dnnOptions + ":Architecture=GPU";
         factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_GPU", gpuOptions);
      }
      // Multi-core CPU implementation.
      if (Use["DNN_CPU"])
      {
         TString cpuOptions = dnnOptions + ":Architecture=CPU";
         factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
      }
   }
   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod(dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"); // n_cycles:#nodes:#nodes:...

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod(dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod(dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                          "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

   if (Use["BDT"]) // Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                          "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

   if (Use["BDTB"]) // Bagging
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
                          "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
                          "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

   if (Use["BDTF"]) // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTF",
                          "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod(dataloader, TMVA::Types::kRuleFit, "RuleFit",
                          "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");

   // For an example of the category classifier usage, see: TMVAClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // exit(0);

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch())
      TMVA::TMVAGui(outfileName);

   return 0;
}

int main(int argc, char **argv)
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i = 1; i < argc; i++)
   {
      TString regMethod(argv[i]);
      if (regMethod == "-b" || regMethod == "--batch")
         continue;
      if (!methodList.IsNull())
         methodList += TString(",");
      methodList += regMethod;
   }
   return TMVAClassificationNew(methodList);
}
