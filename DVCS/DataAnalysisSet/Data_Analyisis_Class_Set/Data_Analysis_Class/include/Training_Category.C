// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers in categorisation mode.
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAClassificationCategory
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables with category (eta) dependent
/// properties.
///
/// For this example, only Fisher and Likelihood are used. Run via:
///
///     root -l TMVAClassificationCategory.C
///
/// The output file "TMVACC.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
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
 
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
 
 
// two types of category methods are implemented
Bool_t UseOffsetMethod = kTRUE;
 
int BDT::TrainingCategory(TCut cutSB, TString MC_DVCS, TString MC_Pi0, vector<TString> vars)
{
  TString myMethodList = "BDT";
  bool useRandomSplitting = true;

   //---------------------------------------------------------------
   // Example for usage of different event categories with classifiers
 
   std::cout << std::endl << "==> Start TMVAClassificationCategory" << std::endl;
 
   // This loads the library
   TMVA::Tools::Instance();
 
   bool batchMode = true;
 
    // Load the signal and background event samples from ROOT trees
  TFile *input1(0);
  TFile *input2(0);

  TString fname_pi = MC_Pi0;
  TString fname_DVCS = MC_DVCS;

  input1 = TFile::Open( fname_pi );
  input2 = TFile::Open( fname_DVCS );

  TTree *signalTree     = (TTree*)input2->Get("pDVCS");
  TTree *background;
  if( input1->GetListOfKeys()->Contains("pDVCS") )
    background     = (TTree*)input1->Get("pDVCS");
  else
    background     = (TTree*)input1->Get("eppi0");
 
   // Global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;


  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  gSystem->ChangeDirectory(Folder);
  TString outfileName( "TMVACC.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
   // Create the factory object (see TMVAClassification.C for more information)
 
  std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D" );
  if (batchMode) factoryOptions += ":!Color:!DrawProgressBar";
 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassificationCategory", outputFile, factoryOptions );
 
   // Create DataLoader
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
 
   // Define the input variables used for the MVA training
  for(int l=0; l<vars.size(); l++)
    {
      dataloader->AddVariable(vars.at(l),'F');
    }
 
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   dataloader->AddSpectator( "strip_Ph_Theta" );
   dataloader->AddSpectator( "strip_Nuc_Theta" );
 
   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight     );
   dataloader->AddBackgroundTree( background, backgroundWeight );
 
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = cutSB; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = cutSB; // for example: TCut mycutb = "abs(var1)<0.5";
/*
   if(!categories)
   { 
    mycuts = cutSB + TCut("strip_Ph_Theta > 5");
    mycutb = cutSB + TCut("strip_Ph_Theta > 5");
   }
 */
   // Tell the factory how to use the training and testing events
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
 
   // Book MVA methods
   TString optBDT="!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:DoBoostMonitor";
   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",optBDT);  

 
   // Categorised classifier
   TMVA::MethodCategory* mcat = 0;
 
   // The variable sets
   TString theCat1Vars="";
   TString theCat2Vars;
   TString theCat3Vars;
   TString theCat4Vars;

   for(int k=0; k<vars.size();k++)
   {
   	theCat1Vars = theCat1Vars + vars.at(k);
   	if(k+1!=vars.size())
 		theCat1Vars = theCat1Vars + TString(":");
   }

   theCat2Vars = theCat1Vars;
   theCat3Vars = theCat1Vars;
   theCat4Vars = theCat1Vars;
  
  if(categories)
  { 
   // Categories
   std::cerr<<"\033[1;31m Creating trainings for photon topologies...\033[0m - "<<nft<<endl;

   TMVA::MethodBase* fiCat = factory->BookMethod( dataloader, TMVA::Types::kCategory, "Ph_Topology","" );
   mcat = dynamic_cast<TMVA::MethodCategory*>(fiCat);
   //mcat->AddMethod( "strip_Ph_Theta < 5 && strip_Nuc_Theta < 35", theCat1Vars, TMVA::Types::kBDT, "BDT_1",optBDT );
   //mcat->AddMethod( "strip_Ph_Theta > 5 && strip_Nuc_Theta < 35", theCat2Vars, TMVA::Types::kBDT, "BDT_2",optBDT );
   //mcat->AddMethod( "strip_Ph_Theta < 5 && strip_Nuc_Theta > 35", theCat3Vars, TMVA::Types::kBDT, "BDT_3",optBDT );
   //mcat->AddMethod( "strip_Ph_Theta > 5 && strip_Nuc_Theta > 35", theCat4Vars, TMVA::Types::kBDT, "BDT_4",optBDT );
   //mcat->AddMethod( "strip_Nuc_Theta < 35", theCat1Vars, TMVA::Types::kBDT, "BDT_1",optBDT );
   //mcat->AddMethod( "strip_Nuc_Theta > 35", theCat2Vars, TMVA::Types::kBDT, "BDT_2",optBDT );
   mcat->AddMethod( "strip_Ph_Theta < 5", theCat1Vars, TMVA::Types::kBDT, "BDT_1",optBDT );
   mcat->AddMethod( "strip_Ph_Theta > 5", theCat2Vars, TMVA::Types::kBDT, "BDT_2",optBDT );
  }
  else
  {
     std::cerr<<"\033[1;31m Creating FD training only. There is not enough background in the FT...\033[0m - "<<nft<<endl;
  }
 
   // Now you can tell the factory to train, test, and evaluate the MVAs
 
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
 
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
 
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
 
   // --------------------------------------------------------------
 
   // Save the output
   outputFile->Close();
 
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassificationCategory is done!" << std::endl;
 
   // Clean up
   delete background;
   input1->Close();
   input2->Close();
   delete input1;
   delete input2;
   delete factory;
   delete dataloader;
 
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
   gSystem->ChangeDirectory(dir);
   gSystem->Exec("pwd");
   
   return 0;
}
