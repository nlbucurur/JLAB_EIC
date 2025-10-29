#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TCut.h>
#include <iostream>

using namespace std;

void CreateWeightFile()
{
    // Load input files
    TFile *fData = TFile::Open("/home/lorena/Documents/Thesis/Data_Analysis_Class/Analysis/Data_NP_Theta_g_5.root");
    TFile *fMC = TFile::Open("/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/1pDVCS_simulation.root");

    // Load trees
    TTree *tData = (TTree*)fData->Get("pDVCS");
    TTree *tMC = (TTree*)fMC->Get("pDVCS_stripped");

    // Histograms
    int nBins = 100;
    float min = -1.0; // Adjust based on _mm2_eNg range
    float max = 5.0;

    TH1F *hData = new TH1F("hData", "", nBins, min, max);
    TH1F *hMC = new TH1F("hMC", "", nBins, min, max);
    TH1F *hWeight = new TH1F("c_mm2_eNg_cuts_weights", "", nBins, min, max);

    // Fill histograms
    tData->Draw("_mm2_eNg >> hData", "", "goff");
    tMC->Draw("_mm2_eNg >> hMC", "", "goff");

    // Calculate weights: Data / MC
    for (int i = 1; i <= nBins; i++)
    {
        double dataCount = hData->GetBinContent(i);
        double mcCount = hMC->GetBinContent(i);

        double weight = (mcCount > 0) ? dataCount / mcCount : 1.0;
        hWeight->SetBinContent(i, weight);
    }

    // Save weight histogram
    TFile *fOut = new TFile("weight_mm2_eNg_1D.root", "RECREATE");
    hWeight->Write();
    fOut->Close();

    cout << "Weight file created: weight_mm2_eNg_1D.root" << endl;
}
