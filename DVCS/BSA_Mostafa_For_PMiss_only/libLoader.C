void declare(std::string lib)
{
    gSystem->Load(lib.c_str());
    gSystem->AddLinkedLibs(lib.c_str());
};

void libCompile()
{
    // compile
    gROOT->ProcessLine(".L TWeighter1D.C+");
    gROOT->ProcessLine(".L RooMyPdf.cxx+");
    gROOT->ProcessLine(".L RooMyPdftemp.cxx+");
    gROOT->ProcessLine(".L RooApollonios.cxx+");
}

void libLoad()
{

    // Load necessary libraries
    gSystem->Load("libRooFit");
    gSystem->Load("libRooStats");
    gSystem->Load("libMinuit");
    // load
    declare("TWeighter1D_C.so");
    declare("RooMyPdf_cxx.so");
    declare("RooMyPdftemp_cxx.so");
    declare("RooApollonios_cxx.so");
    // -----
    using namespace RooFit;
    gSystem->Load("libMinuit");
}

void libLoader()
{
    libCompile();
    libLoad();

    gROOT->ProcessLine(".L chain_pDVCS.C");
    gROOT->ProcessLine(".L chain_pDVCS_MC.C");
    gROOT->ProcessLine(".L chain_eppi0_data.C");
    gROOT->ProcessLine(".L chain_eppi0_MC.C");
    gROOT->ProcessLine(".L chain_eppi01gamma_MC.C");

    return;
}
