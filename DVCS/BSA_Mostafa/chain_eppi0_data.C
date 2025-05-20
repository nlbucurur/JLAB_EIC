//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_DATA_eppi0(std::string name, std::string period, std::string path = "/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");
    if (period == "fall2019" || period == "all")
    {
        chain->AddFile(std::string(path + "0eppi0_fall2019_FTPhotonsCorrected.root").c_str());
    }
    if (period == "spring2019" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0eppi0_spring2019_FTPhotonsCorrected.root").c_str());
    }
    if (period == "spring2020" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0eppi0_spring2020_FTPhotonsCorrected.root").c_str());
    }

    if (period == "fall2018" || period == "rgainbending" || period == "RGA")
    {
        chain->AddFile(std::string("/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/0eppi0_fall2018_rga.root").c_str());
    }

    // just to test flat simulations please ignore for later studies
    if (period == "flat")
    {
        chain->AddFile(std::string(path + "flat10p2_eppi0_simu_new_genpi.root").c_str());
    }
    return chain;
}
