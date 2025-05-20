//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_DATA_enpi0(std::string name, std::string period, std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");
    if (period == "fall2019" || period == "all")
    {
        chain->AddFile(std::string(path + "0enpi0_fall2019_FTPhotonsCorrected_MVA.root").c_str());
    }
    if (period == "spring2019" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0enpi0_spring2019_FTPhotonsCorrected_MVA.root").c_str());
    }
    if (period == "spring2020" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0enpi0_spring2020_FTPhotonsCorrected_MVA.root").c_str());
    }

    return chain;
}
