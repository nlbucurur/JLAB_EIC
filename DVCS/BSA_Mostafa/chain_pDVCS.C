//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_pDVCS(std::string name, std::string period, std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");
    if (period == "fall2019" || period == "all")
    {
        chain->AddFile(std::string(path + "0pDVCS_fall2019_FTPhotonsCorrected.root").c_str());
    }
    if (period == "spring2019" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0pDVCS_spring2019_FTPhotonsCorrected.root").c_str());
    }
    if (period == "spring2020" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0pDVCS_spring2020_FTPhotonsCorrected.root").c_str());
    }
    if (period == "fall2018" || period == "rgainbending" || period == "RGA")
    {
        //chain->AddFile(std::string("/mnt/d/Jlab_analysis/recon_files/0pDVCS_rgaInbending.root").c_str());
        chain->AddFile(std::string("/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/0pDVCS_fall2018_rga.root").c_str());
    }
    if (period == "rgaoutbending" || period == "RGA")
    {
        chain->AddFile(std::string("/mnt/d/Jlab_analysis/recon_files/").c_str());
    }
    if (period == "sebastian")
    {
        chain->AddFile(std::string("/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/stripped_data_Test_8_pDVCS.root").c_str());
    }
    if (period == "rad")
    {
        chain->AddFile(std::string(path + "0rad_rad.root").c_str());
    }
    if (period == "norad")
    {
        chain->AddFile(std::string(path + "0norad_rad.root").c_str());
    }
    if (period == "genrad")
    {
        chain->AddFile(std::string(path + "0genrad_rad.root").c_str());
    }

    return chain;
}
