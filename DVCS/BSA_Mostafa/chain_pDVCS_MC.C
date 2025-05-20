//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_MC_pDVCS(std::string name, std::string energy = "", std::string period = "", std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    // chain->AddFile(std::string(path + "eppi0_Paulpi0simulation.root").c_str());
    if (energy == "10p2" || energy == "")
    {
        chain->AddFile(std::string(path + "1pDVCS_simulation.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        chain->AddFile(std::string(path + "2pDVCS_simulation.root").c_str());
    }
    if (energy == "10p4" || energy == "")
    {
        chain->AddFile(std::string(path + "1pDVCS_simulation.root").c_str());
        chain->AddFile(std::string(path + "2pDVCS_simulation.root").c_str());
    }

    if (period == "fall2018" && (energy == "10p6" || energy == ""))
    {
        path = "/mnt/d/Jlab_analysis/new_recon_files/";
        chain->AddFile(std::string(path + "0pDVCS_simulation_rga.root").c_str());
    }

    return chain;
}
