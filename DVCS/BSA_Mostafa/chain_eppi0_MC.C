//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_MC_eppi0(std::string name, std::string energy = "", std::string period = "", std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    if (energy == "10p2" || energy == "")
    {
        // chain->AddFile(std::string(path + "1eppi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "eppi0_10p2_octobreVersion.root").c_str());
        //chain->AddFile(std::string(path + "10p2_eppi0_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p4_march2023.root").c_str());
        chain->AddFile(std::string(path + "0eppi0_10p2_bkgMerging.root").c_str());
    }

    if (energy == "flat10p2")
    {
        chain->AddFile(std::string(path + "flat10p2_eppi0_simu_new_genpi.root").c_str());
    }
    if (energy == "genflat10p2")
    {
        chain->AddFile(std::string(path + "genflat10p2_eppi0_simu_new_genpi.root").c_str());
    }
    if (energy == "gen10p2")
    {
        chain->AddFile(std::string(path + "gen10p2_eppi0_simu_new_genpi.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        // chain->AddFile(std::string(path + "2eppi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "eppi0_10p6_octobreVersion.root").c_str());
        //chain->AddFile(std::string(path + "10p6_eppi0_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p6_february2023.root").c_str());
        chain->AddFile(std::string(path + "0eppi0_10p6_bkgMerging.root").c_str());
    }
    if (energy == "10p4" || energy == "")
    {
        // chain->AddFile(std::string(path + "1eppi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "eppi0_10p4_octobreVersion.root").c_str());
        //chain->AddFile(std::string(path + "10p4_eppi0_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p4_march2023.root").c_str());
        chain->AddFile(std::string(path + "0eppi0_10p4_bkgMerging.root").c_str());
    }

    if (period == "fall2018" && (energy == "10p6" || energy == ""))
    {
        chain->AddFile(std::string(path + "0eppi0_simulation_rga.root").c_str());
    }

    return chain;
}
