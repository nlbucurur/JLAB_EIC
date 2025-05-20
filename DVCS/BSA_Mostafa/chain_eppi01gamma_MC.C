//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

// simulations previous to new genpi version are in recon_files and not in new_recon_files

TChain *linked_MC_eppi01g(std::string name, std::string energy = "", std::string period = "", std::string path = "/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    if (energy == "10p2" || energy == "")
    {
        // chain->AddFile(std::string(path + "4pDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "pDVCS_10p2_PartialPi0_octobreVersion.root").c_str());
        // chain->AddFile(std::string(path + "10p2pDVCS_dvmpasdvcs_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p2_dvmpasdvcs_march2023.root").c_str());
        chain->AddFile(std::string(path + "0pDVCS_Pi0dataAsDVCS_10p2.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        // chain->AddFile(std::string(path + "5pDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "pDVCS_10p6_PartialPi0_octobreVersion.root").c_str());
        // chain->AddFile(std::string(path + "10p6pDVCS_dvmpasdvcs_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p6_dvmpasdvcs_february2023.root").c_str());
        chain->AddFile(std::string(path + "0pDVCS_Pi0dataAsDVCS_10p6.root").c_str());
    }
    if (energy == "10p4")
    {
        // chain->AddFile(std::string(path + "4pDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "5pDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "10p4pDVCS_dvmpasdvcs_simu_new_genpi.root").c_str());
        //chain->AddFile(std::string(path + "0eppi0_10p4_dvmpasdvcs_march2023.root").c_str());
        chain->AddFile(std::string(path + "0pDVCS_Pi0dataAsDVCS_10p4.root").c_str());
    }

    if (period == "fall2018" && (energy == "10p6" || energy == ""))
    {
        chain->AddFile(std::string(path + "0eppi01gamma_simulation_rga.root").c_str());
    }

    return chain;
}
