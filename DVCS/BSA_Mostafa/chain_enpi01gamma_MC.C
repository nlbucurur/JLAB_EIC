//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_MC_enpi01g(std::string name, std::string energy = "", std::string path = "/projet/PRAE/hoballah/theData/Data/Jlab_analysis/new_recon_files/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    if (energy == "10p2" || energy == "")
    {
        // chain->AddFile(std::string(path + "4nDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "nDVCS_10p2_PartialPi0_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        // chain->AddFile(std::string(path + "10p2nDVCS_dvmpasdvcs_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p2_dvmpasdvcs_march2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0nDVCS_Pi0dataAsDVCS_10p2_MVA.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        // chain->AddFile(std::string(path + "5nDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "nDVCS_10p6_PartialPi0_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        // chain->AddFile(std::string(path + "10p6nDVCS_dvmpasdvcs_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p6_dvmpasdvcs_february2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0nDVCS_Pi0dataAsDVCS_10p6_MVA.root").c_str());
    }
    if (energy == "10p4" || energy == "")
    {
        // chain->AddFile(std::string(path + "nDVCS_10p4_PartialPi0_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        // chain->AddFile(std::string(path + "4nDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "5nDVCS_simulation.root").c_str());
        // chain->AddFile(std::string(path + "10p4nDVCS_dvmpasdvcs_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p4_dvmpasdvcs_march2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0nDVCS_Pi0dataAsDVCS_10p4_MVA.root").c_str());
    }
    return chain;
}
