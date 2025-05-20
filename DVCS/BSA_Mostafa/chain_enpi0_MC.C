//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_MC_enpi0(std::string name, std::string energy = "", std::string path = ""/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/"")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    // chain->AddFile(std::string(path + "enpi0_Paulpi0simulation.root").c_str());
    if (energy == "10p2" || energy == "")
    {
        // chain->AddFile(std::string(path + "1enpi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "enpi0_10p2_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        //chain->AddFile(std::string(path + "10p2_enpi0_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p2_march2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0enpi0_10p2_bkgMerging_MVA.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        // chain->AddFile(std::string(path + "2enpi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "enpi0_10p6_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        //chain->AddFile(std::string(path + "10p6_enpi0_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p6_february2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0enpi0_10p6_bkgMerging_MVA.root").c_str());
    }
    if (energy == "10p4" || energy == "")
    {
        // chain->AddFile(std::string(path + "enpi0_10p4_octobreVersion_ProtonContamination_MVA_response.root").c_str());
        // chain->AddFile(std::string(path + "1enpi0_simu.root").c_str());
        // chain->AddFile(std::string(path + "2enpi0_simu.root").c_str());
        //chain->AddFile(std::string(path + "10p4_enpi0_simu_new_genpi_MVA.root").c_str());
        //chain->AddFile(std::string(path + "0enpi0_10p4_march2023_MVA.root").c_str());
        chain->AddFile(std::string(path + "0enpi0_10p4_bkgMerging_MVA.root").c_str());
    }
    return chain;
}
