//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_MC_nDVCS(std::string name, std::string energy = "", std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    if (energy == "10p2" || energy == "")
    {
        chain->AddFile(std::string(path + "1nDVCS_simulation_ProtonContamination_studies_MVA_response.root").c_str());
    }
    if (energy == "10p6" || energy == "")
    {
        chain->AddFile(std::string(path + "2nDVCS_simulation_ProtonContamination_studies_MVA_response.root").c_str());
    }
    if (energy == "10p4")
    {
        chain->AddFile(std::string(path + "1nDVCS_simulation_ProtonContamination_studies_MVA_response.root").c_str());
        chain->AddFile(std::string(path + "2nDVCS_simulation_ProtonContamination_studies_MVA_response.root").c_str());
    }
    return chain;
}
