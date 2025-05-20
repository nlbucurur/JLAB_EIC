//------------
// usage :
// gROOT->ProcessLine(".L pathtochain/chainname.C");
// TChain* chain = linked(tupleName);
//------------

TChain *linked_nDVCS(std::string name, std::string period, std::string path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/")
{
    cout << " =============================================================================" << endl;
    cout << " === path  = " << path << endl;
    cout << " === tuple = " << name << endl;
    cout << " =============================================================================" << endl;
    TChain *chain = new TChain(name.c_str(), "");

    if (period == "fall2019" || period == "all")
    {
        chain->AddFile(std::string(path + "0nDVCS_fall2019_FTPhotonsCorrected_MVA.root").c_str());
    }

    if (period == "spring2019" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0nDVCS_spring2019_FTPhotonsCorrected_MVA.root").c_str());
    }

    if (period == "spring2020" || period == "all" || period == "inbending")
    {
        chain->AddFile(std::string(path + "0nDVCS_spring2020_FTPhotonsCorrected_MVA.root").c_str());
    }

    if (period == "rgainbending" || period == "RGA")
    {
        chain->AddFile(std::string("/projet/PRAE/hoballah/theData/DataJlab_analysis/NeutronID/first_attempt_CLAS12PID_required/0nDVCS_NID_studies_rgaInbending.root").c_str());
    }
    if (period == "rgaoutbending" || period == "RGA")
    {
        chain->AddFile(std::string("/projet/PRAE/hoballah/theData/Data/Jlab_analysis/NeutronID/first_attempt_CLAS12PID_required/").c_str());
    }

    if (period == "cut_test_spring2019")
    {
        chain->AddFile(std::string(path + "0nDVCS_spring2019_FTPhotonsCorrected_test_cuts_for_note_MVA.root").c_str());
    }

    return chain;
}
