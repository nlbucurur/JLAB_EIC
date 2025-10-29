#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <string>
#include <iostream>
#include <vector>

void rename_trees()
{
    std::string input_file_path = "/home/lorena/Documents/Thesis/Data_Analysis_Class/Analysis/";
    std::string output_file_path = "/home/lorena/Documents/Thesis/JLAB_EIC/DVCS/data/";
    std::vector<std::pair<std::string, std::string>> filenames = {
        {"Data_NP_Theta_g_5.root", "pDVCS_stripped"}};

    for (const auto &[filename, newTreeName] : filenames)
    {
        std::string input_path = input_file_path + filename;
        std::string output_path = output_file_path + filename;

        TFile *input_file = TFile::Open(input_path.c_str(), "READ");
        if (!input_file || input_file->IsZombie()) {
            std::cerr << "Error: Could not open input file: " << input_path << std::endl;
            continue;
        }

        TFile *output_file = TFile::Open(output_path.c_str(), "RECREATE");
        if (!output_file || output_file->IsZombie()) {
            std::cerr << "Error: Could not create output file: " << output_path << std::endl;
            input_file->Close();
            delete input_file;
            continue;
        }

        TIter next(input_file->GetListOfKeys());
        TKey *key;

        while ((key = (TKey *)next()))
        {
            std::string key_name = key->GetName();
            std::string class_name = key->GetClassName();

            if (class_name == "TTree" && key_name.find("pDVCS") == 0)
            {
                TTree *old_tree = dynamic_cast<TTree *>(input_file->Get(key_name.c_str()));
                if (old_tree)
                {
                    TTree *new_tree = old_tree->CloneTree(-1, "fast");  // Use the "fast" option for safer cloning
                    if (!new_tree) {
                        std::cerr << "Error: Could not clone tree '" << key_name << "' in file: " << input_path << std::endl;
                        continue;
                    }
                    new_tree->SetName(newTreeName.c_str());
                    new_tree->Write("", TObject::kOverwrite);
                    delete new_tree;
                    std::cout << "Successfully renamed tree '" << key_name
                              << "' to '" << newTreeName
                              << "' in file: " << filename << std::endl;
                }
            }
            else
            {
                TObject *obj = key->ReadObj();
                obj->Write();
                delete obj;
            }
        }

        output_file->Close();
        input_file->Close();
        delete output_file;
        delete input_file;
    }
}
