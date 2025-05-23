
TH1* BDT::Get_Contamination_Mostafa_Maxime(TCut cut,  vector<double> Pbins)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TChain *Pi01g= new TChain("pDVCS");
  //Pi01g->Add("/home/munoz/Datasets/DVMP/Quality_Pi_as_DVCS_P.root");
  Pi01g->Add(Folder + TString("Tested_1gamma.root"));
  Pi01g->Add(Folder + TString("Tested_1gamma_1.root"));

  TChain *Pi02g= new TChain("eppi0");
  Pi02g->Add("/home/munoz/Datasets/DVMP/P/Quality_Sim_eppi0_P.root");
  Pi02g->Add("/home/munoz/Datasets/DVMP/P/Quality_Sim_eppi0_P_1.root");

  //Use quality 3 as it has the cuts for RGA eppi0 selection
  TChain *Pi02gData= new TChain("eppi0");
  Pi02gData->Add("~/Quality3_Data_eppi0_P.root");

  TChain *Data= new TChain("pDVCS");
  Data->Add(Folder + TData);

  
  //Convert cut to string
  TString String_cut = cut.GetTitle();
  //replace Ph by Pi0
  String_cut.ReplaceAll("Ph ", "Pi0 ");
  String_cut.ReplaceAll("Ph>", "Pi0>");
  String_cut.ReplaceAll("Ph<", "Pi0<");
  String_cut.ReplaceAll("strip_Ph_P", "strip_W");
  String_cut.ReplaceAll("gamma", "Pi0");
  //Create cut for 2gamma case
  TCut cut2g = TCut(String_cut) + TCut("(strip_Ph1_P>2 || strip_Ph2_P>2)");

  int N_Phi=Pbins.size()-1;

  TH1F *Pi01g_p = new TH1F("Pi01g_p","",N_Phi,0,360);
  TH1F *Pi01g_m = new TH1F("Pi01g_m","",N_Phi,0,360);
  TH1F *Pi01g_p_BDT = new TH1F("Pi01g_p_BDT","",N_Phi,0,360);
  TH1F *Pi01g_m_BDT = new TH1F("Pi01g_m_BDT","",N_Phi,0,360);
  
  TH1F *Pi02g_p = new TH1F("Pi02g_p","",N_Phi,0,360);
  TH1F *Pi02g_m = new TH1F("Pi02g_m","",N_Phi,0,360);
  
  TH1F *Pi02g_Data_p = new TH1F("Pi02g_Data_p","",N_Phi,0,360);
  TH1F *Pi02g_Data_m = new TH1F("Pi02g_Data_m","",N_Phi,0,360);

  TH1F *Data_p = new TH1F("Data_p","",N_Phi,0,360);
  TH1F *Data_m = new TH1F("Data_m","",N_Phi,0,360);
  TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",N_Phi,0,360);
  TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",N_Phi,0,360);

  Pi01g_p->Sumw2();
  Pi01g_m->Sumw2();
  Pi01g_p_BDT->Sumw2();
  Pi01g_m_BDT->Sumw2();
  Pi02g_p->Sumw2();
  Pi02g_m->Sumw2();
  Pi02g_Data_p->Sumw2();
  Pi02g_Data_m->Sumw2();
  Data_p->Sumw2();
  Data_m->Sumw2();

  Pi01g->Project("Pi01g_p", "Phi_Ph", cut);
  Pi01g->Project("Pi01g_m", "Phi_Ph", cut);
  Pi01g->Project("Pi01g_p_BDT", "Phi_Ph", cut + TCut("_strip_Nuc_BDT > 0.0"));
  Pi01g->Project("Pi01g_m_BDT", "Phi_Ph", cut + TCut("_strip_Nuc_BDT > 0.0"));

  Pi02g->Project("Pi02g_p", "Phi_Pi0", cut2g + TCut("strip_Pi0_2DChi2 < 0.04 && abs(mm2_eNgg)<0.01"));
  Pi02g->Project("Pi02g_m", "Phi_Pi0", cut2g + TCut("strip_Pi0_2DChi2 < 0.04 && abs(mm2_eNgg)<0.01"));
  
  Pi02gData->Project("Pi02g_Data_p", "Phi_Pi0", cut2g + TCut("Helicity>0") + TCut("strip_Pi0_2DChi2 < 0.04 && abs(mm2_eNgg)<0.01"));
  Pi02gData->Project("Pi02g_Data_m", "Phi_Pi0", cut2g + TCut("Helicity<0") + TCut("strip_Pi0_2DChi2 < 0.04 && abs(mm2_eNgg)<0.01"));

  Data->Project("Data_p", "Phi_Ph", cut + TCut("Helicity>0"));
  Data->Project("Data_m", "Phi_Ph", cut + TCut("Helicity<0"));
  Data->Project("Data_p_BDT", "Phi_Ph", cut + TCut("Helicity>0 && _strip_Nuc_BDT > 0.0"));
  Data->Project("Data_m_BDT", "Phi_Ph", cut + TCut("Helicity<0 && _strip_Nuc_BDT > 0.0"));

  Pi01g_p->Divide(Pi02g_p);
  for(int i=1; i<=N_Phi; i++)
    {
      std::cout<<Pi01g_p->GetBinContent(i)<<" "<<Pi01g_p->GetBinError(i)<<endl;
    }
  std::cout<<endl;
  Pi01g_p->Multiply(Pi02g_Data_p);

  Pi01g_p_BDT->Divide(Pi02g_p);
  for(int i=1; i<=N_Phi; i++)
    {
      std::cout<<Pi01g_p_BDT->GetBinContent(i)<<" "<<Pi01g_p_BDT->GetBinError(i)<<endl;
    }
  std::cout<<endl;
  Pi01g_p_BDT->Multiply(Pi02g_Data_p);


  Pi01g_m->Divide(Pi02g_m);
  Pi01g_m->Multiply(Pi02g_Data_m);

  Pi01g_m_BDT->Divide(Pi02g_m);
  Pi01g_m_BDT->Multiply(Pi02g_Data_m);


  double est1=(Pi01g_p->Integral() + Pi01g_m->Integral())*1.0/(Data_p->Integral() + Data_m->Integral());
  double est2=(Pi01g_p_BDT->Integral() + Pi01g_m_BDT->Integral())*1.0/(Data_p_BDT->Integral() + Data_m_BDT->Integral());
  std::cout<<"Mostafa estimaation "<<est1<<" "<<est2<<endl;
  boundaries.push_back(est1);
  boundaries.push_back(est2);
  
  //TH1F BSA = -1.0*((*Data_p-*Pi01g_p) - (*Data_m-*Pi01g_m))/((*Data_p-*Pi01g_p) + (*Data_m-*Pi01g_m));
  Pi01g_p_BDT->Scale(-1.0);
  Pi01g_m_BDT->Scale(-1.0);

  Data_p_BDT->Add(Pi01g_p_BDT);
  Data_m_BDT->Add(Pi01g_m_BDT);
  
  TH1 *BA= Data_m_BDT->GetAsymmetry(Data_p_BDT);
  TCanvas* c2 = new TCanvas("c2","Histograms");

  
  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf = new TF1("fitf","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf->SetParameter(0,0.1);
  fitf->SetParameter(1,-0.3);
  fitf->SetParLimits(0,0,0.3);
  fitf->SetParLimits(1,-1.,1.);
  fitf->SetLineColor(kBlue);
	  
  BA->Fit("fitf","Q");
  BA->Scale(1.0/Bpol);
  BA->SetTitle("Hall B method");
  BA->SetLineColor(kBlue);
  BA->SetMarkerColor(kBlue);
  BA->SetAxisRange(-0.3, 0.3,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi (deg)");
  BA->Draw();

  c2->Print(Folder + TString("Mostafa.pdf"));

  delete Pi01g;
  delete Pi02g;
  delete Pi02gData;
  delete Data;

  delete Pi01g_p;
  delete Pi01g_m;
  delete Pi01g_p_BDT;
  delete Pi01g_m_BDT;
  delete Pi02g_p;
  delete Pi02g_m;
  delete Pi02g_Data_p;
  delete Pi02g_Data_m;
  delete Data_p;
  delete Data_m;
  delete Data_p_BDT;
  delete Data_m_BDT;
  delete c2;

  return BA;

}
