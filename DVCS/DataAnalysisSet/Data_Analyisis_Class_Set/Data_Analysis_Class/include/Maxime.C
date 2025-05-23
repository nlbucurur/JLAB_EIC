double BDT::GetPhiDistance(TLorentzVector& gamma){
  double phi_dist=atan2(gamma.Py(),gamma.Px())*TMath::RadToDeg();
  if (phi_dist<0) phi_dist+=360;
  if (fabs(phi_dist-TMath::Floor((phi_dist+30)/60.)*60.)==0) return 0.001;
  return fabs(phi_dist-TMath::Floor((phi_dist+30)/60.)*60.);
}

//Method to determine whether it is an exclusive photon
double BDT::CheckPi0Thr(TLorentzVector photon){
  if (photon.E()<2||photon.E()>9.0) return -1;
  double phi_temp=atan2(photon.Py(),photon.Px());
  double photon_angle=acos(photon.Pz()/photon.E())*TMath::RadToDeg();
  if (photon_angle<4.25){
    if (phi_temp<0) phi_temp=phi_temp+2*TMath::Pi();
    return MapSector[0][(int)floor((photon.E()-2)/(9.-2.)*70+0.5)]->Interpolate(photon.Pz()/photon.E(),phi_temp);
  }
  else{
    if (phi_temp<-30*TMath::DegToRad()) phi_temp=phi_temp+2*TMath::Pi();
    int num_sector=(phi_temp*TMath::RadToDeg()+30)/60+1;
    return MapSector[num_sector][(int)floor((photon.E()-2)/(9.-2.)*70+0.5)]->Interpolate(photon.Pz()/photon.E(),phi_temp);
  }
  return -1;
}


void BDT::GeneratePhotonMomentum(TLorentzVector& pionMomentum, TLorentzVector& photon1, TLorentzVector& photon2) {
  // define the particle masses and energies
  double masses[2] = { 0., 0.} ;
 
  // generate the decay
  TGenPhaseSpace decay;
  decay.SetDecay(pionMomentum, 2, masses);
  Double_t weight = decay.Generate();
  photon1 = *decay.GetDecay(0);
  photon2 = *decay.GetDecay(1);
}

int BDT::Is_DVCS(TTree*& pDVCS_tree, TVector3*& vertex, TLorentzVector*& electron, TLorentzVector*& photon1, TLorentzVector*& photon2, TLorentzVector*& Nuc)
{
  double thr_acc=0.05;

  bool fiducial1 = Fiducial_cut(*vertex,photon1->Vect());
  bool fiducial2 = Fiducial_cut(*vertex,photon2->Vect());

  bool ph1 = Check_DVCS_cuts(electron, photon1, Nuc);
  bool ph2 = Check_DVCS_cuts(electron, photon2, Nuc);

  bool acc1 = CheckPi0Thr(*photon1) > thr_acc;
  bool acc2 = CheckPi0Thr(*photon2) > thr_acc;
  /*
    if(photon1->Theta()*TMath::RadToDeg() <5 && photon1->Theta()*TMath::RadToDeg() > 2.5 && photon1->E() > 2 && ph1)
    std::cout<<"lala "<<fiducial1<<" "<<ph1<<" "<<acc1<<" "<<photon1->E()<<endl;

    if(photon2->Theta()*TMath::RadToDeg() <5 && photon2->Theta()*TMath::RadToDeg() > 2.5 && photon2->E() > 2 && ph2)
    std::cout<<"lala "<<fiducial2<<" "<<ph2<<" "<<acc2<<" "<<photon2->E()<<endl;

  */
  if(ph1 && fiducial1 && acc1 && ph2 && fiducial2 && acc2)
    {
      if(( *beam + *target - *electron - *photon1 - *Nuc).M2() < ( *beam + *target - *electron - *photon2 - *Nuc).M2())
	{
	  Build_DVCS_Tree(pDVCS_tree, electron, photon1, Nuc);
	  return 1;
	}
      else
	{
	  Build_DVCS_Tree(pDVCS_tree, electron, photon2, Nuc);
	  return 2;
	}
    }

  if(ph1 && fiducial1 && acc1 && !(ph2 && fiducial2 && acc2))
    {
      Build_DVCS_Tree(pDVCS_tree, electron, photon1, Nuc);
      return 1;
    }
  if(!(ph1 && fiducial1 && acc1 ) && ph2 && fiducial2 && acc2)
    {
      Build_DVCS_Tree(pDVCS_tree, electron, photon2, Nuc);
      return 2;
    }
  return 0;  
}

bool BDT::Is_DVMP(TVector3*& vertex, TLorentzVector*& photon1, TLorentzVector*& photon2)
{
  double thr_acc=0.05;
  double thr_pi0acc=0.01;

  bool angles1 = photon1->Theta() > 2*TMath::Pi()/180. && photon1->Theta() < 65*TMath::Pi()/180.;
  bool angles2 = photon2->Theta() > 2*TMath::Pi()/180. && photon2->Theta() < 65*TMath::Pi()/180.;
  
  bool energies = photon1->E() > 0.65 && photon2->E() > 0.65;

  bool fiducial1 = Fiducial_cut(*vertex,photon1->Vect());
  bool fiducial2 = Fiducial_cut(*vertex,photon2->Vect());

  bool acc = CheckPi0Thr(*photon1 + *photon2) > thr_pi0acc;
  
  return angles1 && angles2 && energies && fiducial1 && fiducial2 && acc;
}

bool BDT::Check_DVCS_cuts(TLorentzVector*& electron, TLorentzVector*& photon, TLorentzVector*& Nuc)
{
  double Pmass=0.938;
  double Ebeam=10.6;

  double number=gRandom->Rndm();
  if (photon->E()<2 || (photon->Theta() < 5*TMath::Pi()/180 && number<0.1) || (photon->Theta() > 5*TMath::Pi()/180 && number<0.045))
    {
      return false;
    }
  
  //Kinematics
  double Q2 = 4 * 10.6 * electron->P() * TMath::Power(TMath::Sin(electron->Theta() / 2), 2.);
  double W = TMath::Sqrt(0.938*0.938 + 2 * 0.938 * (10.6 - electron->P()) - Q2);
  double Xbj = Q2 / (2 * 0.938 * (10.6 - electron->P()));
  
      
  //delta_phi 
  TVector3 VelectronIn = beam->Vect();
  TVector3 VelectronOut = electron->Vect();
  TVector3 VnucleonOut = Nuc->Vect();
  TVector3 VphotonOut = photon->Vect();
  TVector3 Vvirtualphoton = (*beam - *electron).Vect();
  
  TVector3 Vlepto = VelectronIn.Cross(VelectronOut);
  TVector3 Vhadro = VnucleonOut.Cross(Vvirtualphoton);
  TVector3 VhadroPP = VnucleonOut.Cross(VphotonOut);
  
  double Phi_Nuc_temp = 180. / TMath::Pi() * Vlepto.Angle(Vhadro);
  double Phi_Ph_temp = 180. / TMath::Pi() * Vlepto.Angle(VhadroPP);
  
  if (Vlepto.Dot(VnucleonOut) > 0.)
    Phi_Nuc_temp = 360. - Phi_Nuc_temp;
  if (Vlepto.Dot(VphotonOut) < 0.)
    Phi_Ph_temp = 360. - Phi_Ph_temp;
  
  double delta_phi=Phi_Nuc_temp - Phi_Ph_temp;
  
  //delta_t
  double t_Nuc=(*Nuc - *target).M2();

  double ratio= Q2/(2.0*Pmass*Xbj);
  double cos=(-electron->Px()*photon->Px() - electron->Py()*photon->Py() + (Ebeam - electron->Pz())*photon->Pz())/(photon->P()*sqrt( pow(electron->P(),2) + Ebeam*Ebeam - 2*Ebeam*electron->P()*TMath::Cos(electron->Theta()) ) );
  double t_Ph = -(Q2*Pmass + (Q2/Xbj)*( ratio - cos*sqrt(Q2 + pow(ratio,2)) ) )/(Pmass + ratio - cos*sqrt(Q2 + pow(ratio,2)) );
  //std::cout<<t_Ph<<endl;
  double delta_t=t_Nuc - t_Ph;
  
  //Missing momentum
  TLorentzVector BalV = *beam + *target - *photon  - *electron - *Nuc;
  double miss_mom_eNg = TMath::Sqrt(pow(BalV.X(), 2) + pow(BalV.Y(), 2) + pow(BalV.Z(), 2));
  
  //theta_gamma_e
  double theta_gamma_e=180. / TMath::Pi() * TMath::ACos((VphotonOut.Dot(VelectronOut)) / (VphotonOut.Mag() * VelectronOut.Mag()));

  /*
    if(photon->Theta() < 5*TMath::Pi()/180 && photon->Theta() > 2.Z5*TMath::Pi()/180 && photon->E() > 2)
    {
    std::cout<<delta_t<<" "<<delta_phi<<" "<<" "<<miss_mom_eNg<<" "<<t_Ph<<" "<<Q2<<" "<<W<<" "<<Xbj<<" "<<theta_gamma_e<<endl;
    }
  */
  return abs(delta_phi)<2 && abs(delta_t)<2 && miss_mom_eNg<1 && t_Ph<0 && Q2>1 && W>2 && Xbj > 0 && Xbj < 1 && theta_gamma_e>5;

}

void BDT::fill_DVCS_histograms(TLorentzVector*& electron, TLorentzVector*& photon, TLorentzVector*& Nuc, TH1F*& mm2_eg, TH1F*& Phi_Ph)
{
  //mm2_eg
  double mm2 = (*beam + *target - *photon  - *electron).M2();

  //delta_phi 
  TVector3 VelectronIn = beam->Vect();
  TVector3 VelectronOut = electron->Vect();
  TVector3 VnucleonOut = Nuc->Vect();
  TVector3 VphotonOut = photon->Vect();
  
  TVector3 Vlepto = VelectronIn.Cross(VelectronOut);
  TVector3 VhadroPP = VnucleonOut.Cross(VphotonOut);
  
  double Phi = 180. / TMath::Pi() * Vlepto.Angle(VhadroPP);
  
  if (Vlepto.Dot(VphotonOut) < 0.)
    Phi = 360. - Phi;

  //fill
  mm2_eg->Fill(mm2);
  Phi_Ph->Fill(Phi);
}

TH1* BDT::Maxime(TCut cut, double BDT_cut, int bin_number, int Nphi)
{
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  TChain *ch1= new TChain("eppi0");
  ch1->Add(eppi0);

  ch1->SetBranchStatus("*",0);
  ch1->SetBranchStatus("Helicity",1);
  ch1->SetBranchStatus("EventNumber",1);
  ch1->SetBranchStatus("RunNumber",1);
  ch1->SetBranchStatus("bestCandidateFlag",1);

  ch1->SetBranchStatus("strip_Pi0_px",1);
  ch1->SetBranchStatus("strip_Pi0_py",1);
  ch1->SetBranchStatus("strip_Pi0_pz",1);
  ch1->SetBranchStatus("strip_Pi0_E",1);
  ch1->SetBranchStatus("strip_Pi0_IM2",1);
  
  ch1->SetBranchStatus("strip_El_px",1);
  ch1->SetBranchStatus("strip_El_py",1);
  ch1->SetBranchStatus("strip_El_pz",1);
  ch1->SetBranchStatus("strip_El_E",1);

  ch1->SetBranchStatus("strip_El_vx",1);
  ch1->SetBranchStatus("strip_El_vy",1);
  ch1->SetBranchStatus("strip_El_vz",1);
  
  ch1->SetBranchStatus("strip_Nuc_px",1);
  ch1->SetBranchStatus("strip_Nuc_py",1);
  ch1->SetBranchStatus("strip_Nuc_pz",1);
  ch1->SetBranchStatus("strip_Nuc_E",1);
  
  
  int helicity, EvNmbr, RnNmbr, N=1500;
  double N2gamma=0, N1gamma=0;
  TH1* BAerr;

  static vector<int> *flag;
  static vector<double> *Pi0_px;
  static vector<double> *Pi0_py;
  static vector<double> *Pi0_pz;
  static vector<double> *Pi0_E;
  static vector<double> *Pi0_Mass;
  
  static vector<double> *El_px;
  static vector<double> *El_py;
  static vector<double> *El_pz;
  static vector<double> *El_E;

  static vector<double> *El_vx;
  static vector<double> *El_vy;
  static vector<double> *El_vz;

  static vector<double> *Nuc_px;
  static vector<double> *Nuc_py;
  static vector<double> *Nuc_pz;
  static vector<double> *Nuc_E;
  
  static vector<double> *t_Pi0;
  static vector<double> *strip_Q2;
  static vector<double> *strip_Xbj;

  
  TLorentzVector *Pi0 = new TLorentzVector();
  TLorentzVector *Ph1 = new TLorentzVector();
  TLorentzVector *Ph2 = new TLorentzVector();
  TLorentzVector *El = new TLorentzVector();
  TLorentzVector *Nuc = new TLorentzVector();
  TVector3 *Vertex = new TVector3();

  ch1->SetBranchAddress("bestCandidateFlag", &flag);
  ch1->SetBranchAddress("Helicity", &helicity);
  ch1->SetBranchAddress("EventNumber", &EvNmbr);
  ch1->SetBranchAddress("RunNumber", &RnNmbr);
  
  ch1->SetBranchAddress("t_Pi0", &t_Pi0);
  ch1->SetBranchAddress("strip_Q2", &strip_Q2);
  ch1->SetBranchAddress("strip_Xbj", &strip_Xbj);

  ch1->SetBranchAddress("strip_Pi0_px", &Pi0_px);
  ch1->SetBranchAddress("strip_Pi0_py", &Pi0_py);
  ch1->SetBranchAddress("strip_Pi0_pz", &Pi0_pz);
  ch1->SetBranchAddress("strip_Pi0_E", &Pi0_E);
  ch1->SetBranchAddress("strip_Pi0_IM2", &Pi0_Mass);

  ch1->SetBranchAddress("strip_El_px", &El_px);
  ch1->SetBranchAddress("strip_El_py", &El_py);
  ch1->SetBranchAddress("strip_El_pz", &El_pz);
  ch1->SetBranchAddress("strip_El_E", &El_E);

  ch1->SetBranchAddress("strip_El_vx", &El_vx);
  ch1->SetBranchAddress("strip_El_vy", &El_vy);
  ch1->SetBranchAddress("strip_El_vz", &El_vz);

  ch1->SetBranchAddress("strip_Nuc_px", &Nuc_px);
  ch1->SetBranchAddress("strip_Nuc_py", &Nuc_py);
  ch1->SetBranchAddress("strip_Nuc_pz", &Nuc_pz);
  ch1->SetBranchAddress("strip_Nuc_E", &Nuc_E);

  //Loading 2g-inefficiency from GEMC
  TFile* ineff2g=new TFile(maps_path + TString("Neutron_diag_2D.root"));
  TH3F *corr2g=(TH3F*)gROOT->FindObject("eff");

  //Loading histos for pi0 acceptance
  TFile* MapE[7];
  for (int sec=0;sec<7;sec++){
    MapE[sec]=new TFile(maps_path + Form("MapSector_%d.root",sec));
    for (int eneb=0;eneb<71;eneb++) MapSector[sec][eneb]=(TH2F*)gROOT->FindObject(Form("Map_%d",eneb));
  }
    
  if(generate)
    {
      //TFile* out_sim=new TFile(Form("Maxime_contamination/Maxime_pi0_%i.root",bin_number), "RECREATE");
      TFile* out_sim=new TFile(Maxime_bkg + Form("Maxime_pi0_%i.root",bin_number), "RECREATE");
      TTree *tree = new TTree("pDVCS","pDVCS");
      tree->SetMaxTreeSize(100000000000LL);
      init_tree(tree);
      TBranch* Wb = tree->Branch("Weight", &weight);

      int *DVCS = new int(0);
      double ineff;
      int k0=0;
      double number=gRandom->Rndm();

      for(int i=0; i<ch1->GetEntries(); i++)
	{
	  printProgress(i*1.0/ch1->GetEntries());
	  //std::cout<<i<<" "<<ch1->GetEntries()<<endl;
	  ch1->GetEntry(i);
	  EventNumber=EvNmbr;
	  RunNumber=RnNmbr;
	  Helicity=helicity;
      
	  for(int j=0; j<flag->size();j++)
	    {
	      if(flag->at(j)==1 && t_Pi0->at(j) > boundaries.at(0) && t_Pi0->at(j) < boundaries.at(1) && strip_Q2->at(j) > boundaries.at(2) && strip_Q2->at(j) < boundaries.at(3) && strip_Xbj->at(j) > boundaries.at(4) && strip_Xbj->at(j) < boundaries.at(5))
		{
		  Pi0->SetXYZT(Pi0_px->at(j),Pi0_py->at(j),Pi0_pz->at(j),Pi0_E->at(j));
		  El->SetXYZT(El_px->at(j),El_py->at(j),El_pz->at(j),El_E->at(j));
		  Nuc->SetXYZT(Nuc_px->at(j),Nuc_py->at(j),Nuc_pz->at(j),Nuc_E->at(j));
		  Vertex->SetXYZ(El_vx->at(j),El_vy->at(j),El_vz->at(j));
		  strip_El_vz=El_vz->at(j);

		  N1gamma=0;
		  N2gamma=0;
		  //Generate 1500 decays of Pi0
		  for(int k=0; k<N; k++)
		    {
		      GeneratePhotonMomentum(*Pi0, *Ph1, *Ph2);
		      //Check if the decay photons would have been detected and the scenario
		      //Check if both photons are detected
		      if(Is_DVMP(Vertex, Ph1, Ph2))
			{
			  if (Ph1->Theta() > 5*TMath::Pi()/180 && Ph2->Theta() > 5*TMath::Pi()/180)//Both photons in FD.
			    {
			      if (Pi0->Theta()*TMath::RadToDeg()>5)
				ineff=corr2g->Interpolate(GetPhiDistance(*Pi0),Pi0->Theta()*TMath::RadToDeg(),Pi0->E());
			      else ineff=1;

			      if (ineff!=0) N2gamma+=ineff;
			      else N2gamma+=0.91;
			    } 
			  else if (Ph1->Theta() < 5*TMath::Pi()/180 && Ph2->Theta() < 5*TMath::Pi()/180) N2gamma+=0.81; //Both photon in FT
			  else if (Ph1->Theta() > 5*TMath::Pi()/180 || Ph2->Theta() > 5*TMath::Pi()/180) N2gamma+=0.85; //One Photon in FD.
			  else std::cout<<"Error in the Pion data"<<endl;
			}  
		      //Check if the event pass the DVCS analysis		  
		  
		      *DVCS=Is_DVCS(tree, Vertex, El, Ph1, Ph2, Nuc);
		      if(*DVCS==1 || *DVCS==2) N1gamma+=1;
		    }

		  if(N2gamma*1.0/N>0.01)
		    {
		      for(int k=k0; k<N1gamma;k++)
			{
			  tree->GetEntry(k);
			  weight = 1.0/N2gamma;
			  Wb->Fill();
			}
		      k0=N1gamma;
		    }
		  else
		    {
		      for(int k=k0; k<N1gamma;k++)
			{
			  tree->GetEntry(k);
			  weight = 0.0;
			  Wb->Fill();
			}
		      k0=N1gamma;
		    }	      
		}
	    }
	}

      tree->Write();
  
      delete flag;
      delete Pi0_px;
      delete Pi0_py;
      delete Pi0_pz;
      delete Pi0_E;
      delete Pi0_Mass;
      
      delete El_px;
      delete El_py;
      delete El_pz;
      delete El_E;
    
      delete El_vx;
      delete El_vy;
      delete El_vz;
    
      delete Nuc_px;
      delete Nuc_py;
      delete Nuc_pz;
      delete Nuc_E;
      
      delete t_Pi0;
      delete strip_Q2;
      delete strip_Xbj;
    
      
      delete Pi0;
      delete Ph1;
      delete Ph2;
      delete El;
      delete Nuc;
      delete Vertex;
      
      delete ch1;
      delete ineff2g;
      delete tree;
      out_sim->Close();
      delete out_sim;

      std::cerr << "Stopped execution due to event generation, set generate=true and try again for getting results" << std::endl;
      return BAerr;


    }

  TChain *chain= new TChain("pDVCS");
  for (int i = 1; i <= int(sizeof(Mbins) / sizeof(Mbins[0])) ; i++) {
    // Construct the filename for each file
    
    //TString filename = Form("Maxime_contamination/Maxime_pi0_%i.root",i);
    //TString filename_r = Form("Maxime_contamination/Maxime_pi0_recast_%i.root",i);
    TString filename = Maxime_bkg + Form("Maxime_pi0_%i.root",i);
    TString filename_r = Maxime_bkg + Form("Maxime_pi0_recast_%i.root",i);

    // Check if the file exists
    std::ifstream fileStream(filename);
    std::ifstream fileStream_r(filename_r);
    if (!fileStream.is_open() && !fileStream_r.is_open()) 
      {
	std::cerr << "Error: " << filename << " does not exist." << std::endl;
	std::cerr << "Please generate all Maxime_Pi0_i.root files before continuing" << std::endl;
	std::cerr << "You can do it by launching this Maxime method with the option generate=true" << std::endl;
	std::cerr << "Then, execute it again with generate=false, otherwise no BSA will be generated" << std::endl;

	return BAerr;
      }
    chain->AddFile(filename);

  }   

  if(recast)
    {
      TTree* Tree = chain ;
      TFile f(Maxime_bkg + Form("Maxime_pi0_recast_%i.root",bin_number),"RECREATE");
      Tree->SetMaxTreeSize(40000000000LL);        
      printf("... copying tree\n");
      TTree* sTree = Tree->CopyTree(cut);
      sTree->SetMaxTreeSize(40000000000LL);    
      printf("... tree copied ... \n");
      sTree->Write();
      f.Close();

      return BAerr;
    }
  
  if(add_BDT_Max)
    {
      Add_BDT_var_float(cut, Maxime_bkg + Form("Maxime_pi0_recast_%i.root",bin_number), TString("TMaxime_pi0.root"));
      //gSystem->Exec(TString("rm ") + Folder + TString("Maxime_pi0.root"));
    }
  
  TChain *Data= new TChain("pDVCS");
  Data->Add(Folder + TData);

  TChain *bkg= new TChain("pDVCS");
  bkg->Add(Folder + TString("TMaxime_pi0.root"));




  TCut cut_phi;
  double est1;
  double est1_FT;
  double est1_FD;
  double est2;
  double est2_FT;
  double est2_FD;
  std::ofstream outFile;
  std::ofstream outFile2;

  TString mvars[4] = {"Phi_Ph","t_Ph","strip_Q2","strip_Xbj"};
  double Lbound[4] = {0.  , boundaries.at(0), boundaries.at(2), boundaries.at(4)};
  double Ubound[4] = {360., boundaries.at(1), boundaries.at(3), boundaries.at(5)};

//Mean on the 3D bin
for(int j=1; j<4; j++)
  {
    TH1F *Data_BDT = new TH1F("Data_BDT","",1000,Lbound[j], Ubound[j]);
    Data_BDT->Sumw2();

    TH1F *Phi_BDT = new TH1F("Phi_BDT","",1000,Lbound[j], Ubound[j]);
    Phi_BDT->Sumw2();

    Data->Project("Data_BDT", mvars[j], cut + cut_phi + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)));
    bkg->Project("Phi_BDT", mvars[j], (cut + cut_phi + TCut(Form("_strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));

    if(Data_BDT->Integral() < Phi_BDT->Integral())
      {
	for(int k=1; k<=Data_BDT->GetNbinsX(); k++)
	  Phi_BDT->SetBinContent(k, Data_BDT->GetBinContent(k));
      }

    Data_BDT->Add(Data_BDT, Phi_BDT,1 ,-1);
    std::cout<<mvars[j]<<" (Maxime) mean is: "<<Data_BDT->GetMean()<<endl;
	  
    switch(j)
      {
      case 1:
	tmean2 = Data_BDT->GetMean();
	break;
      case 2:
	Qmean2 = Data_BDT->GetMean();
	break;
      case 3:
	xmean2 = Data_BDT->GetMean();
	break;
      }
    delete Phi_BDT;
    delete Data_BDT;

  }
//end of mean on the 3D bin

  if(means_maxi)
    {
      TString mvars[4] = {"Phi_Ph","t_Ph","strip_Q2","strip_Xbj"};
      double Lbound[4] = {0., -10., 0., 0.};
      double Ubound[4] = {360., 0., 10., 1.0};
      outFile.open(Folder + TString("means_maxi.txt"));
      outFile2.open(Folder + TString("contamination_maxi.txt"));
      outFile <<"mean_phi error mean_t error mean_Q2 error mean_xB error"<<endl;
      outFile2 <<"bin_number entries_bef entries_aft before after"<<endl;

      for(int i=0; i<Nphi; i++)
	{
	  cut_phi =TCut(Form("Phi_Ph > %f && Phi_Ph < %f", (360./Nphi)*i, (360./Nphi)*(i+1)));
	  for(int j=0; j<4; j++)
	    {

	      TH1F *Data_p = new TH1F("Data_p","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_m = new TH1F("Data_m","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",1000,Lbound[j], Ubound[j]);
	      TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",1000,Lbound[j], Ubound[j]);

	      Data_p->Sumw2();
	      Data_m->Sumw2();

	      Data->Project("Data_p", mvars[j], cut + cut_phi + TCut("Helicity>0"));
	      Data->Project("Data_m", mvars[j], cut + cut_phi + TCut("Helicity<0"));
	      Data->Project("Data_p_BDT", mvars[j], cut + cut_phi + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
	      Data->Project("Data_m_BDT", mvars[j], cut + cut_phi + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));

	      //****************To get estimation before BDT*****************
	      TH1F *Phi_p = new TH1F("Phi_p","",1000,Lbound[j], Ubound[j]);
	      TH1F *Phi_m = new TH1F("Phi_m","",1000,Lbound[j], Ubound[j]);
	      Phi_p->Sumw2();
	      Phi_m->Sumw2();

	      bkg->Project("Phi_p", mvars[j], (cut + cut_phi + TCut("Helicity>0"))*TCut("Weight"));
	      bkg->Project("Phi_m", mvars[j], (cut + cut_phi + TCut("Helicity<0"))*TCut("Weight"));
	      
	      if(Data_p->Integral() < Phi_p->Integral())
		{
		  for(int k=1; k<=Data_p->GetNbinsX(); k++)
		    Phi_p->SetBinContent(k, Data_p->GetBinContent(k));
		}

	      if(Data_m->Integral() < Phi_m->Integral())
		{
		  for(int k=1; k<=Data_m->GetNbinsX(); k++)
		    Phi_m->SetBinContent(k, Data_m->GetBinContent(k));
		}

	      if(Data_p->Integral() + Data_m->Integral() ==0)
	      	est1=0;
	      else
	      	est1=(Phi_p->Integral() + Phi_m->Integral())*1.0/(Data_p->Integral() + Data_m->Integral());  

	      delete Phi_p;
	      delete Phi_m;
	      //*************************************************************
	      TH1F *Phi_p_BDT = new TH1F("Phi_p_BDT","",1000,Lbound[j], Ubound[j]);
	      TH1F *Phi_m_BDT = new TH1F("Phi_m_BDT","",1000,Lbound[j], Ubound[j]);

	      Phi_p_BDT->Sumw2();
	      Phi_m_BDT->Sumw2();

	      bkg->Project("Phi_p_BDT", mvars[j], (cut + cut_phi + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
	      bkg->Project("Phi_m_BDT", mvars[j], (cut + cut_phi + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));

	      if(Data_p_BDT->Integral() < Phi_p_BDT->Integral())
		{
		  for(int k=1; k<=Data_p_BDT->GetNbinsX(); k++)
		    Phi_p_BDT->SetBinContent(k, Data_p_BDT->GetBinContent(k));
		}

	      if(Data_m_BDT->Integral() < Phi_m_BDT->Integral())
		{
		  for(int k=1; k<=Data_m_BDT->GetNbinsX(); k++)
		    Phi_m_BDT->SetBinContent(k, Data_m_BDT->GetBinContent(k));
		}

	      entries_bef_maxi=Data_p_BDT->Integral() + Data_m_BDT->Integral();
	      if(Data_p_BDT->Integral() + Data_m_BDT->Integral() ==0)
	      	est2=0;
	      else
	      	est2=(Phi_p_BDT->Integral() + Phi_m_BDT->Integral())*1.0/(Data_p_BDT->Integral() + Data_m_BDT->Integral());

	      Data_p_BDT->Add(Data_p_BDT, Phi_p_BDT,1 ,-1);
	      Data_m_BDT->Add(Data_m_BDT, Phi_m_BDT,1 ,-1);
	      
	      entries_aft_maxi=Data_p_BDT->Integral() + Data_m_BDT->Integral();

	      TH1F *mean_final= new TH1F("mean_final","mean_final",1000, Lbound[j], Ubound[j]);
	      mean_final->Add(Data_p_BDT, Data_m_BDT, 1, 1);  
	  
	      if(j==0)
		{
		  std::cout<<"\n Maxime estimation: Phi bin "<<i+1<<" : "<<est1*100<<"% "<<est2*100<<"%"<<endl;
		  outFile2<<i+1<<" "<<entries_bef_maxi<<" "<<entries_aft_maxi<<" "<<est1*100<<"% "<<est2*100<<"%"<<endl;

		  if(entries_aft_maxi>1)
		    est1=mean_final->GetMean();
		  else
		    est1=0;

		  std::cout<<mvars[j]<<" mean in bin "<<i+1<<" is: "<<est1<<endl;
		  outFile<<" "<<est1<<" "<<mean_final->GetStdDev();
		}
	      else
		{
		  if(entries_aft_maxi>1)
		    est1=mean_final->GetMean();
		  else
		    est1=0;

		  std::cout<<mvars[j]<<" mean in bin "<<i+1<<" is: "<<est1<<endl;
		  outFile<<" "<<est1<<" "<<mean_final->GetStdDev();
		}
	  
	      delete Phi_p_BDT;
	      delete Phi_m_BDT;
	      delete Data_p;
	      delete Data_m;
	      delete Data_p_BDT;
	      delete Data_m_BDT;
	      delete mean_final;

	    }
	  outFile<<" "<<endl;
	}
      outFile.close();
    }






  TH1F *Data_p = new TH1F("Data_p","",Nphi,0,360);
  TH1F *Data_m = new TH1F("Data_m","",Nphi,0,360);
  TH1F *Data_p_FT = new TH1F("Data_p_FT","",Nphi,0,360);
  TH1F *Data_m_FT = new TH1F("Data_m_FT","",Nphi,0,360);
  TH1F *Data_p_FD = new TH1F("Data_p_FD","",Nphi,0,360);
  TH1F *Data_m_FD = new TH1F("Data_m_FD","",Nphi,0,360);
  TH1F *Data_p_BDT = new TH1F("Data_p_BDT","",Nphi,0,360);
  TH1F *Data_m_BDT = new TH1F("Data_m_BDT","",Nphi,0,360);
  TH1F *Data_p_BDT_FT = new TH1F("Data_p_BDT_FT","",Nphi,0,360);
  TH1F *Data_m_BDT_FT = new TH1F("Data_m_BDT_FT","",Nphi,0,360);
  TH1F *Data_p_BDT_FD = new TH1F("Data_p_BDT_FD","",Nphi,0,360);
  TH1F *Data_m_BDT_FD = new TH1F("Data_m_BDT_FD","",Nphi,0,360);

  Data_p->Sumw2();
  Data_m->Sumw2();

  Data->Project("Data_p", "Phi_Ph", cut + TCut("Helicity>0"));
  Data->Project("Data_m", "Phi_Ph", cut + TCut("Helicity<0"));
  Data->Project("Data_p_FT", "Phi_Ph", cut + TCut("strip_Ph_Theta < 5 && Helicity>0"));
  Data->Project("Data_m_FT", "Phi_Ph", cut + TCut("strip_Ph_Theta < 5 && Helicity<0"));
  Data->Project("Data_p_FD", "Phi_Ph", cut + TCut("strip_Ph_Theta > 5 && Helicity>0"));
  Data->Project("Data_m_FD", "Phi_Ph", cut + TCut("strip_Ph_Theta > 5 && Helicity<0"));
  Data->Project("Data_p_BDT", "Phi_Ph", cut + TCut(Form("Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT", "Phi_Ph", cut + TCut(Form("Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_p_BDT_FT", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta < 5 && Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT_FT", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta < 5 && Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_p_BDT_FD", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta > 5 && Helicity>0 && _strip_Nuc_BDT > %f",BDT_cut)));
  Data->Project("Data_m_BDT_FD", "Phi_Ph", cut + TCut(Form("strip_Ph_Theta > 5 && Helicity<0 && _strip_Nuc_BDT > %f",BDT_cut)));

  //****************To get estimation before BDT*****************

  TH1F *Phi_p = new TH1F("Phi_p","",Nphi,0,360);
  TH1F *Phi_m = new TH1F("Phi_m","",Nphi,0,360);
  TH1F *Phi_p_FT = new TH1F("Phi_p_FT","",Nphi,0,360);
  TH1F *Phi_m_FT = new TH1F("Phi_m_FT","",Nphi,0,360);
  TH1F *Phi_p_FD = new TH1F("Phi_p_FD","",Nphi,0,360);
  TH1F *Phi_m_FD = new TH1F("Phi_m_FD","",Nphi,0,360);

  Phi_p->Sumw2();
  Phi_m->Sumw2();
  Phi_p_FT->Sumw2();
  Phi_m_FT->Sumw2();
  Phi_p_FD->Sumw2();
  Phi_m_FD->Sumw2();

  bkg->Project("Phi_p", "Phi_Ph", (cut + TCut("Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m", "Phi_Ph", (cut + TCut("Helicity < 0"))*TCut("Weight"));
  bkg->Project("Phi_p_FT", "Phi_Ph", (cut + TCut("strip_Ph_Theta < 5 && Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m_FT", "Phi_Ph", (cut + TCut("strip_Ph_Theta < 5 && Helicity < 0"))*TCut("Weight"));
  bkg->Project("Phi_p_FD", "Phi_Ph", (cut + TCut("strip_Ph_Theta > 5 && Helicity > 0"))*TCut("Weight"));
  bkg->Project("Phi_m_FD", "Phi_Ph", (cut + TCut("strip_Ph_Theta > 5 && Helicity < 0"))*TCut("Weight"));

  for(int k=1; k<=Data_p->GetNbinsX(); k++)
    {
      if(Data_p->GetBinContent(k) < Phi_p->GetBinContent(k))
	{
	  Phi_p->SetBinContent(k, Data_p->GetBinContent(k));
	    }
      if(Data_m->GetBinContent(k) < Phi_m->GetBinContent(k))
	{
	  Phi_m->SetBinContent(k, Data_m->GetBinContent(k));
	    }
      if(Data_p_FT->GetBinContent(k) < Phi_p_FT->GetBinContent(k))
	{
	  Phi_p_FT->SetBinContent(k, Data_p_FT->GetBinContent(k));
	    }
      if(Data_m_FT->GetBinContent(k) < Phi_m_FT->GetBinContent(k))
	{
	  Phi_m_FT->SetBinContent(k, Data_m_FT->GetBinContent(k));
	    }
      if(Data_p_FD->GetBinContent(k) < Phi_p_FD->GetBinContent(k))
	{
	  Phi_p_FD->SetBinContent(k, Data_p_FD->GetBinContent(k));
	    }
      if(Data_m_FD->GetBinContent(k) < Phi_m_FD->GetBinContent(k))
	{
	  Phi_m_FD->SetBinContent(k, Data_m_FD->GetBinContent(k));
	    }
    }

  if(Data_p->Integral() + Data_m->Integral() ==0)
 	est1=0;
  else
  	est1=(Phi_p->Integral() + Phi_m->Integral())*1.0/(Data_p->Integral() + Data_m->Integral());  

  if(Data_p_FT->Integral() + Data_m_FT->Integral() ==0)
 	est1_FT=0;
  else
  	est1_FT=(Phi_p_FT->Integral() + Phi_m_FT->Integral())*1.0/(Data_p_FT->Integral() + Data_m_FT->Integral());  

  if(Data_p_FD->Integral() + Data_m_FD->Integral() ==0)
 	est1_FD=0;
  else
  	est1_FD=(Phi_p_FD->Integral() + Phi_m_FD->Integral())*1.0/(Data_p_FD->Integral() + Data_m_FD->Integral());  


  delete Phi_p;
  delete Phi_m;
  delete Phi_p_FT;
  delete Phi_m_FT;
  delete Phi_p_FD;
  delete Phi_m_FD;
  delete Data_p_FT;
  delete Data_m_FT;
  delete Data_p_FD;
  delete Data_m_FD;
  //*************************************************************
  TH1F *Phi_p_BDT = new TH1F("Phi_p_BDT","",Nphi,0,360);
  TH1F *Phi_m_BDT = new TH1F("Phi_m_BDT","",Nphi,0,360);
  TH1F *Phi_p_BDT_FT = new TH1F("Phi_p_BDT_FT","",Nphi,0,360);
  TH1F *Phi_m_BDT_FT = new TH1F("Phi_m_BDT_FT","",Nphi,0,360);
  TH1F *Phi_p_BDT_FD = new TH1F("Phi_p_BDT_FD","",Nphi,0,360);
  TH1F *Phi_m_BDT_FD = new TH1F("Phi_m_BDT_FD","",Nphi,0,360);

  Phi_p_BDT->Sumw2();
  Phi_m_BDT->Sumw2();

  bkg->Project("Phi_p_BDT", "Phi_Ph", (cut + TCut(Form("Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT", "Phi_Ph", (cut + TCut(Form("Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_p_BDT_FT", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta < 5 && Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT_FT", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta < 5 && Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_p_BDT_FD", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta > 5 && Helicity > 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));
  bkg->Project("Phi_m_BDT_FD", "Phi_Ph", (cut + TCut(Form("strip_Ph_Theta > 5 && Helicity < 0 && _strip_Nuc_BDT > %f",BDT_cut)))*TCut("Weight"));

  for(int k=1; k<=Data_p_BDT->GetNbinsX(); k++)
    {
      if(Data_p_BDT->GetBinContent(k) < Phi_p_BDT->GetBinContent(k))
	{
	  Phi_p_BDT->SetBinContent(k, Data_p_BDT->GetBinContent(k));
	    }
      if(Data_m_BDT->GetBinContent(k) < Phi_m_BDT->GetBinContent(k))
	{
	  Phi_m_BDT->SetBinContent(k, Data_m_BDT->GetBinContent(k));
	    }
      if(Data_p_BDT_FT->GetBinContent(k) < Phi_p_BDT_FT->GetBinContent(k))
	{
	  Phi_p_BDT_FT->SetBinContent(k, Data_p_BDT_FT->GetBinContent(k));
	    }
      if(Data_m_BDT_FT->GetBinContent(k) < Phi_m_BDT_FT->GetBinContent(k))
	{
	  Phi_m_BDT_FT->SetBinContent(k, Data_m_BDT_FT->GetBinContent(k));
	    }
      if(Data_p_BDT_FD->GetBinContent(k) < Phi_p_BDT_FD->GetBinContent(k))
	{
	  Phi_p_BDT_FD->SetBinContent(k, Data_p_BDT_FD->GetBinContent(k));
	    }
      if(Data_m_BDT_FD->GetBinContent(k) < Phi_m_BDT_FD->GetBinContent(k))
	{
	  Phi_m_BDT_FD->SetBinContent(k, Data_m_BDT_FD->GetBinContent(k));
	    }
    }

  entries_bef_maxi=Data_p_BDT->Integral() + Data_m_BDT->Integral();
  entries_bef_maxi_FT=Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral();
  entries_bef_maxi_FD=Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral();
  
  if(Data_p_BDT->Integral() + Data_m_BDT->Integral() ==0)
  	est2=0;
  else
  	est2=(Phi_p_BDT->Integral() + Phi_m_BDT->Integral())*1.0/(Data_p_BDT->Integral() + Data_m_BDT->Integral());

  if(Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral() ==0)
  	est2_FT=0;
  else
  	est2_FT=(Phi_p_BDT_FT->Integral() + Phi_m_BDT_FT->Integral())*1.0/(Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral());

  if(Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral() ==0)
  	est2_FD=0;
  else
  	est2_FD=(Phi_p_BDT_FD->Integral() + Phi_m_BDT_FD->Integral())*1.0/(Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral());
  std::cout<<"Maxime estimation "<<est1<<" "<<est2<<endl;
  std::cout<<"Maxime estimation FT "<<est1_FT<<" "<<est2_FT<<endl;
  std::cout<<"Maxime estimation FD "<<est1_FD<<" "<<est2_FD<<endl;
  
  Data_p_BDT->Add(Data_p_BDT, Phi_p_BDT,1 ,-1);
  Data_m_BDT->Add(Data_m_BDT, Phi_m_BDT,1 ,-1);
  
  entries_aft_maxi=Data_p_BDT->Integral() + Data_m_BDT->Integral();
  entries_aft_maxi_FT=Data_p_BDT_FT->Integral() + Data_m_BDT_FT->Integral();
  entries_aft_maxi_FD=Data_p_BDT_FD->Integral() + Data_m_BDT_FD->Integral();

  std::cout<<"Overall "<<entries_bef_maxi<<" "<<entries_aft_maxi<<" "<<entries_bef_maxi_FT<<" "<<entries_bef_maxi_FD<<" (Total "<<entries_bef_maxi_FT+entries_bef_maxi_FD<<") "<<est1*100<<"% "<<est2*100<<"% "<<est1_FT*100<<"% "<<est2_FT*100<<"% "<<est1_FD*100<<"% "<<est2_FD*100<<"% "<<endl;
  boundaries.push_back(est1);
  boundaries.push_back(est2);
  boundaries.push_back(est1_FT);
  boundaries.push_back(est2_FT);
  boundaries.push_back(est1_FD);
  boundaries.push_back(est2_FD);

  if(means_maxi)
    {
      outFile2<<"Overall "<<entries_bef_maxi<<" "<<entries_aft_maxi<<" "<<est1*100<<"% "<<est2*100<<"%"<<endl;
      outFile2.close();
    }

  
  TH1 *BA= Data_m_BDT->GetAsymmetry(Data_p_BDT);
  TCanvas* c2 = new TCanvas("c2","Histograms");

  
  //The fit gets attached to BA, so if BA is deleted, everything is deleted
  //That is why I plot a new "fit" function
  TF1 *fitf = new TF1("fitf","[0]*sin(x*TMath::Pi()/180)/(1+[1]*cos(x*TMath::Pi()/180))",0,360);
  fitf->SetParameter(0,0.1);
  fitf->SetParameter(1,-0.3);
  fitf->SetParLimits(0,0,1.0);
  fitf->SetParLimits(1,-1.,1.);
  
  BA->SetAxisRange(-1., 1.,"Y");  
  BA->Scale(1.0/Bpol);
  //BA->Fit("fitf","Q");
  BA->SetTitle("RG-A note");


  // Create an output file to save the histogram
  TFile *outputFile = new TFile(Folder + TString("Maxime_Clean.root"), "RECREATE");
  fitf->SetLineColor(kBlack);
  BA->SetLineColor(kBlack);
  BA->SetMarkerColor(kBlack);

  // Write the histogram to the output file
  BA->Write();

  // Close the output file
  outputFile->Close();

  std::ofstream outFile6(Folder + TString("BSA_Maxi_Values.txt"));
  std::ofstream outFile7(Folder + TString("entries_maxi_p.txt"));
  std::ofstream outFile8(Folder + TString("entries_maxi_m.txt"));
  for(int k=1; k<=BA->GetNbinsX(); k++)
    {
      outFile6<<Data_p_BDT->GetBinContent(k) + Data_m_BDT->GetBinContent(k)<<", "<<BA->GetBinContent(k)<<", "<<BA->GetBinError(k)<<endl;
      outFile7<<Data_p_BDT->GetBinContent(k)<<endl;
      outFile8<<Data_m_BDT->GetBinContent(k)<<endl;
    }
  
  outFile6.close();
  outFile7.close();
  outFile8.close();


  fitf->SetLineColor(kRed);
  BA->SetLineColor(kRed);
  BA->SetMarkerColor(kRed);
  BA->SetAxisRange(-1.0, 1.0,"Y");
  BA->SetAxisRange( 0. ,360.,"X");
  BA->GetXaxis()->SetTitle("#phi(deg)");
  BA->GetYaxis()->SetTitle("BSA");
  BA->SetTitle("Method 2 estimation");
  BA->GetXaxis()->SetTitleSize(0.06);
  BA->GetYaxis()->SetTitleSize(0.06);
  BA->GetXaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetLabelSize(0.05);
  BA->GetYaxis()->SetTitleOffset(0.5);
  BA->GetXaxis()->SetTitleOffset(0.5);
  BA->GetYaxis()->SetNdivisions(6);
  BA->GetXaxis()->SetNdivisions(4);
  BA->Draw();

  c2->Print(Folder + TString("Maxime_BSA.pdf"));
  BSA_Amplitude_maxi=BA->GetBinContent(BA->GetMaximumBin());
  BSA_Amplitude_maxi_fit=fitf->GetParameter(0);
  BSA_Error_maxi_fit=fitf->GetParError(0);

  std::cout<<"Final amplitudes: value - fit "<<BSA_Amplitude_maxi<<" "<<BSA_Amplitude_maxi_fit<<endl;



  // Clean up memory
  delete outputFile;
      
  delete c2;
  delete Phi_p_BDT;
  delete Phi_m_BDT;
  delete Phi_p_BDT_FT;
  delete Phi_m_BDT_FT;
  delete Phi_p_BDT_FD;
  delete Phi_m_BDT_FD;
  delete Data_p;
  delete Data_m;
  delete Data_p_BDT;
  delete Data_m_BDT;
  delete Data_p_BDT_FT;
  delete Data_m_BDT_FT;
  delete Data_p_BDT_FD;
  delete Data_m_BDT_FD;

  delete Data;
  delete bkg;
  
  return BA;
  
}
