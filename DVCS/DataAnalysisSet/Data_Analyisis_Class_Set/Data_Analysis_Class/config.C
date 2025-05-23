void BDT::config()
{
    Nbins=64;

	//Working directory
    dir="/work/clas12/jsalvg/pass2_RGA-Analysis/inb/P/";
	//External BSA from analysis note -> To compare with a previous analysis
    extBSA="/work/clas12/jsalvg/pass2_RGA-Analysis/Maxime_BSA/"; 
	//External xsec from analysis note -> To compare with a previous analysis
    extXSEC="/work/clas12/jsalvg/pass2_RGA-Analysis/inb/P/Theory_xsec/"; 
	//Folder to store the results
	Folder="Analysis/";
    gSystem->Exec(TString("mkdir -p ") + Folder);
    
	//Compute Q2, t and xB on each phi bin ?
    means_most=false;
    means_maxi=false;

    //To generate the contamination files. 
    //Only one of these can be set to true. 
    //Thus, you need to run the script 3 times before being able to get any physics results
    generate=false;		//Generate pi0 decays
    recast=false;   	//Recombine 1-gamma pi0 events into files for each bin
    add_BDT_Max=true;  //Add BDT variable to the contamination files.

	//Estimate contamination from method 1 ? It is not needed to be computed everytime
    generate_most=true;

	//Compute cross-section?
	xsection=false;
	
	//Data samples for training and experimental data
    DVCS="/lustre24/expphy/volatile/clas12/jsalvg/simulation/dvcsgen/inb/withRC/Quality_DVCS_Train.root";    
    Pi0 ="/lustre24/expphy/volatile/clas12/jsalvg/simulation/aaogen/P/inb/train/Quality_Pi_as_DVCS_Train.root";
    //Pi0 ="/work/clas12/jsalvg/RGA-Simulation/pass2_inb/Pi0/Quality_Pi_as_DVCS_Train.root";
    Data="/work/clas12/jsalvg/Data/DVCS/pass2_RGA/inb/Quality2_Data_P.root"; //Original, keep after testing //gcorrected
    TData="Tested_Quality_Data.root";
    TDVCS="Tested_DVCS.root";
    TPi0="Tested_Pi0.root";

	//For eta Background estimation
    epeta = "/work/clas12/jsalvg/Data/eta/pass2_RGA/inb/Quality_Data_epeta_P.root";	
    sim_eta_as_dvcs = "/lustre24/expphy/volatile/clas12/jsalvg/simulation/eta/inb/P/1gamma/Quality_1gamma_eta.root";
    sim_epeta="/lustre24/expphy/volatile/clas12/jsalvg/simulation/eta/inb/P/2gamma/Quality_Sim_epeta_P.root";
    
	//For Pi0 Background subtraction
    eppi0_name = "Quality_Data_eppi0_P.root";	
    eppi0 = "/work/clas12/jsalvg/Data/DVMP/pass2_RGA/inb/P/" + eppi0_name;	
    maps_path = "../../maps/";
    sim_eppi0="/lustre24/expphy/volatile/clas12/jsalvg/simulation/aaogen/P/inb/bkg_sub/2gamma/Quality_Sim_eppi0_P.root";
    sim_eppi0_1="/lustre24/expphy/volatile/clas12/jsalvg/simulation/aaogen/P/inb/bkg_sub/2gamma/Quality_Sim_eppi0_P.root";
    sim_pi_as_dvcs = "/lustre24/expphy/volatile/clas12/jsalvg/simulation/aaogen/P/inb/bkg_sub/1gamma/Quality_Pi_as_DVCS_P_1.root";
    sim_pi_as_dvcs_1 = "/lustre24/expphy/volatile/clas12/jsalvg/simulation/aaogen/P/inb/bkg_sub/1gamma/Quality_Pi_as_DVCS_P_1.root";

	//RC effects dvcs sample
	RC_Sim=DVCS;
	//RC_Sim="/volatile/clas12/jsalvg/simulation/dvcsgen/inb/Quality_BM_DVCS.root"; 
    TRC_Sim="Tested_BM_Sim.root";
	//The MC version of the RC sample
	MC_RC_Sim="/volatile/clas12/jsalvg/simulation/dvcsgen/inb/MCgen/Quality_MC_DVCS.root"; 

    //Acc and BM corrections. Need a dvcs sample with RC effects. 
    //RC_can be used if it has enough statistics in all bins.
	BM_Sim=DVCS;
	//BM_Sim="/volatile/clas12/jsalvg/simulation/dvcsgen/inb/Quality_BM_DVCS.root"; 
    TBM_Sim="Tested_BM_Sim.root";
	//The MC version of the BM sample
	MC_BM_Sim="/volatile/clas12/jsalvg/simulation/dvcsgen/inb/MCgen/Quality_MC_DVCS.root"; 

	//Direct output from dvcsgen to compute the RC correction factor
	MC_DVCS_RC="/work/clas12/jsalvg/RGA-Simulation/dvcsgen_RC_gen.root"; 

	//Directory to temporaly store the contamination files from method 2
	Maxime_bkg = "/lustre24/expphy/volatile/clas12/jsalvg/DVCS_analysis/inb/P/";

	//Basic selection/exclusivity cuts
	cut="bestCandidateFlag==1 && \
    strip_Xbj <1 && strip_Xbj >0 && t_Ph <0 && strip_Q2 > 1.0 && \
    strip_W > 2 && strip_El_P > 1.0 && strip_Ph_P>2  && strip_El_vz < 10 && strip_El_vz > -12 && \
    theta_gamma_e > 5 && abs(delta_t)<2 && abs(delta_Phi)%180 < 2 && TMath::Sqrt(Xbal * Xbal + Ybal*Ybal + Zbal*Zbal) <1 && abs(mm2_ep)<0.5";
  
	//Cut to estimate systematic error due to cuts
    cut_sys="bestCandidateFlag==1 &&\
strip_Xbj <1 && strip_Xbj >0 && t_Ph <0 && strip_Q2 > 1.0 && \
strip_W > 2 && strip_El_P > 1.0 && strip_Ph_P>2  && strip_El_vz < 10 && strip_El_vz > -12 && \
theta_gamma_e > 5 && abs(delta_t)<1.5 && abs(delta_Phi)%180 < 1.5 && TMath::Sqrt(Xbal * Xbal + Ybal*Ybal + Zbal*Zbal) <0.8";


	//Cut to estimate systematic error due pid selection
    cut_pid="strip_El_vz > -8 && strip_El_vz<5 && \
    strip_Nuc_vz > -8 && strip_Nuc_vz<5 && \
    strip_Ph_beta>0.9 && strip_Ph_beta<1.1 && \
    (strip_El_P<4.5 || (strip_El_P>4.5 && strip_El_ECin_energy/strip_El_P > 0.2 - strip_El_PCAL_energy/strip_El_P))";
    //cut_pid="abs(strip_El_chi2pid)<3 && strip_El_vz > -8 && strip_El_vz<5 && \
    abs(strip_Nuc_chi2pid)<3 && strip_Nuc_vz > -8 && strip_Nuc_vz<5 && \
    strip_Ph_beta>0.9 && strip_Ph_beta<1.1 && \
    (strip_El_P<4.5 || (strip_El_P>4.5 && strip_El_ECin_energy/strip_El_P > 0.2 - strip_El_PCAL_energy/strip_El_P))";

	//Refinemnet cuts, if needed after BDT classification.
    //cut_ref="delta_t<0.1 && miss_mom_eNg<0.4 && theta_gamma_X < 0.6";
    //cut_ref="bestCandidateFlag==1 && theta_gamma_X<1.0";
    cut_ref="bestCandidateFlag==1";

    //Training variables
    Vars.push_back(TString("mm2_eNg"));
    Vars.push_back(TString("mm2_eg"));
    Vars.push_back(TString("delta_Phi"));
    Vars.push_back(TString("delta_t"));
    Vars.push_back(TString("theta_gamma_X"));

	//Beam polarization
    Bpol=0.86;
	//BDT cut
    BDT_value=0.0;
	//BDT cut for systematic error
    BDT_value_sys=0.04;
    
	beam->SetXYZT(0.0, 0.0, 10.6, 10.6);
	target->SetXYZT(0.0, 0.0, 0.0, 0.938);

	//|t|>|tmin| cut. It is good to have it coded
	//(t_Ph<-(strip_Q2*0.938 + (strip_Q2/strip_Xbj)*( (strip_Q2/(2*0.938*strip_Xbj)) - sqrt(strip_Q2 + pow((strip_Q2/(2*0.938*strip_Xbj)),2))))/(0.938 + (strip_Q2/(2*0.938*strip_Xbj)) - sqrt(strip_Q2 + pow((strip_Q2/(2*0.938*strip_Xbj)),2))))


  Mbins[0]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[1]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[2]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[3]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[4]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[5]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[6]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[7]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[8]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[9]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[10]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[11]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1 && strip_Q2 < 1.4 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");

  Mbins[12]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[13]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[14]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[15]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[16]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[17]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[18]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[19]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[20]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");
  Mbins[21]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0    && strip_Xbj < 0.13");
  Mbins[22]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.13 && strip_Xbj < 0.21");
  Mbins[23]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.4 && strip_Q2 < 1.8 && strip_Xbj > 0.21 && strip_Xbj < 1.0 ");

  Mbins[24]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  Mbins[25]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  Mbins[26]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  Mbins[27]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  Mbins[28]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  Mbins[29]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  Mbins[30]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  Mbins[31]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  Mbins[32]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");
  Mbins[33]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0    && strip_Xbj < 0.16");
  Mbins[34]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.16 && strip_Xbj < 0.26");
  Mbins[35]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 1.8 && strip_Q2 < 2.4 && strip_Xbj > 0.26 && strip_Xbj < 1.0 ");

  Mbins[36]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  Mbins[37]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  Mbins[38]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  Mbins[39]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  Mbins[40]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  Mbins[41]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  Mbins[42]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  Mbins[43]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  Mbins[44]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  Mbins[45]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0    && strip_Xbj < 0.21");
  Mbins[46]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.21 && strip_Xbj < 0.33");
  Mbins[47]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 2.4 && strip_Q2 < 3.25 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");

  Mbins[48]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33");
  Mbins[49]= TCut("bestCandidateFlag==1 && t_Ph<0.0    && t_Ph>-0.2 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  Mbins[50]= TCut("bestCandidateFlag==1 && t_Ph<0.0  && t_Ph>-0.4 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33 ");
  Mbins[51]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0  ");
  Mbins[52]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.   && strip_Xbj < 0.33");
  Mbins[53]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");
  Mbins[54]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0    && strip_Xbj < 0.33");
  Mbins[55]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 3.25 && strip_Q2 < 5.0 && strip_Xbj > 0.33 && strip_Xbj < 1.0 ");

  Mbins[56]= TCut("bestCandidateFlag==1 && t_Ph<0.0 && t_Ph>-0.2 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  Mbins[57]= TCut("bestCandidateFlag==1 && t_Ph<0.0 && t_Ph>-0.2 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  Mbins[58]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  Mbins[59]= TCut("bestCandidateFlag==1 && t_Ph<-0.2 && t_Ph>-0.4 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  Mbins[60]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  Mbins[61]= TCut("bestCandidateFlag==1 && t_Ph<-0.4 && t_Ph>-0.8 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  Mbins[62]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0    && strip_Xbj < 0.55");
  Mbins[63]= TCut("bestCandidateFlag==1 && t_Ph<-0.8 && t_Ph>-100 && strip_Q2 > 5.0 && strip_Q2 < 15.0 && strip_Xbj > 0.55 && strip_Xbj < 1.0 ");
  
  		bins ={{
  			{-0.2, 0, 1, 1.4, 0, 0.13},
			{-0.2, 0, 1, 1.4, 0.13, 0.21},
			{-0.2, 0, 1, 1.4, 0.21, 1},
			{-0.4, -0.2, 1, 1.4, 0, 0.13},
			{-0.4, -0.2, 1, 1.4, 0.13, 0.21},
			{-0.4, -0.2, 1, 1.4, 0.21, 1},
			{-0.8, -0.4, 1, 1.4, 0, 0.13},
			{-0.8, -0.4, 1, 1.4, 0.13, 0.21},
			{-0.8, -0.4, 1, 1.4, 0.21, 1},
			{-100, -0.8, 1, 1.4, 0, 0.13},
			{-100, -0.8, 1, 1.4, 0.13, 0.21},
			{-100, -0.8, 1, 1.4, 0.21, 1},
			{-0.2, 0, 1.4, 1.8, 0, 0.13},
			{-0.2, 0, 1.4, 1.8, 0.13, 0.21},
			{-0.2, 0, 1.4, 1.8, 0.21, 1},
			{-0.4, -0.2, 1.4, 1.8, 0, 0.13},
			{-0.4, -0.2, 1.4, 1.8, 0.13, 0.21},
			{-0.4, -0.2, 1.4, 1.8, 0.21, 1},
			{-0.8, -0.4, 1.4, 1.8, 0, 0.13},
			{-0.8, -0.4, 1.4, 1.8, 0.13, 0.21},
			{-0.8, -0.4, 1.4, 1.8, 0.21, 1},
			{-100, -0.8, 1.4, 1.8, 0, 0.13},
			{-100, -0.8, 1.4, 1.8, 0.13, 0.21},
			{-100, -0.8, 1.4, 1.8, 0.21, 1},
			{-0.2, 0, 1.8, 2.4, 0, 0.16},
			{-0.2, 0, 1.8, 2.4, 0.16, 0.26},
			{-0.2, 0, 1.8, 2.4, 0.26, 1},
			{-0.4, -0.2, 1.8, 2.4, 0, 0.16},
			{-0.4, -0.2, 1.8, 2.4, 0.16, 0.26},
			{-0.4, -0.2, 1.8, 2.4, 0.26, 1},
			{-0.8, -0.4, 1.8, 2.4, 0, 0.16},
			{-0.8, -0.4, 1.8, 2.4, 0.16, 0.26},
			{-0.8, -0.4, 1.8, 2.4, 0.26, 1},
			{-100, -0.8, 1.8, 2.4, 0, 0.16},
			{-100, -0.8, 1.8, 2.4, 0.16, 0.26},
			{-100, -0.8, 1.8, 2.4, 0.26, 1},
			{-0.2, 0, 2.4, 3.25, 0, 0.21},
			{-0.2, 0, 2.4, 3.25, 0.21, 0.33},
			{-0.2, 0, 2.4, 3.25, 0.33, 1},
			{-0.4, -0.2, 2.4, 3.25, 0, 0.21},
			{-0.4, -0.2, 2.4, 3.25, 0.21, 0.33},
			{-0.4, -0.2, 2.4, 3.25, 0.33, 1},
			{-0.8, -0.4, 2.4, 3.25, 0, 0.21},
			{-0.8, -0.4, 2.4, 3.25, 0.21, 0.33},
			{-0.8, -0.4, 2.4, 3.25, 0.33, 1},
			{-100, -0.8, 2.4, 3.25, 0, 0.21},
			{-100, -0.8, 2.4, 3.25, 0.21, 0.33},
			{-100, -0.8, 2.4, 3.25, 0.33, 1},
			{-0.2, 0, 3.25, 5, 0, 0.33},
			{-0.2, 0, 3.25, 5, 0.33, 1},
			{-0.4, -0.2, 3.25, 5, 0, 0.33},
			{-0.4, -0.2, 3.25, 5, 0.33, 1},			
			{-0.8, -0.4, 3.25, 5, 0, 0.33},
			{-0.8, -0.4, 3.25, 5, 0.33, 1},
			{-100, -0.8, 3.25, 5, 0, 0.33},
			{-100, -0.8, 3.25, 5, 0.33, 1},
			{-0.2, 0, 5, 15, 0, 0.55},
			{-0.2, 0, 5, 15, 0.55, 1},
			{-0.4, -0.2, 5, 15, 0, 0.55},
			{-0.4, -0.2, 5, 15, 0.55, 1},
			{-0.8, -0.4, 5, 15, 0, 0.55},
			{-0.8, -0.4, 5, 15, 0.55, 1},
			{-100, -0.8, 5, 15, 0, 0.55},
			{-100, -0.8, 5, 15, 0.55, 1}
  			}};  
	Nphibins = {
		      29,
		      22,
		      18,
		      30,
		      23,
		      13,
		      27,
		      26,
		      9,
		      14,
		      12,
		      12,
		      23,
		      24,
		      21,
		      24,
		      23,
		      15,
		      17,
		      22,
		      11,
		      18,
		      14,
		      23,
		      24,
		      25,
		      23,
		      24,
		      28,
		      15,
		      19,
		      27,
		      12,
		      27,
		      17,
		      10,
		      21,
		      21,
		      18,
		      23,
		      25,
		      15,
		      19,
		      22,
		      12,
		      13,
		      16,
		      10,
		      19,
		      19,
		      23,
		      18,
		      21,
		      17,
		      15,
		      13,
		      12,
		      12,
		      13,
		      12,
		      13,
		      12,
		      13,
		      18};  

  std::cout<<"Configuration loaded !"<<endl;		      
		      
}

