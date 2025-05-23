void BDT::Build_DVCS_Tree(TTree*& pDVCS_tree, TLorentzVector*& electron, TLorentzVector*& photon, TLorentzVector*& Nuc)
{
  double Pmass=0.938;
  double Ebeam=10.6;
  
  TLorentzVector *beam = new TLorentzVector();
  TLorentzVector *target = new TLorentzVector();
  beam->SetXYZT(0.0, 0.0, 10.6, 10.6);
  target->SetXYZT(0.0, 0.0, 0.0, Pmass);

  // for exclusivity variables
  strip_Q2 = 4 * 10.6 * electron->P() * TMath::Power(TMath::Sin(electron->Theta() / 2), 2.);
  strip_W = TMath::Sqrt(Pmass*Pmass + 2 * Pmass * (10.6 - electron->P()) - strip_Q2);
  strip_Xbj =  strip_Q2 / (2 * Pmass * (10.6 - electron->P())) ;

  strip_El_px = electron->Px();
  strip_El_py = electron->Py();
  strip_El_pz = electron->Pz();
  strip_El_E  = electron->E();
  strip_El_P  = electron->P();
  strip_El_Theta = electron->Theta()*180/TMath::Pi();
  strip_El_Phi = electron->Phi()*180/TMath::Pi();

  strip_Ph_px = photon->Px();
  strip_Ph_py = photon->Py();
  strip_Ph_pz = photon->Pz();
  strip_Ph_E  = photon->E();
  strip_Ph_P  = photon->P();
  strip_Ph_Theta = photon->Theta()*180/TMath::Pi();
  strip_Ph_Phi = photon->Phi()*180/TMath::Pi();

  strip_Nuc_px = Nuc->Px();
  strip_Nuc_py = Nuc->Py();
  strip_Nuc_pz = Nuc->Pz();
  strip_Nuc_E  = Nuc->E();
  strip_Nuc_P  = Nuc->P();
  strip_Nuc_Theta = Nuc->Theta()*180/TMath::Pi();
  strip_Nuc_Phi = Nuc->Phi()*180/TMath::Pi();
  
  VelectronIn = beam->Vect();
  VelectronOut = electron->Vect();
  VnucleonOut = Nuc->Vect();
  VphotonOut = photon->Vect();
  Vvirtualphoton = (*beam - *electron).Vect();
  Vlepto = VelectronIn.Cross(VelectronOut);
  Vhadro = VnucleonOut.Cross(Vvirtualphoton);
  VhadroPP = VnucleonOut.Cross(VphotonOut);
  PhV_Vec = (*beam - *electron);
  pp2 = (PhV_Vec.E() - (*photon).E() +2*Pmass)*(PhV_Vec.E() - (*photon).E());;
  cos2theta=pow((PhV_Vec.P()*PhV_Vec.P() +pp2 - (*photon).E()*(*photon).E()),2)/(4*pp2*pow(PhV_Vec.P(),2));
  cos2theta_exp = pow(TMath::Cos(VnucleonOut.Angle(Vvirtualphoton)),2);
  dcos2theta = cos2theta - cos2theta_exp;
  Ph_E_Th = 0.5*(2*Pmass*( beam->E() - electron->E()) - strip_Q2)/(beam->E() - electron->E() + Pmass - Vvirtualphoton.Mag()*TMath::Cos(VphotonOut.Angle(Vvirtualphoton)));
  deltaE = photon->E() - Ph_E_Th;

 
  Phi_Nuc = 180. / TMath::Pi() * Vlepto.Angle(Vhadro);
  Phi_Ph = 180. / TMath::Pi() * Vlepto.Angle(VhadroPP);
  if (Vlepto.Dot(VnucleonOut) > 0.)
    Phi_Nuc = 360. - Phi_Nuc;
  if (Vlepto.Dot(VphotonOut) < 0.)
    Phi_Ph = 360. - Phi_Ph;
  delta_Phi=Phi_Nuc - Phi_Ph;

  t_Nuc=(*Nuc - *target).M2();
  double ratio= strip_Q2/(2.0*Pmass*strip_Xbj);
  double cos=(-electron->Px()*photon->Px() - electron->Py()*photon->Py() + (Ebeam - electron->Pz())*photon->Pz())/(photon->P()*sqrt( pow(electron->P(),2) + Ebeam*Ebeam - 2*Ebeam*electron->P()*TMath::Cos(electron->Theta()) ) );
  t_Ph = -(strip_Q2*Pmass + (strip_Q2/strip_Xbj)*( ratio - cos*sqrt(strip_Q2 + pow(ratio,2)) ) )/(Pmass + ratio - cos*sqrt(strip_Q2 + pow(ratio,2)) );

  delta_t=t_Nuc - t_Ph;

  BalV = *beam + *target - *photon  - *electron - *Nuc;
  miss_mom_eNg = TMath::Sqrt(pow(BalV.X(), 2) + pow(BalV.Y(), 2) + pow(BalV.Z(), 2));
  Xbal = BalV.X();
  Ybal = BalV.Y();
  Zbal = BalV.Z();
  Ebal = BalV.T();
  p_perp = TMath::Sqrt(pow(BalV.X(), 2) + pow(BalV.Y(), 2));
  bestCandidateFlag = 1;
  
  mm2_eNg = (*beam + *target - *photon  - *electron - *Nuc).M2();
  mm2_eNg_N = (*beam + *target - *photon  - *electron - *Nuc).M2();
  mm2_eg = (*beam + *target - *photon  - *electron).M2();
  mm2_ep = (*beam + *target - *electron - *Nuc).M2();
  mm2_gp = (*beam + *target - *photon - *Nuc).M2();
  mm2_e = (*beam + *target- *electron).M2();
  mm2_g = (*beam + *target - *photon).M2();
  mm2_p = (*beam + *target - *Nuc).M2();

  theta_gamma_X = 180. / TMath::Pi() * (*beam + *target - *electron - *Nuc).Angle(photon->Vect());
  theta_gamma_e = 180. / TMath::Pi() * TMath::ACos((VphotonOut.Dot(VelectronOut)) / (VphotonOut.Mag() * VelectronOut.Mag()));
  theta_N_e = 180. / TMath::Pi() * TMath::ACos((VnucleonOut.Dot(VelectronOut)) / (VnucleonOut.Mag() * VelectronOut.Mag()));

  pDVCS_tree->Fill(); 

  delete beam;
  delete target;
}



