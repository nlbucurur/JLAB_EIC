TVector3 BDT::PCAL_PointFromLabToLocal(TVector3 Lab, int nsector){
  TVector3 Loc(0,0,0);
  Loc.SetXYZ(Lab.X(),Lab.Y(),Lab.Z());
  Double_t phi=(nsector-1)*60*TMath::DegToRad();
  Loc.RotateZ(-phi);
  Loc.SetXYZ(Loc.X()-xO,Loc.Y()-yO,Loc.Z()-zO);
  Double_t xloc=Loc.Y();
  Double_t yloc=TMath::Cos(beta)*Loc.X()-TMath::Sin(beta)*Loc.Z();
  Double_t zloc=TMath::Sin(beta)*Loc.X()+TMath::Cos(beta)*Loc.Z();
  Loc.SetXYZ(xloc,yloc,zloc);
  return Loc;
}


TVector3 BDT::PCAL_VectorFromLabToLocal(TVector3 Lab, int nsector){
  TVector3 Loc(0,0,0);
  Loc.SetXYZ(Lab.X(),Lab.Y(),Lab.Z());
  Double_t phi=(nsector-1)*60*TMath::DegToRad();
  Loc.RotateZ(-phi);
  Loc.SetXYZ(Loc.X(),Loc.Y(),Loc.Z());
  Double_t xloc=Loc.Y();
  Double_t yloc=TMath::Cos(beta)*Loc.X()-TMath::Sin(beta)*Loc.Z();
  Double_t zloc=TMath::Sin(beta)*Loc.X()+TMath::Cos(beta)*Loc.Z();
  Loc.SetXYZ(xloc,yloc,zloc);
  return Loc;
}


TVector3 BDT::FromLocalXYZtoUVW(TVector3 XYZ){
  double u=(XYZ.Y()-yoU)/TMath::Sin(alpha);
  double v=XYZ.X()-xoV-(yh-XYZ.Y())/TMath::Tan(alpha);
  double w=-XYZ.X()+xoW+(XYZ.Y()-yh)/TMath::Tan(alpha);
  return TVector3(u,v,w);
}

///////////////////////////////////////////////// Intersection with Calorimter planes
void BDT::GetIntersectionUVW(TVector3 vertex,TVector3 gamma, TLorentzVector& intersection){
  intersection.SetPxPyPzE(0,0,0,0);
  double theta=acos(gamma.Z()/gamma.Mag());
  double phi=atan2(gamma.Y(),gamma.X());

  if (theta*TMath::RadToDeg()<4.5){//Then it is FTCal
    double radius=TMath::Tan(theta)*(z_FTCAL-vertex.Z());
    intersection.SetPxPyPzE(radius*TMath::Cos(phi)+vertex.X(),radius*TMath::Sin(phi)+vertex.Y(),z_FTCAL,-666);   
  }
  else{
    double phi_temp=phi*TMath::RadToDeg();
    if (phi_temp<-30) phi_temp=phi_temp+360;
    int nsector=(phi_temp+30)/60+1;
    TVector3 vertex_loc=PCAL_PointFromLabToLocal(vertex, nsector);
    TVector3 gamma_loc=PCAL_VectorFromLabToLocal(gamma, nsector);
    double thetaloc=acos(gamma_loc.Z()/gamma_loc.Mag());
    double philoc=atan2(gamma_loc.Y(),gamma_loc.X());
    double radius=TMath::Tan(thetaloc)*(-vertex_loc.Z());
    TVector3 intersectionloc(radius*TMath::Cos(philoc)+vertex_loc.X(),radius*TMath::Sin(philoc)+vertex_loc.Y(),0);
    TVector3 intersectionXYZ=FromLocalXYZtoUVW(intersectionloc);
    intersection.SetPxPyPzE(intersectionXYZ.X(),intersectionXYZ.Y(),intersectionXYZ.Z(),0);
  }
  
}


/// //////////////////////////////////////////////////////////////////////////////////////////////
/// 1. Fiducial cuts for the PCAL (required for electrons and photons)from RGA S.Diehl
///

/// 1.1 A homgenous cut for all sectors which does not consider variations of the sampling fraction
///     This cut assumes, that for the electron ID only a proper cluster formation is required, which
///     is given, if the center of the cluster has enough distance to the edges.
///     loose: cuts 1 bar from the outer side (4.5 cm)
///     medium: cuts 2 bars from the outer side (9.0 cm)
///     tight: cuts 3 bars from the outer side (13.5 cm)

bool BDT::EC_hit_position_fiducial_cut_homogeneous(TLorentzVector uvw)
{

  ///////////////////////////
  bool tight = true;
  bool medium = false;
  bool loose = false;

  // Cut using the natural directions of the scintillator bars/ fibers:

  if (uvw.X() == 0 && uvw.Y() == 0 && uvw.Z() == 0)
    return true;

  double u = uvw.X();
  double v = uvw.Y();
  double w = uvw.Z();

  /// v + w is going from the side to the back end of the PCAL, u is going from side to side
  /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.

  ///////////////////////////////////////////////////////////////////
  //
  double min_u_tight_inb[] = {39.0, 39.0, 39.0, 39.0, 39.0, 39.0};
  double min_u_med_inb[] = {29.0, 29.0, 29.0, 29.0, 29.0, 29.0};
  double min_u_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
  //
  double max_u_tight_inb[] = {400, 400, 400, 400, 400, 400};
  double max_u_med_inb[] = {420, 420, 420, 420, 420, 420};
  double max_u_loose_inb[] = {420, 420, 420, 420, 420, 420};
  //
  double min_v_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
  double min_v_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
  double min_v_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
  //
  double max_v_tight_inb[] = {420, 420, 420, 420, 420, 420};
  double max_v_med_inb[] = {420, 420, 420, 420, 420, 420};
  double max_v_loose_inb[] = {420, 420, 420, 420, 420, 420};
  //
  double min_w_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
  double min_w_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
  double min_w_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
  //
  double max_w_tight_inb[] = {420, 420, 420, 420, 420, 420};
  double max_w_med_inb[] = {420, 420, 420, 420, 420, 420};
  double max_w_loose_inb[] = {420, 420, 420, 420, 420, 420};

  double min_u = 0;
  double max_u = 0;
  double min_v = 0;
  double max_v = 0;
  double min_w = 0;
  double max_w = 0;

  for (Int_t k = 0; k < 1; k++) //Do not need to find the sector...cuts are the same
    {
	  if (tight == true)
            {
	      min_u = min_u_tight_inb[k];
	      max_u = max_u_tight_inb[k];
	      min_v = min_v_tight_inb[k];
	      max_v = max_v_tight_inb[k];
	      min_w = min_w_tight_inb[k];
	      max_w = max_w_tight_inb[k];
            }
	  if (medium == true)
            {
	      min_u = min_u_med_inb[k];
	      max_u = max_u_med_inb[k];
	      min_v = min_v_med_inb[k];
	      max_v = max_v_med_inb[k];
	      min_w = min_w_med_inb[k];
	      max_w = max_w_med_inb[k];
            }
	  if (loose == true)
            {
	      min_u = min_u_loose_inb[k];
	      max_u = max_u_loose_inb[k];
	      min_v = min_v_loose_inb[k];
	      max_v = max_v_loose_inb[k];
	      min_w = min_w_loose_inb[k];
	      max_w = max_w_loose_inb[k];
            }
    }

  bool edges = v > min_v && v < max_v && w > min_w && w < max_w;

  if(edges)
    return true;
  else
    return false;
}


bool BDT::Fiducial_cut(TVector3 vertex,TVector3 gamma)
{
if(gamma.Theta()*TMath::RadToDeg() > 5)
{
  TLorentzVector uvw_vector;
  GetIntersectionUVW(vertex, gamma, uvw_vector);
  return EC_hit_position_fiducial_cut_homogeneous(uvw_vector);
  }
  else
  {
  return true;
  }
}
