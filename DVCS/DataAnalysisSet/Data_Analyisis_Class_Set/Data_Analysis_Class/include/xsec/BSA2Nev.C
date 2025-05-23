TH1* BDT::BSA2Nev(TH1* BSA){

  TH1F *Nev = new TH1F("Nev","",BSA->GetNbinsX(),0,360);
  double bsa,err, nev;
  
for(int i=1; i<=BSA->GetNbinsX(); i++)
{
bsa=BSA->GetBinContent(i);
err=BSA->GetBinError(i);
if(err==0)
  nev=0;
else
  nev=(1-bsa*bsa)/pow(err,2);

Nev->SetBinContent(i, nev);
Nev->SetBinError(i,sqrt(nev));

}

Nev->SetName(TString("N_") +BSA->GetName());
Nev->SetTitle(TString("N_") +BSA->GetTitle());

return Nev;
}

