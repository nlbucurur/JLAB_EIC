void BDT::Get_BDT_Score()
{
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT");
   gStyle->SetOptStat(0);

     
   TFile *File = new TFile(Folder + TString("TMVACC.root"));   
   TTree *training = (TTree *)(File->Get("dataset/TrainTree"));

   TH1D *MVA_BDT_Train_S__4 = new TH1D("MVA_BDT_Train_S__4","Signal",50,-0.5, 0.5);
   
   TString bdtvar;
   if(categories)
   	bdtvar="Ph_Topology";
   else
   	bdtvar="BDT";
   
   training->Project("MVA_BDT_Train_S__4", bdtvar,"classID==0");
   MVA_BDT_Train_S__4->Scale(1.0/MVA_BDT_Train_S__4->Integral());


   MVA_BDT_Train_S__4->SetLineColor(kBlue);
   MVA_BDT_Train_S__4->SetMarkerStyle(20);
   MVA_BDT_Train_S__4->SetMarkerColor(kBlue);
   MVA_BDT_Train_S__4->SetMarkerSize(0.7);
   
   TH1D *MVA_BDT_Train_B__5 = new TH1D("MVA_BDT_Train_B__5","Background",50,-0.5, 0.5);
   training->Project("MVA_BDT_Train_B__5", bdtvar,"classID==1");
   //MVA_BDT_Train_B__5->Scale(1.0/MVA_BDT_Train_B__5->Integral());
   //MVA_BDT_Train_B__5->Add(MVA_BDT_Train_S__4);
   MVA_BDT_Train_B__5->Scale(1.0/MVA_BDT_Train_B__5->Integral());

   MVA_BDT_Train_B__5->SetLineColor(2);
   MVA_BDT_Train_B__5->SetMarkerColor(kRed);
   MVA_BDT_Train_B__5->SetMarkerStyle(20);
   MVA_BDT_Train_B__5->SetMarkerSize(0.7);


   //**************************************************************

   
   TFile *File_DVCS = new TFile(Folder + TString("Tested_Quality_Data.root"));   
   TTree *pDVCS = (TTree *)(File_DVCS->Get("pDVCS"));

   TH1D *hData = new TH1D("hData","Data",50,-0.5, 0.5);
   TCut TheCut =TCut("bestCandidateFlag==1 && strip_Xbj <1 && strip_Xbj >0 && t_Ph <0 && strip_Q2 > 1.0 && strip_W > 2 && strip_Nuc_P > 0.35 && strip_El_P > 1.0 && strip_Ph_P>2  && strip_El_vz < 10 && strip_El_vz > -12 && TMath::Abs(Phi_Nuc - Phi_Ph) < 2 && TMath::Abs(t_Nuc - t_Ph) < 2 && TMath::Sqrt(Xbal * Xbal + Ybal*Ybal + Zbal*Zbal) <1");
   pDVCS->Project("hData", "_strip_Nuc_BDT",TheCut);

   hData->SetMarkerColor(kBlack);
   hData->SetMarkerStyle(20);
   hData->SetMarkerSize(0.7);


   auto hs = new THStack("hs","");
   hs->Add(MVA_BDT_Train_S__4);
   hs->Add(MVA_BDT_Train_B__5);

   hData->Scale(1.0/hData->Integral());
   hData->SetLineColor(kBlack);
   hData->SetFillColor(kBlack);
   hData->SetFillStyle(3554);

   hs->SetHistogram(new TH1F("hstot","",50,-0.5, 0.5));
   hs->GetHistogram()->GetXaxis()->SetTitle("BDT Score");
   hs->GetHistogram()->GetYaxis()->SetTitle("counts/total events");
   hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
   hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
   hs->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
   hs->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
   hs->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
   hs->GetHistogram()->GetXaxis()->SetTitleOffset(0.8);
   hs->GetHistogram()->GetYaxis()->SetNdivisions(8);
   hs->GetHistogram()->GetXaxis()->SetNdivisions(8);
   hs->GetHistogram()->GetYaxis()->SetMaxDigits(2);


   hs->Draw("nostack");
   hData->Draw("hist same");
canvas1->BuildLegend();

/*      
   TLegend *leg = new TLegend(0.105,0.78,0.505,0.9);
   leg->AddEntry("MVA_BDT_S","Signal","F");
   leg->AddEntry("MVA_BDT_B","Background","F");
   leg->AddEntry("Data","Data","F");
   leg->Draw();
*/


  canvas1->Print(Folder + TString("BDT_Score.pdf"));
}
