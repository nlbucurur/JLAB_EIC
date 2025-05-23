void BDT::Write_Null(int bin_number)
{
      int Nphi=Nphibins[bin_number-1];
      std::ofstream outFile1(Folder + TString("entries_most_p.txt"));
      std::ofstream outFile2(Folder + TString("entries_most_m.txt"));
      std::ofstream outFile3(Folder + TString("entries_maxi_p.txt"));
      std::ofstream outFile4(Folder + TString("entries_maxi_m.txt"));

      std::ofstream outFile5(Folder + TString("means_most.txt"));
      std::ofstream outFile6(Folder + TString("means_maxi.txt"));

      std::ofstream outFile7(Folder + TString("Amplitudes.txt"));


	for(int i=0; i<Nphi; i++)
	{
	outFile1<<0<<endl;
	outFile2<<0<<endl;
	outFile3<<0<<endl;
	outFile4<<0<<endl;
	
	outFile5<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
	outFile6<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;

      outFile7<<"type value fit error on fit"<<endl;	      
      outFile7<<"Raw "<<"0 0 0"<<endl;	      
      outFile7<<"Mostafa "<<"0 0 0"<<endl;	      
      outFile7<<"Maxime "<<"0 0 0"<<endl;	      
      outFile7<<"Entries before BDT (All/FT/FD): "<<"0 0 0"<<endl;
      outFile7<<"Mostafa entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<"0 0 0 0 0% 0% 0% 0% 0% 0%"<<endl;  
      outFile7<<"Maxime  entries/estimation (bef/aft/bef_FT/aft_FT/bef_FD/aft_FD) on bin "<<bin_number<< " before/after: "<<"0 0 0 0 0% 0% 0% 0% 0% 0%"<<endl;  
      
	}
      outFile1.close();
      outFile2.close();
      outFile3.close();
      outFile4.close();

      outFile5.close();
      outFile6.close();
      outFile7.close();

}
