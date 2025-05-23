rm entries.txt
rm Mostafa.txt
rm Maxime.txt
rm amplitudes_most.txt
rm amplitudes_maxi.txt
rm amplitudes_full_most.txt
rm amplitudes_full_maxi.txt
rm entries_most_p.csv
rm entries_most_m.csv
rm entries_maxi_p.csv
rm entries_maxi_m.csv
rm Mostafa_eta.txt

touch entries.txt
touch Mostafa.txt
touch Maxime.txt
touch amplitudes_most.txt
touch amplitudes_maxi.txt
touch amplitudes_full_most.txt
touch amplitudes_full_maxi.txt
touch entries_most_p.csv
touch entries_most_m.csv
touch entries_maxi_p.csv
touch entries_maxi_m.csv
touch Mostafa_eta.txt

for ((i=1; i<=64; i++))
do
        tail -n 3 bin_$i/Amplitudes.txt | head -n 1 >> entries.txt 
        tail -n 2 bin_$i/Amplitudes.txt | head -n 1 >> Mostafa.txt 
        tail -n 1 bin_$i/Amplitudes.txt  >> Maxime.txt 

         tail -n 5 bin_$i/Amplitudes.txt | head -n 1 >> amplitudes_most.txt 
         tail -n 4 bin_$i/Amplitudes.txt | head -n 1 >> amplitudes_maxi.txt 

         cat bin_$i/BSA_Most_Values.txt >> amplitudes_full_most.txt 
         cat bin_$i/BSA_Maxi_Values.txt >> amplitudes_full_maxi.txt 

         cat bin_$i/entries_most_p.txt >> entries_most_p.csv 
         cat bin_$i/entries_most_m.txt >> entries_most_m.csv 
         cat bin_$i/entries_maxi_p.txt >> entries_maxi_p.csv 
         cat bin_$i/entries_maxi_m.txt >> entries_maxi_m.csv 
         
         cat bin_$i/eta_contamination.txt >> Mostafa_eta.txt
done


rm -r All
mkdir -p All
cd All
mkdir -p To_bin_or_not_to_bin
mkdir -p Effect_of_background_subtraction
mkdir -p RGA_vs_This
mkdir -p BDT
mkdir -p Excl_vars
mkdir -p Train_vars
mkdir -p Kin_vars_bef_BDT
mkdir -p Kin_vars_aft_BDT
mkdir -p eppi0
cd ..

for ((i=1; i<=64; i++))
do
    # Checking if the bin_i directory exists
    if [ -d "bin_$i" ]; then
        cp bin_$i/BSA_Global_vs_Bin.png All/To_bin_or_not_to_bin/Bin_$i.png
        cp bin_$i/Background_subtraction.pdf All/Effect_of_background_subtraction/bin_$i.pdf
        cp bin_$i/RGA_vs_This.pdf All/RGA_vs_This/bin_$i.pdf                
        cp bin_$i/BDT_Score.pdf All/BDT/bin_$i.pdf
        cp bin_$i/Excl_vars.pdf All/Excl_vars/bin_$i.pdf
        cp bin_$i/Training_vars.pdf All/Train_vars/bin_$i.pdf
        cp bin_$i/Kin_Vars_bef_BDT.pdf All/Kin_vars_bef_BDT/bin_$i.pdf
        cp bin_$i/Kin_Vars_aft_BDT.pdf All/Kin_vars_aft_BDT/bin_$i.pdf
        cp bin_$i/eppi0_Data_vs_MC.pdf All/eppi0/bin_$i.pdf
    fi
done

