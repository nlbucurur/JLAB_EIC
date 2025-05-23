(* ::Package:: *)

Matrix = Import["/home/munoz/Documents/Important Codes/bin_mig/BM_Matrix.csv", "Data"];
Print["Round"]
InvM=Matrix;

vecIpMost = Flatten[Import["/home/munoz/Documents/Important Codes/bin_mig/entries_most_p.csv", "Data"]];
vecImMost = Flatten[Import["/home/munoz/Documents/Important Codes/bin_mig/entries_most_m.csv", "Data"]];
vecIpMaxi = Flatten[Import["/home/munoz/Documents/Important Codes/bin_mig/entries_maxi_p.csv", "Data"]];
vecImMaxi = Flatten[Import["/home/munoz/Documents/Important Codes/bin_mig/entries_maxi_m.csv", "Data"]];

AMost = ((vecIpMost - vecImMost)/vecIpMost + vecImMost)/. {Indeterminate -> 0.0, ComplexInfinity->0.0};
AMaxi = ((vecIpMaxi - vecImMaxi)/vecIpMaxi + vecImMaxi)/. {Indeterminate -> 0.0, ComplexInfinity->0.0};
(*
AMost = If[(vecIpMost + vecImMost)>0,(vecIpMost - vecImMost)/(vecIpMost + vecImMost),0.0];
AMaxi = If[(vecIpMaxi + vecImMaxi),(vecIpMaxi - vecImMaxi)/(vecIpMaxi + vecImMaxi),0.0];
*)
(*Remember that Helicities are flipped in the data*)
vecOMost=InvM . AMost;
vecOMaxi=InvM . AMaxi;

(*To obtain a BM factor between 0 and 1*)
vecOMost = vecOMost/AMost/. {Indeterminate -> 0.0, ComplexInfinity->0.0};
vecOMaxi = vecOMaxi/AMaxi/. {Indeterminate -> 0.0, ComplexInfinity->0.0};
Export["factor_most_mig.csv", vecOMost]
Export["factor_maxi_mig.csv", vecOMaxi]
Clear[]


Quit

