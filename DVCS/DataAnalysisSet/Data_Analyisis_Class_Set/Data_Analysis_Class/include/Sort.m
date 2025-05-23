$HistoryLength = 0

Lala = Import["Data_P.csv", "Data"];
Lala2 = Flatten[Lala];
Lele = Import["Data_NP.csv", "Data"];
Lele2 = Flatten[Lele];

Print["Number of events in P dataset"]
Length[Lala2]
Print["Number of events in NP dataset"]
Length[Lele2]
Print["Number of events in common"]
Length[Intersection[Lala2, Lele2]]
Print["Raw percentage of P in NP (relative to P)"]
Length[Intersection[Lala2, Lele2]]/Length[Lala2] //N
Print["Minimum percentage of P in NP (relative to P)"]
(Length[Intersection[Lala2, Lele2]]-Length[Lala2]*0.12)/Length[Lala2] //N

Print["Maximum Confirmed DVCS on NP dataset (relative to NP)"]
(Length[Intersection[Lala2, Lele2]])/Length[Lele2] //N
Print["Minimum confirmed DVCS on NP dataset (relative to NP)"]
(Length[Intersection[Lala2, Lele2]]-Length[Lala2]*0.12)/Length[Lele2] //N
Clear[]

Quit
