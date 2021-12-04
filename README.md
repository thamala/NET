# NET
Network Enrichment Test, program used in Hämälä et al. 2020 MBE (https://doi.org/10.1093/molbev/msz206)

Compiling:  
gcc NET.c -o NET -lm

The program can be used for two types of analyses: finding co-expression modules that are overenriched among a set of candidate genes (e.g. *F*<sub>ST</sub> outliers) and finding modules that are enriched for high test scores (e.g. *F*<sub>ST</sub>).

Example for conducting the candidate gene test:  
./NET -m example.modules.txt -i example.candidates.txt -t 0 > candidate.output.txt

Modules are defined (-m) with a tab-delimited text file (UNIX line endings), where each line corresponds to a gene and the module of which the gene is part of. The input parameter (-i) is followed by a list of candidate genes (single gene per line). The choice of test (-t 0) defines that the program searches for enrichment among candidate genes. The program then counts genes belonging to each module and performs Fisher's exact tests to find overenriched modules. If the network contains less than 100 modules, chi-square test is performed to asses overall deviation from expected proportions.

Examples for conducting the whole-module enrichment test:  
./NET -m example.modules.txt -i example.allgenes.txt -t 1 -c 10000 > all.output.txt  
./NET -m example.modules.txt -i example.allgenes.txt -n example.neutral.txt -t 1 -c 1e6 > all.output.txt

Here, the input list (-i) is tab-delimited file, where the gene at each line is followed by a test score, such as *F*<sub>ST</sub>. In the first example, *P*-values are estimated by randomizing (-c) the observed estimates 10,000 times across modules. In the second example, the null distribution is constructed by sampling estimates from a list of neutral values (-n). This is repeated 1 million times (-c 1e6). By default, *P*-values are estimated as the proportion of repeats that produce the same or higher value than the observed one. This can be changed to lower than the observed by adding '-s 1' parameter.

Use './NET -h' to see all options.
