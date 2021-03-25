# fcfdr_manuscript
Repo to reproduce results from flexible cFDR manuscript

## simulations/

The simulations/ directory contains code to reproduce the simulation analysis, including code to generate Figure 1 and Figure 2.  
Please see seperate README file in this directory for further details. 

## GenoCanyon_analysis/

The GenoCanyon_analysis/ directory contains code to perform IHW, Boca and Leek's FDR regression and flexible cFDR analysis leveraging GenoCanyon scores.  
Please see the GenoWAP software (https://github.com/rlpowles/GenoWAP-V1.2) for instructions to run GenoWAP. 

## generate_H3K27ac_data/

The generate_H3K27ac_data/ directory contains code to download the raw H3K27ac fold change values, convert them to .bed format,  
map the values to the SNPs in the asthma analysis and generate q1 and q2 (the auxiliary data vectors to iterate over using flexible cFDR). 
See generate_H3K27ac_data/makeH3K27acdata_code for details. 

## H3K27ac_cFDR/

The H3K27ac_cFDR/ directory contains code to perform the flexible cFDR analysis leveraging H3K27ac fold change values.  
Please see the FINDOR software (https://github.com/gkichaev/FINDOR) for instructions to run FINDOR. 


