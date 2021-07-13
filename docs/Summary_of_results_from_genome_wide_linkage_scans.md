# Genome-wide linkage analysis
We performed genome-wide linkage (haplotype) scans for each of the blood traits at each timepoint, as well as lifespan. For these scans sex and generation were modeled as additive fixed-effect covariates and a kinship matrix was included as an additive random effect to account for genetic relatedness (additive gnetic relatedness) across the population. Additionally, three lifespan scans were performed for which the 7-month value of either RDW, HDW, or NLR were also included as a covariate. 

To determine thresholds for statistically significant QTL while accounting for the multiple testing burden, we used permutation analysis (Churchill & Doerge, 1994). Specifically, we permuted the genotypes, breaking the potential associations at QTL while maintaining the effects of sex, generation, and kinship. For lifespan and each blood trait at each timepoint, we performed 10,000 permutations and estimated the 95% permutation threshold for identifying significant QTL and the 80% permutation threshold for identifying suggestive QTL. We did not perform permuations for the lifespan scans conditioned on RDW, HDW, or NLR.

For significant or suggestive QTL, we calculate the 1.5-LOD drop support interval (Dupuis & Siegmund 1999).

## Lifespan peaks

We did not identify any significant (LOD > 95% threshold) QTL for lifespan, but we did identify 3 suggestive (LOD > 80% threshold) QTL. The p-values in the table below were calculated using the observed LOD score and the maximum LOD scores from the 10,000 permutations

| Chromosome | Position (Mb) | Support Interval (Mb)   | P-value | LOD  | LOD (RDW) | LOD (HDW) | LOD (NLR) |
| -----------| ------------- | ----------------------- | ------- | ---- | --------- | --------- | --------- |  
|  2         | 128.764153    | 128.372768 - 129.364448 | 0.0648  | 7.22 | 7.35      | 6.64      | 7.26      | 	
|  6         |  25.22833     |  23.581197 -  27.652608 | 0.0904  | 7.02 | 6.67      | 6.67      | 6.47      |					
| 18         |  10.777436    |   8.523683 -  14.397343 | 0.1853  | 6.59 | 5.29      | 4.68      | 6.84      |		

The suggestive lifepan QTL on chromosome 18 drops by more than 1 LOD when conditioned on either RDW (LOD drop of 1.30) or HDW (LOD drop of 1.91) but not NLR (LOD increases by 0.25). Moreover RDW and HDW (see below) both have significant QTL whose 1.5-LOD drop support intervals overlap the support interval of the suggestive lifespan QTL.

<br>

## Peaks of aging hematology traits: RDW, HDW, and NLR

The table below shows the significant and suggestive QTL for RDW, HDW, and NLR.


| Trait | Timepoint | Chromosome | Position (Mb) | Support Interval (Mb)   | P-value | LOD   |
| ----- | --------- | -----------| ------------- | ----------------------- | ------- | ----- |
| HDW   | 7 months  | 7          | 105.602197    | 101.661608 - 105.713902 | < 1e-04 | 20.41 |  	
| HDW   | 13 months | 7          | 102.389926    |  99.623551 - 109.078324 | 0.0002  | 10.31 |  				
| HDW   | 19 months | 7          | 102.76635     |  98.351788 - 107.372019 | 0.0115  | 8.21  |  
| HDW   | 7 months  | 9          | 73.877038     |  69.492256 - 78.334238  | 0.0043  | 8.69  |  	
| HDW   | 7 months  | 9          | 108.092239    | 107.582543 - 108.344858 | < 1e-04 | 13.43 |  				
| HDW   | 13 months | 9          | 74.903813     |  71.485485 - 78.243724  | 0.0153  | 8.10  |  
| HDW   | 13 months | 9          | 108.064434    | 107.582543 - 108.321865 | < 1e-04 | 15.15 |  	
| HDW   | 19 months | 9          | 108.144158    | 107.606558 - 108.321865 | 0.0001  | 11.73 |  				
| NLR   | 7 months  | 13         | 34.308266     |  20.139385 - 35.088922  | 0.1879  | 6.60  |  
| HDW   | 7 months  | 18         | 17.063425     |  13.462651 - 24.624215  | 0.0960  | 6.99  |  	
| HDW   | 13 months | 18         | 16.121207     |  13.405937 - 23.773718  | 0.0320  | 7.64  |  				
| RDW   | 7 months  | 18         | 15.831359     |  14.053427 - 18.914412  | 0.0436  | 7.49  |  				
| RDW   | 19 months | 18         | 17.645947     |  13.462651 - 20.605282  | 0.0248  | 7.81  | 

 <br>

The table below is a simplified version of the full peak table shown above. We merged (taking lowest lowerbound of the support intervals and highest upper bound of the support intervals) QTL that overlapped for the same trait measured at different timepoints.

| Trait | Timepoints            | Chromosome | Position (Mb) | Support Interval (Mb)   | Min P-value | Max LOD   | Overlap Lifespan QTL |
| ----- | --------------------- | -----------| ------------- | ----------------------- | ----------- | --------- | -------------------- |
| HDW   | 7, 13, and 19 months  | 7          | 105.602197    |  98.351788 - 109.078324 | < 1e-04     | 20.41     | No  |  
| HDW   | 7 and 13 months       | 9          | 73.877038     |  69.492256 - 78.334238  | 0.0043      | 8.69      | No  |  
| HDW   | 7, 13, and 19 months  | 9          | 108.092239    | 107.582543 - 108.344858 | < 1e-04     | 15.15     | No  |  	
| NLR   | 7 months              | 13         | 34.308266     |  20.139385 - 35.088922  | 0.1879      | 6.60      | No  |  
| HDW   | 7 and 13 months       | 18         | 17.063425     |  13.405937 - 24.624215  | 0.0320      | 7.64      | Yes | 
| RDW   | 7 and 19 months       | 18         | 15.831359     |  13.462651 - 20.605282  | 0.0248      | 7.81      | Yes | 

RDW and HDW both have significant QTL on chromosome 18 that overlap the suggestive lifespan QTL on chromosome 18. Moreover, as shown above, the LOD score for lifespan at the suggestive QTL drops notably when conditioned on the 7-month measurement of either RDW or HDW. This could be a locus that effects lifespan through biological pathways that influence RDW and HDW. To verify this, we need to look further into the haplotype effects and the SNP associations at the locus.

<br>



## Peaks of other hematology traits

We will not go into detail regarding the results of the linkage analysis for the other traits, since our interest lies in using the aging hematology traits (RDW, HDW, and NLR) to improve our ability to detect lifespan QTL. 

The table below shows the significant and suggestive QTL for the other blood traits.

| Trait          | Timepoint | Chromosome | Position (Mb) | Support Interval (Mb)  | P-value | LOD   |
| -------------- | --------- | -----------| ------------- | ---------------------- | ------- | ----- |
| Eos. count     | 13 months |  1         | 158.231393   | 156.497611 - 166.883277 | 0.0760  |  7.14 |  	
| Neut. count    | 19 months |  1         |  78.649814   |  76.523368 - 84.231746  | 0.1103  |  6.93 |  				
| WBC count      |  7 months |  1         | 143.739911   | 135.816613 - 147.08021  | 0.0775  |  7.15 |  
| CHCM           | 13 months |  2         | 165.498634   | 165.289259 - 165.880893 | 0.0772  |  7.13 |  	
| Hgb            |  7 months |  2         | 110.263784   | 108.384206 - 168.928408 | 0.1865  |  6.58 |  				
| MPV            |  7 months |  3         |  57.159313   |  53.040565 - 63.039987  | 0.0001  | 11.00 |  
| MPV            | 13 months |  3         |  61.393825   |  53.040565 - 63.039987  | 0.1606  |  6.67 |  	
| MPV            | 19 months |  3         |  61.393825   |  56.339548 - 63.039987  | 0.0013  |  9.43 |  				
| Hct            |  7 months |  4         |  58.603055   |  46.537864 - 63.474008  | 0.1339  |  6.81 |  
| Neut. count    |  7 months |  4         |  30.786898   |  24.59641  - 33.871097  | 0.1033  |  6.97 |  	
| Neut. count    |  7 months |  4         | 138.46331    | 133.650139 - 140.291392 | 0.1918  |  6.58 |  				
| Platelet count | 19 months |  4         |  28.234346   |  17.671082 - 29.210592  | 0.1380  |  6.82 |  				
| WBC count      | 13 months |  4         | 132.955036   | 132.879727 - 134.543155 | 0.0925  |  7.00 | 			
| MCV            | 13 months |  5         | 137.58353    | 103.33775  - 139.304843 | 0.0794  |  7.16 | 			
| CHCM           |  7 months |  7         | 104.23904    | 102.533667 - 105.602197 | < 1e-04 | 37.87 | 
| CHCM           | 13 months |  7         | 104.148971   | 102.389926 - 105.602197 | < 1e-04 | 31.71 |  	
| CHCM           | 19 months |  7         | 104.23904    | 101.775441 - 105.705868 | < 1e-04 | 22.84 |  					
| CH             | 13 months |  8         |  12.56355    |  10.738292 - 13.476821  | 0.1075  |  6.93 |  
| CHCM           |  7 months |  9         | 107.576242   | 106.563624 - 108.519018 | < 1e-04 | 14.33 |  	
| CHCM           | 13 months |  9         | 111.177662   | 107.278677 - 111.484667 | 0.0011  |  9.72 |  				
| CHCM           | 19 months |  9         | 111.43643    | 108.973633 - 111.562813 | < 1e-04 | 11.60 |  		
| MCV            |  7 months |  9         | 108.944904   | 108.015342 - 110.917445 | 0.0007  |  9.45 | 
| RBC count      |  7 months |  9         |   6.25458    |   3.369866 - 17.898502  | 0.0314  |  7.64 |  	
| RBC count      |  7 months |  9         | 111.884934   | 110.505368 - 113.016242 | 0.0184  |  7.97 |  				
| RBC count      | 13 months |  9         |  21.497315   |  15.17283  - 23.20059   | 0.0422  |  7.45 |  
| CH             | 13 months | 10         |  65.768251   |  58.908569 - 68.867113  | 0.0664  |  7.21 |  	
| MCV            | 13 months | 10         |  61.989393   |  52.886204 - 69.259139  | 0.1574  |  6.71 |  				
| CH             |  7 months | 11         |  88.933542   |  82.146142 - 97.069663  | 0.0413  |  7.47 |  
| CH             |  7 months | 11         | 111.296906   | 109.296395 - 112.921042 | 0.1749  |  6.58 |  	
| CH             | 13 months | 11         |  91.749319   |  83.656832 - 97.147277  | 0.0285  |  7.67 |  				
| CH             | 13 months | 11         | 111.296906   | 108.172117 - 112.462401 | 0.1145  |  6.90 |  
| CH             | 19 months | 11         |  92.20418    |  88.921751 - 97.027016  | 0.0124  |  8.25 |  	
| MCV            |  7 months | 11         |  88.933542   |  87.857447 - 92.74347   | 0.0466  |  7.39 |  				
| MCV            |  7 months | 11         | 111.296906   | 111.241337 - 112.424492 | 0.0034  |  8.70 |  				
| MCV            | 13 months | 11         | 111.296906   |  80.937399 - 111.610202 | 0.1448  |  6.76 | 			
| mpm            |  7 months | 11         |  21.663668   |  18.533106 - 40.837269  | 0.1067  |  6.90 | 			
| WBC count      | 13 months | 12         |  77.763943   |  77.023767 - 83.04458   | 0.1291  |  6.82 | 	
| Platelet count |  7 months | 16         |  14.783585   |  13.47913  - 16.41683   | 0.0336  |  7.62 |  				
| RBC count      | 19 months | 16         |  96.122784   |  93.074097 - 96.618853  | 0.1151  |  6.93 |  
| CHCM           | 13 months | 17         |  66.177265   |  65.396598 - 67.409284  | 0.0662  |  7.22 |  	
| WBC count      |  7 months | 17         |  32.297558   |  31.329496 - 43.50301   | 0.0100  |  8.29 |  			
| Eos. count     | 19 months | 18         |  10.842164   |  10.63157  - 12.402095  | 0.1951  |  6.60 |  				
| Lymph. count   | 13 months | 18         |  15.834131   |  13.462651 - 17.063425  | 0.0887  |  7.06 |  			
| Mono. count    |  7 months | 19         |  10.167257   |   7.347888 - 10.320975  | 0.0083  |  8.38 |  				
| CH             |  7 months |  X         | 142.880525   | 141.241721 - 153.843307 | 0.0712  |  7.15 | 			
| Neut. count    | 13 months |  X         |  95.903369   |  93.901906 - 99.154293  | 0.1978  |  6.55 | 			
| WBC count      |  7 months |  X         | 140.105449   | 136.599742 - 140.938315 | 0.1893  |  6.62 | 

<br>
