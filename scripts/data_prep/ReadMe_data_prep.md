# Samples excluded
We excluded all mice that died before 8 months (43 mice: 19 females, 24 males), missed their first CBC assay (an additional 6 mice: 3 females, 3 males), are missing their genotype data (an additional 26 mice: 19 females, 7 males), or are missing their date of death (an additional 3 mice: 1 female, 2 males). In total we excluded 79 mice (43 females and 36 males), leaving us with 521 mice (257 females and 264 males) for analysis.

<br>

# Description of hematology data

We used an Advia 2120 automated blood counter equipped with software specialized for mice to perform blood counts at three different timepoints: 7, 13, and 19 months

For each hematology assay, we analyzed 16 traits automatically reported by the automated blood counter, as well as the neutrophil to lymphocyte ratio (NLR) which we manually calculated from the automatically reported neutrophil and lymphocyte counts:  
&emsp;&emsp;NLR = # Neutrophils / # Lymphocytes  
In addition to NLR, the automatically reported trait hematocrit is also a derived measure, calculated from RBC count and mean red blood cell volume (MCV):  
&emsp;&emsp;Hct = (# RBC * MCV) / 10  
Additionally, in order to match the units of red cell distribution width (RDW, the coefficient of variation of red blood cell volume, a percent), we converted hemoglobin distribution width (HDW, the standard deviation of red blood cell hemoglobin concentration) from a standard deviation to a percent coefficient of variation by dividing by the cellular hemoglobin concentration mean (CHCM, mean hemoglobin concentration of red blood cells) and then multiplying by 100%:  
&emsp;&emsp;HDW[%] = (HDW[SD] / CHCM) * 100%

The table below lists the peripheral blood traits included in our analyis:                                                                              
| Name												| Description											       | Units                        |
| ------------------------------------------------- | ------------------------------------------------------------ | ---------------------------- |
| Red blood cell (RBC) count						| Concentration of RBC in peripheral blood				       | 1000000 cells per microliter |		
| Hematocrit										| Proportion of peripheral blood consisting of RBCs		       | percent 					  |						
| Hemoglobin										| Hemoglobin concentration in peripheral blood			       | grams per deciliter          |			
| Corpuscular hemoglobin (CH)						| Average hemoglobin content of RBCs					       | picograms                    |		
| Corpuscular hemoglobin concentration mean (CHCM)	| Average hemoglobin concentration in RBCs				       | grams per deciliter          |	
| Hemoglobin distribution width (HDW)				| Coefficient of variation of hemoglobin concentration in RBCs | percent                      |				
| Mean corpuscular volume (MCV)						| Average volume of RBCs							 	       | femtoliters                  |
| Red cell distribution width (RDW)					| Coefficient of variation of RBC volume				       | percent                      |	
| Platelet count									| Concentration of platelets in peripheral blood		       | 1000 cells per microliter    |															
| Mean platelet volume (MPV)						| Average volume of platelets							       | femtoliters                  |			
| Mean platelet mass (MPM)							| Average dry mass of platelets							       | picograms                    |			
| White blood cell (WBC) count						| Concentration of WBCs in peripheral blood				       | 1000 cells per microliter    |						
| Lymphocyte count									| Concentration of lymphocytes in peripheral blood		       | 1000 cells per microliter    |								
| Neutrophil count									| Concentration of neutrophils in peripheral blood		       | 1000 cells per microliter    |								
| Monocyte count									| Concentration of monocytes in peripheral blood		       | 1000 cells per microliter    |								
| Eosinophil count									| Concentration of eosinophils in peripheral blood		       | 1000 cells per microliter    |								
| Lymphocyte to neutrophil ratio (LNR)				| Ratio of lymphocyte count to neutrophil count			       | none (ratio)                 |		

<br>

## Sample counts by timepoint

 - 21 mice (9 females and 12 males) have hematology data only from 7 months 
 - 60 mice (26 females and 34 males) have hematology data from 7 and 13 months, but not 19 months. 
 - 439 mice (222 females and 217 males) have hematology data from all three timepoints. 
 - The remaining mouse, a female, has hematology data from 7 and 19 months, but its 13-month hematology data is missing.					

<br>

## Missing data

There are no missing values in the 7-month hematology data. In the 13-month hematology data, there are two samples (from males) that are missing their eosinophil counts. In the 18-month hematology data, there is one sample (a male) that is missing data on the WBC differential (lymphocyte, neutrophil, monocyte, and eosinophil counts), and there are an additional two samples (one female, one male) that are missing only their eosinophil counts.

<br>

## Adjustment for batch and clumping effects

Eosinophil count and all the platelet traits (platelet count, mean platelet volume, mean platelet mass) were adjusted for batch effects and clumping effects. A mixed effects linear model (Eq1) was used to estimate the batch and clumping effects on each trait. The estimated effects were subtracted from the trait values and then the resulting values were rescaled to have the same range as the original data.  
&emsp;&emsp;Eq1: Trait ~ (1|Timepoint) + (1|Sex) + (1|Mouse) + (1|Batch) + Clump Count  
 1) Subtract batch and clump effects  
    - Trait* = Trait - Batch Effect - Clump Effect * Clump Count  
 2) Then rescale
    - Traitadj = ((Trait* - min(Trait*))/(max(Trait*) - min(Trait*))) * (max(Trait) - min(Trait)) + (min(Trait))  

<br>

The rest of the phenotypes were adjusted only for batch effects.  
&emsp;&emsp;Eq2: Trait ~ (1|Timepoint) + (1|Sex) + (1|Mouse) + (1|Batch)  
 1) Subtract batch and clump effects
    - Trait* = Trait - Batch Effect  
	
 2) Then rescale
    - Traitadj = ((Trait* - min(Trait*))/(max(Trait*) - min(Trait*))) * (max(Trait) - min(Trait)) + (min(Trait))  


