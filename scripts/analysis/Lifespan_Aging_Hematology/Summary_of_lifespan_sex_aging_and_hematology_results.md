# Sex and lifespan
The lifespans of our 521 mice (257 females and 264 males) are approximately normally distributed with a mean of 27.5 months (median 27.9) and a standard deviation of 8.9 months. The 25th and 75th percentiles are 21.6 and 33.4 months, and the maximum lifespan is 55.6 months (a male mouse).

The mean lifespan of females (28.0 months) is one month longer than the mean of males (27.0 months), but the difference in not statistically significant (t-test p-value = 0.154; F-test conditioned on generation p-value = 0.155). Moreover, the five longest living mice are all males.

<br>

# Age, sex, hematology, and lifespan

We estimated the age and sex effects on the blood traits, and the effects of the blood traits on lifespan using mixed effect models. Blood traits that are consistently associated with lifespan and change with age in a manner concordant with their effect on lifespan may be useful indicators of underlying processes that affect aging.

## Effects of blood traits on lifespan

The table below lists the effect of each trait at each timepoint on lifespan. We adjusted for multiple testing (Holm FWER, 17 tests) separately at each timepoint. For traits that did not have a significant effect (FWER >= 0.05) on lifespan at a given timepoint, the effect on lifespan is listed as "NA"

| Name												| 7 months         | 13 months        | 19 months        |
| ------------------------------------------------- | ---------------- | ---------------- | ---------------- |
| Red blood cell (RBC) count						| NA               | NA               | Longer lifespan  |		
| Hematocrit										| NA               | NA               | Longer lifespan  |						
| Hemoglobin										| NA               | Longer lifespan  | Longer lifespan  |			
| Corpuscular hemoglobin (CH)						| NA               | NA               | NA               |		
| Corpuscular hemoglobin concentration mean (CHCM)	| NA               | NA               | Longer lifespan  |	
| Hemoglobin distribution width (HDW)				| Shorter lifespan | Shorter lifespan | Shorter lifespan |				
| Mean corpuscular volume (MCV)						| NA               | NA               | NA               |
| Red cell distribution width (RDW)					| Shorter lifespan | Shorter lifespan | Shorter lifespan |	
| Platelet count									| NA               | NA               | NA               |															
| Mean platelet volume (MPV)						| NA               | NA               | NA               |			
| Mean platelet mass (MPM)							| NA               | Shorter lifespan | Shorter lifespan |			
| White blood cell (WBC) count						| NA               | NA               | NA               |						
| Lymphocyte count									| NA               | NA               | NA               |								
| Neutrophil count									| NA               | Shorter lifespan | Shorter lifespan |								
| Monocyte count									| NA               | NA               | NA               |								
| Eosinophil count									| NA               | NA               | NA               |								
| Neutrophil to lymphocte ratio (NLR)				| Shorter lifespan | Shorter lifespan | Shorter lifespan |	


Many traits are associated with lifespan when measured later in life (19 months). However, the only traits to affect lifespan at all three timepoints are RDW, HDW, and NLR. Higher levels of these traits at any of the timepoints are associated with shorter lifespans.

<br>

## Sex and age effects on blood traits


The table below lists the direct and interaction effects of sex and age on each of the blood traits. For testing statistical significance we used and alpha level of 0.05 and the Holm adjustment (17 tests for each trait) for family-wise error rate (FWER). Effects with an FWER >= 0.05 are listed as "NA"

| Name												| Age effect (direct) | Sex effect (direct) | Age-sex interaction effect                                    |
| ------------------------------------------------- | ------------------- | ------------------- | ------------------------------------------------------------- |
| Red blood cell (RBC) count						| NA                  | NA                  | NA                                                            |		
| Hematocrit										| Decrease            | NA                  | Decreases more rapidly in males                               |						
| Hemoglobin										| Decrease            | Lower in males      | NA                                                            |			
| Corpuscular hemoglobin (CH)						| Decrease            | NA                  | Decreases more rapidly in males                               |		
| Corpuscular hemoglobin concentration mean (CHCM)	| NA                  | Lower in males      | NA                                                            |	
| Hemoglobin distribution width (HDW)				| Increase            | Lower in males      | NA                                                            |				
| Mean corpuscular volume (MCV)						| Decrease            | NA                  | Decreases more rapidly in males                               |
| Red cell distribution width (RDW)					| Increase            | Lower in males      | NA                                                            |	
| Platelet count									| Increase            | Higher in males     | Increases more rapidly in males (little/no change in females) |															
| Mean platelet volume (MPV)						| Decrease            | Higher in males     | NA                                                            |			
| Mean platelet mass (MPM)							| NA                  | NA                  | Decreases more rapidly in males                               |			
| White blood cell (WBC) count						| NA                  | Lower in males      | NA                                                            |						
| Lymphocyte count									| Decrease            | NA                  | NA                                                            |								
| Neutrophil count									| Increase            | Higher in males     | Increases more rapidly in females                             |								
| Monocyte count									| Increase            | Higher in males     | NA                                                            |								
| Eosinophil count									| NA                  | NA                  | NA                                                            |								
| Neutrophil to lymphocte ratio (NLR)				| Increase            | Higher in males     | NA                                                            |		


The majority of the blood traits change signficantly with age. Of the traits for which the age effect differs significantly by sex (i.e. a significant age-sex interaction effect), the majority differ such that the age-related changes occur more rapidly in males than females. The excpetion to this pattern is neutrophil count, which increases more rapidly in females than males; however, on average males have higher neutrophil counts than females (direct sex effect).

Notably, RDW, HDW, and NLR (the three traits that are associated with lifespan at every timepoint) all increase with age, which is consistent with their effects on lifespan. These traits may serve as proxy measures of underlying processes that influence aging.

<br>
