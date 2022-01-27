The following is a brief description of the files/folders contained within this folder


- 01_sex_effect_on_lifespan.R
	- generates figures:
		- Distribution_of_lifespan.pdf
		- Sex_surv_curve.pdf

- 02a_sex_effects_on_hematology.R
	- estimates the effect of sex on each hematology trait using a mixed effects model
	- p-value of sex effect obtained using both parametric boot strap and likelihood ratio test
	- general model form:	CBC Trait ~ (1|Generation) + (1|Date of Test) + (1|MouseID) + Age + Sex
	- Output saved as: sex_effects_on_hematology_traits.csv

- 02a_sex_effects_on_hematology.sh
	- runs 02a_sex_effects_on_hematology.R on JAX HPC

- 02b_age_effects_on_hematology.R
	- estimates the effect of age on each hematology trait using a mixed effects model
	- p-value of age effect obtained using both parametric boot strap and likelihood ratio test
	- general model form:	CBC Trait ~ (1|Generation) + (1|Date of Test) + (1|MouseID) + Age + Sex
	- Output saved as: age_effects_on_hematology_traits.csv

- 02b_age_effects_on_hematology.sh
	- runs 02b_age_effects_on_hematology.R on JAX HPC

- 02c_sexspecific_age_effects_on_hematology.R
	- estimates the age-sex interaction effect on each hematology trait using a mixed effects model
	- p-value of age-sex interaction effect obtained using both parametric boot strap and likelihood ratio test
	- general model form:	CBC Trait ~ (1|Generation) + (1|Date of Test) + (1|MouseID) + Age + Sex + Age:Sex
	- Output saved as: agesex_interaction_effects_on_hematology_traits.csv

- 02c_sexspecific_age_effects_on_hematology.sh
	- runs 02c_sexspecific_age_effects_on_hematology.R on JAX HPC

- 03_hema_effects_on_lifespan_by_timepoint.R
	- estimates the effect of each hematology trait at each timepoint on lifespan using a mixed effects model
	- p-value of effects are obtained using both parametric boot strap and likelihood ratio test
	- general model form:	Lifespan ~ (1|Generation) + Sex + CBC Trait
	- Output saved as: hematology_effects_on_lifespan.csv

- 03_hema_effects_on_lifespan_by_timepoint.sh
	- runs 03_hema_effects_on_lifespan_by_timepoint.R on JAX HPC

- hemaLS_coeff_plots.R
	- generates a plot of the regression coefficient (effect of CBC trait on lifespan) for each phenotype at each timepoint. The y-axis is the effect of the trait on lifespan (months of lifespan per standard deviations of trait)
		- HemaLS_Coefficients.pdf

- figures_phenotype_specific_aging_and_lifespan_plots.R
	- generates summary age, sex, and lifespan plots for each CBC trait.
	- figures saved in:
		- aging_and_lifespan_pheno_plots/
