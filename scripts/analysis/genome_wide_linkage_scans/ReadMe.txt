The following is a brief description of the files/folders contained within this folder


- 01_haplotype_scans.R
	- performs haplotype (linkage) scans for lifespan and all the CBC traits at each timepoint conditioned on sex and generation
	- performs 3 more haplotype scans for lifespan:
		- conditioned on sex, generation, and 7-month RDW
		- conditioned on sex, generation, and 7-month HDW
		- conditioned on sex, generation, and 7-month NLR

- 01_haplotype_scans.sh
	- runs 01_haplotype_scans.R on JAX HPC

- haploscan_perm_chunks/
	- This folder contains the scripts that perform the permutation scans for the genome-wide haplotype associations.
		- 100 R-scripts that each perform 100 permutation scans for lifespan and all the CBC traits at each timepoint conditioned on sex and generation. Output saved in:
			- haploscan_perm_chunks/
		- 100 shell scripts for submitting the 100 permutation scans to the JAX HPC
		- 1 shell script for looping through the 100 shell scripts to submit all the permutation scans to the JAX HPC

- 02_compile_haploscan_perm_chunks.R
	- Compiles the 100 separate permutation scans (each of 100 permutations) into a one large permutation scan of 10,000 permutations. Output saved as:
		- LifespanHema_10000Perms.rds

- 03_peak_tables_from_haploscans.R
	- create peak tables (1.5 LOD-drop support intervals) for lifespan and the hematology traits
	- output saved as:
		- Lifespan_Peak_Table.csv
		- Aging_Biomarker_Hema_Peak_Table.csv (RDW, HDW, and NLR)
		- Other_Hema_Peak_Table.csv

- figures_plots_for_linkage_scans.R
	- generates plots for haplotyp association scans
		- saved in: genome_wide_linkage_scans/
