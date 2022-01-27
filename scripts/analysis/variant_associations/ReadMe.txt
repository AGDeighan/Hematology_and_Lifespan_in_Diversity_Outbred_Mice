The following is a brief description of the files/folders contained within this folder


- 01_variant_association_scans.R
	- performs genome-wide variant association scans. Output saved as:
		- LifespanHema_GenomeWideVariantScans.rds

- 01_variant_association_scans.sh
	- runs 01_variant_association_scans.R on JAX HPC

- varscan_perm_chunks/
	- This folder contains the scripts that perform the permutation scans for the genome-wide variant associations.
		- 100 R-scripts that each perform 100 permutation scans for lifespan, HDW at 7 months, and RDW at 7 months. Output saved in:
			- varscan_perm_chunks/
		- 100 shell scripts for submitting the 100 permutation scans to the JAX HPC
		- 1 shell script for looping through the 100 shell scripts to submit all the permutation scans to the JAX HPC

- 02_compile_varscan_perm_chunks.R
	- Compiles the 100 separate permutation scans (each of 100 permutations) into a one large permutation scan of 10,000 permutations. Output saved as:
		- LifespanHema_10000Perms.rds

- figures_genomewide_varscans.R
	- generates figures:
		- GenomeWide_VarScan_Lifespan_CondSexGen.tiff
		- GenomeWide_VarScan_Lifespan_CondSexGenBlood.tiff
		- GenomeWide_VarScan_RDW_CondSexGen.tiff
		- GenomeWide_VarScan_HDW_CondSexGen.tiff

- figures_chr18_varasso_plots_w_genes
	- generates figures:
		- LifespanRDW_unconditioned_variant_scan.pdf
		- LifespanHDW_unconditioned_variant_scan.pdf

- linkage_plots_correlation_and_distance.R
	- generates figures:
		- Variant_correlation_and_distance_plot_CAST.pdf
		- Variant_correlation_and_distance_plot_129-WSB.pdf
