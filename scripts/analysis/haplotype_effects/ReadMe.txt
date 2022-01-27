The following is a brief description of the files/folders contained within this folder


- 01_sex_effect_on_lifespan.R
	- performs best linear unbiased prediction (BLUP) scans across chromosome 18 for estimating allele effects on lifespan, RDW at 7 months, and HDW at 7 months. Output saved as:
		- Chr18_BLUPs_Lifespan.rds
		- Chr18_BLUPs_RDW07.rds
		- Chr18_BLUPs_HDW07.rds
	- Generates full chromosome BLUP-scan plot:
		- full_chr18_BLUP_scans_LS_RDW_HDW.pdf

- 02_fit1_coeffs_and_correlations.R
	- Generates plots of the allele effects on lifespan (x-axis) and RDW or HDW (y-axis) at each of the peaks on chromosome 18
		- Lifespan peak: allele_coefs_at_lifespan_locus_peak.pdf
		- 7-month RDW peak: allele_coefs_at_RDW07_locus_peak.pdf
		- 19-month RDW peak: allele_coefs_at_RDW19_locus_peak.pdf
		- 7-month HDW peak: allele_coefs_at_HDW07_locus_peak.pdf
		- 13-month HDW peak: allele_coefs_at_HDW13_locus_peak.pdf

