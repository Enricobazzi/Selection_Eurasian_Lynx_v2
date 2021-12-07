---
title: "Identify_Candidate_Loci"
author: "Enrico"
date: "23/11/2021"
output: html_document
editor_options:
  chunk_output_type: console
---

The strategy adopted to identify loci under selection is a Genome Environment Association (GEA) analysis.

## Redundancy Analysis (RDA)

Of the different methods reviewed to perform GEA, it appears that redundancy analysis (RDA) is one of the most reliable (Forester et al. 2018). RDA is a multivariate method for exploring association of different environmental predictors with genetic data. By exploring which SNPs more strongly correlate with the different costrained axes of RDA ordination, we can identify candidates for selection.

I followed a few tutorials made to guide through the different steps of GEA and RDA analysis. Links of the tutorials are:

https://popgen.nescent.org/2018-03-27_RDA_GEA.html
https://bookdown.org/hhwagner1/LandGenCourse_book/WE-11.html

See RDA.md for details on how the analysis was run and its outputs.

## Univariate GEA with BayPass

In order to reduce false positive rates and add a control for population structure (RDA does not control for it), I have run univariate analysis using the software BayPass v2.1. See BayPass.md for details on how the analysis was run and its outputs.

BayPass results, in the form of SNPs Bayes Factors, have been used to generate a set of candidate genomic windows, using the software GenWin. See GenWin.md for details on how the analysis was run and its outputs.

## Intersecting results to obtain candidate windows and SNPs

By intersecting our RDA results with the genomic windows identified by BayPass, we can point which genomics windows and which SNPs are under selection and run downstream analysis to answer different questions. What genes are being selected and what function do they have for the organism? What effect does the environment have on the genetic structure of the populations? Can we identify an environmental gradient that is directly related to the genetic structure of the populations?

```{bash}
# On genomics-a server
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

for var in bio2 bio5 bio6 bio8 bio13 jan_depth snow_days
 do
  echo "${var}"
  bedtools intersect -a ../RDA/rda_candidate_snps.bed -b ../GenWin/${var}_GenWin_windows_outliers.bed \
   > ${var}_intersect_candidate_snps.bed
  bedtools intersect -wa -a ../GenWin/${var}_GenWin_windows_outliers.bed -b ../RDA/rda_candidate_snps.bed | uniq \
   > ${var}_intersect_candidate_windows.bed
  nsnps=($(wc -l <${var}_intersect_candidate_snps.bed))
  nwindows=($(wc -l <${var}_intersect_candidate_windows.bed))
  winlength=($(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${var}_intersect_candidate_windows.bed))
  echo "${var} has a total of ${nsnps} candidate SNPs, a total of ${nwindows} candidate genomic windows \
  spanning a total of ${winlength} bps"
done
cat *_intersect_candidate_snps.bed | sort -k 1,1 -k2,2n | uniq > total_intersect_candidate_snps.bed
cat *_intersect_candidate_windows.bed | sort -k 1,1 -k2,2n | uniq | bedtools merge -i - > total_intersect_candidate_windows.bed
```
bio2 has a total of 244 candidate SNPs, a total of 35 candidate genomic windows spanning a total of 800000 bps
bio5 has a total of 226 candidate SNPs, a total of 43 candidate genomic windows spanning a total of 1040000 bps
bio6 has a total of 429 candidate SNPs, a total of 55 candidate genomic windows spanning a total of 1220000 bps
bio8 has a total of 234 candidate SNPs, a total of 25 candidate genomic windows spanning a total of 670000 bps
bio13 has a total of 79 candidate SNPs, a total of 34 candidate genomic windows spanning a total of 610000 bps
jan_depth has a total of 1109 candidate SNPs, a total of 82 candidate genomic windows spanning a total of 1870000 bps
snow_days has a total of 894 candidate SNPs, a total of 63 candidate genomic windows spanning a total of 1460000 bps

total candidate SNPs : 2637
total candidate genomic windows : 220
total length of candidate genomic windows : 6.63 Mbp
