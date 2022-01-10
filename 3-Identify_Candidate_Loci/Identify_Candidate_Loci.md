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
bio2 has a total of 1696 candidate SNPs, a total of 326 candidate genomic windows spanning a total of 7240000 bps
bio5 has a total of 1754 candidate SNPs, a total of 334 candidate genomic windows spanning a total of 7720000 bps
bio6 has a total of 2193 candidate SNPs, a total of 369 candidate genomic windows spanning a total of 8210000 bps
bio8 has a total of 1567 candidate SNPs, a total of 307 candidate genomic windows spanning a total of 6740000 bps
bio13 has a total of 824 candidate SNPs, a total of 274 candidate genomic windows spanning a total of 6610000 bps
jan_depth has a total of 4540 candidate SNPs, a total of 622 candidate genomic windows spanning a total of 13950000 bps
snow_days has a total of 3518 candidate SNPs, a total of 505 candidate genomic windows spanning a total of 10910000 bps

total candidate SNPs : 10791
total candidate genomic windows : 1620
total length of candidate genomic windows : 49.19 Mbp

To get the TOPSNP of each candidate window:
```{bash}
# on genomics-a:
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

# get topsnps from all uncorrelated predictors
rm Intersect/uncorvars_topsnps.range
for var in bio2 bio5 bio6 bio8 bio13 jan_depth snow_days
 do
  echo "${var}"
  cat GenWin/${var}_topsnps.range >> Intersect/uncorvars_topsnps.range
done
```
Problem = more than one topsnp for each candidate window (window is same topsnp is different because from different variable or originally window was split in 2)
Solution = get topsnps from each candidate window with only one topnsp and then get one random topsnp from each window with more than one
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

# get list of candidate windows with ONE topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -u > topsnps_nodupwindows.bed
 
# get topsnp from each window with ONE topsnp
bedtools intersect -a <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) \
 -b topsnps_nodupwindows.bed > topsnps_nodupwindows_onesnp.bed

# get list of candidate windows with MORE than one topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -d > topsnps_dupwindows.bed

# get only one topsnp for each of the windows with duplicates
while read p; do
  bedtools intersect -a <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) -b <(echo "$p") | shuf -n 1
done < topsnps_dupwindows.bed > topsnps_dupwindows_onerand.bed

# get full list of topsnps (only one per window)
cat topsnps_nodupwindows_onesnp.bed topsnps_dupwindows_onerand.bed | sort -k 1,1 -k2,2n \
 > total_intersect_candidate_windows_topsnps.bed
```
Total of 1610 topsnps found from 1610 unique candidate windows - although initially 1620 candidate windows, 11 of them are from the inversion so only 1 SNP will be extracted from them
