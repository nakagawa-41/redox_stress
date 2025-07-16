# redox_stress

# Source Code for Nakagawa et al.

This repository contains code which was used for "Bulk RNA sequencing reveals redox imbalance and impaired protein clearance as key drivers of chronic liver injury in biliary atresia". 

## Quality Check of RNA library. 

**scripts/01_preprocessing/fastp.sh**. 

Check library quality check by using fastp. 

## RNA library: gene mapping and read counts. 

**scripts/02_alignment/01_hisat2_align.sh**. 

Annotation file was obtained from the hisat2 homepage (https://daehwankimlab.github.io/hisat2/download/#h-sapiens). 

H.sapiens GRC38, genome_snp_tran was selected. 
 
 

**scripts/02_alignment/02_intersectBed.sh**. 

To remove remaining rRNA and tRNA, use intersectBed. 

bed.file was obtained from https://genome.ucsc.edu/cgi-bin/hgTables.

Setting (UCSC browser) was "group = variation and repeats, track =repeatmasker, table=rmsk, filter repclass=rRNA"
 
 
**scripts/02_alignment/03_featureCounts.R**. 

For annotation data, the comprehensive gene annotation for Release 46 (GRCh38.p14) was used. 


## Analysis. 

**scripts/03_analysis/visualization.R**. 

DESeq2 was used to identify the DEGs.

source("scripts/03_analysis/visualization.R") returns results. 


**scripts/03_analysis/kegg_pathview.R**. 

KEGG pathview was used to visualize "hsa00190" and "hsa04146".
