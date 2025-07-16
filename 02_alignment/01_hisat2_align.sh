#!/bin/sh

# set directories
WORKING_DIR="/path/to/working_dir"
GENOME_INDEX="/grch38_snp_tran/genome_snp_tran"
FILE_DIR="/path/to/file_dir"
OUT_DIR="/path/to/output_dir"

# genome_index
# hisat2-build ./genome/GRCh38.p14.genome.fa ./genome/genome_index

# sample list
FILES_LIST="H6 H15 H16 H17 H20 H21 H30 H38 H39 H42 H43 H44"

# move to working directory
cd $WORKING_DIR

# mapping & bam
for FILE in $FILES_LIST; do
    hisat2 \
        --summary-file ${FILE}_hisat2 \
        --new-summary \
        -x $GENOME_INDEX \
        -1 ${FILE_DIR}/${FILE}_P_1.fq.gz \
        -2 ${FILE_DIR}/${FILE}_P_2.fq.gz \
        -k 3 \
        -p 4 \
    | samtools sort -@ 4 -O BAM - > ${OUT_DIR}/${FILE}.bam && \
      samtools index -@ 4 ${OUT_DIR}/${FILE}.bam;\
done

# run multiqc
multiqc $OUT_DIR