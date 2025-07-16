#!/bin/sh

#set directories
WORKING_DIR="/path/to/working_dir"
FILE_DIR="/path/to/file_dir"
OUT_DIR="/path/to/output_dir"
BED_FILE="/tRNArRNA.bed"

# move to working directory
cd $WORKING_DIR

# define sample list
FILES_LIST="H6 H15 H16 H17 H20 H21 H30 H38 H39 H42 H43 H44"

# create  if not exists
mkdir -p $OUT_DIR

# Run intersectBed to remove rRNA reads
for FILE in $FILES_LIST; do
    echo "Processing; $FILE"
    intersectBed \
        -abam ${FILE_DIR}/${FILE}.bam \
        -b  $BED_FILE \
        -v > ${OUT_DIR}/${FILE}.bam
done

echo "Processing completed"