#!/bin/sh

#set directories
WORKING_DIR="/path/to/working_dir"
FILE_DIR="/path/to/file_dir"
OUT_DIR="/path/to/output"

# change directory
cd $WORKING_DIR

# list of sample files
FILES_LIST="H6 H15 H16 H17 H20 H21 H30 H38 H39 H42 H43 H44"

# loop through each sample file
for FILE in $FILES_LIST; do
    fastp --detect_adapter_for_pe \
    -i ${FILE_DIR}/${FILE}/${FILE}_1.fq.gz \
    -I ${FILE_DIR}/${FILE}/${FILE}_2.fq.gz \
    -o ${OUT_DIR}/${FILE}_P_1.fq.gz \
    -O ${OUT_DIR}/${FILE}_P_2.fq.gz \
    -h ${OUT_DIR}/${FILE}_report.html \
    -w 4
    done

echo "Fastp trimming and quality reports completed for all samples."