# import library
library(Rsubread)
library(here)

# set directory
FILE_DIR <- "/path/to/file_dir"
SAVE_DIR <- "/path/to/save_dir"

ANNOTATION_FILE <- here("reference/gencode.v46.chr_patch_hapl_scaff.annotation.gtf")

# make directory (if not present)
if(!dir.exists(SAVE_DIR)) {
  dir.create(SAVE_DIR, recursive = TRUE)
}

# get BAM file lists
fl <- list.files(FILE_DIR, pattern = "\\.bam$")
print(fl) # check data

# run featurecount
fc <- featureCounts(
  files = paste(FILE_DIR, file, sep="/"),
  isPairedEnd = T,
  annot.ext = ANNOTATION_FILE,
  isGTFAnnotationFile = T,
  nthreads = 4,
  allowMultiOverlap = T
  )

# Define splitLeft function 
removeEnsemblVersion <- function(x) {
  return(strsplit(x, "[.]")[[1]][[1]])
  }

# count data
dt <- fc$counts
rownames(dt) <- sapply(rownames(dt), removeEnsemblVersion) # delete ensembl version

# set sample name as column
colnames(dt) <- gsub(".bam$", "", file)

# save count data
write.table(dt, file = file.path(SAVE_DIR, "counts.txt"), sep="\t", quote=F, col.names=NA)
