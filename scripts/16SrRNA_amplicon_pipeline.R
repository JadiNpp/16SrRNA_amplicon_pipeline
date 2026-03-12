############################################################
# 16S rRNA Amplicon Processing Pipeline
# Author: Jadranka Nappi
# Year: 2026
#
# Description:
# This script processes paired-end 16S amplicon data using:
# - Trimmomatic for quality trimming
# - USEARCH for merging, filtering, dereplication, denoising,
#   chimera removal, and OTU/ASV table generation
# - GTDB r226 as reference database
# - BLCA-based taxonomic assignment
#
# Notes:
# - Taxonomy is assigned against a filtered version of GTDB
#   release 226, in which all sequences shorter than 1000 bp
#   were removed before taxonomic assignment.
# - Update the GTDB files when newer database releases become
#   available.
# - Folder names are kept similar to the original working script.
############################################################

############################
# LOAD LIBRARIES
############################

library(tidyverse)
library(reshape2)
library(seqinr)
library(dplyr)

############################
# PREPARE USEARCH
############################

system("ln -s DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667") # MacOS
#system("ln -s DataNeeded/usearch/usearch11.0.667_i86linux64 usearch11.0.667 ; chmod +x usearch11.0.667") # Linux / Katana
system('./usearch11.0.667') # Test

############################
# PRIMER SETTINGS
############################

# Default primers: V4 region of the 16S rRNA gene
# Primer pair: 515F / 806R
forward_primer <- "^[ATGC]{0,2}GTG[CT]CAGC[AC]GCCGCGGTAA"
reverse_primer <- "ATTAGA[AT]ACCC{0,2}$"

# If using the V3–V4 region instead (314F / 785R),
# replace the primers above with:
# forward_primer <- "^[ATGC]{0,2}CCTACGGG[ATGC]GGC[AT]GCAG"
# reverse_primer <- "GGATTAGATACCC[CGT][AGT]GTAGTC[ATGC]{0,2}$"

############################
# MERGE SETTINGS
############################

# These settings should be modified according to the
# 16S region under study.
#
# For example:
# - V4 (515F / 806R) generally produces shorter amplicons
# - V3–V4 (314F / 785R) generally produces longer amplicons
#
# Adjust the minimum and maximum merged length accordingly.

merge_maxdiffs <- 5
merge_pctid    <- 80
merge_minlen   <- 200
merge_maxlen   <- 500

############################
# CREATE OUTPUT FOLDERS
############################

dir.create("02_PreProcessedData", showWarnings = FALSE)
dir.create("03_MergedData", showWarnings = FALSE)
dir.create("04_FilteredData", showWarnings = FALSE)
dir.create("05_TrimmedData", showWarnings = FALSE)
dir.create("06_DereplicatedData", showWarnings = FALSE)
dir.create("07_JoinedData", showWarnings = FALSE)
dir.create("08_UniqueSequences", showWarnings = FALSE)
dir.create("09_DenoisedSequences", showWarnings = FALSE)
dir.create("10_UchimeReference", showWarnings = FALSE)
dir.create("11_OtuTable", showWarnings = FALSE)
dir.create("12_TaxAssignmentGTDB_BLCA_GTDB226", showWarnings = FALSE)
dir.create("13_FinalOtuTableGTDB_BLCA_GTDB226", showWarnings = FALSE)

############################
# STEP 1. PREPROCESS EACH SAMPLE
############################

files = dir(path = '01_RawData/', pattern = '_R1_')
files

for (file1 in files)
{
  file2 = gsub(x = file1, pattern = '_R1_', replacement = '_R2_')

  # General sample name extraction.
  # Adjust this if your file naming convention differs.
  filename = gsub(x = file1, pattern = "_R1_.*$", replacement = '', perl = TRUE)

  # Quality trimming
  command <- paste(
    'java -jar DataNeeded/trimmomatic-0.38/trimmomatic-0.38.jar PE 01_RawData/', file1,
    ' 01_RawData/', file2,
    ' 02_PreProcessedData/', filename, '_pF.fastq',
    ' 02_PreProcessedData/', filename, '_upF.fastq',
    ' 02_PreProcessedData/', filename, '_pR.fastq',
    ' 02_PreProcessedData/', filename, '_upR.fastq',
    ' HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:100',
    sep = ""
  )
  system(command)

  # Merge paired reads
  command = paste(
    "./usearch11.0.667 -fastq_mergepairs 02_PreProcessedData/", filename,
    "_pF.fastq -reverse 02_PreProcessedData/", filename,
    "_pR.fastq -fastqout 03_MergedData/", filename,
    ".fastq -relabel @ -fastq_maxdiffs ", merge_maxdiffs,
    " -fastq_pctid ", merge_pctid,
    " -fastq_minmergelen ", merge_minlen,
    " -fastq_maxmergelen ", merge_maxlen,
    " -sample ", filename,
    sep = ""
  )
  system(command)

  # Quality filtering
  command = paste(
    "./usearch11.0.667 -fastq_filter 03_MergedData/", filename,
    ".fastq -fastaout 04_FilteredData/", filename,
    ".fasta -fastq_maxns 1 -fastq_maxee 1",
    sep = ""
  )
  system(command)

  # Check and remove primers
  FQ = read.fasta(file = paste0("04_FilteredData/", filename, ".fasta"), as.string = TRUE)

  IDs = names(FQ)
  SEQs = as.character(FQ)

  SEQs = gsub(pattern = forward_primer, replacement = "", SEQs, ignore.case = TRUE)
  SEQs = gsub(pattern = reverse_primer, replacement = "", SEQs, ignore.case = TRUE)

  OUT = file(description = paste0("05_TrimmedData/", filename, ".fasta"), open = "w")
  for(i in 1:length(IDs)) write(x = paste0(">", IDs[i], "\n", SEQs[i]), file = OUT)
  close(OUT)

  # Dereplication
  command = paste(
    "./usearch11.0.667 -fastx_uniques 05_TrimmedData/", filename,
    ".fasta -fastaout 06_DereplicatedData/", filename,
    ".fasta -sizeout",
    sep = ""
  )
  system(command)
}

############################
# STEP 2. COMPARE FILTERED VS TRIMMED READS
############################

# This section compares 04_FilteredData vs 05_TrimmedData
# to confirm whether the number of sequences stayed the same
# and to show where trimming happened.

files <- dir("04_FilteredData", pattern = "\\.fasta$")

all_results <- list()

for (f in files) {

  # Read original filtered file
  orig <- read.fasta(
    file = file.path("04_FilteredData", f),
    as.string = TRUE
  )

  # Read trimmed file
  trim <- read.fasta(
    file = file.path("05_TrimmedData", f),
    as.string = TRUE
  )

  # Convert to data frames
  df_orig <- data.frame(
    file = f,
    ID = names(orig),
    seq_orig = as.character(orig),
    stringsAsFactors = FALSE
  )

  df_trim <- data.frame(
    file = f,
    ID = names(trim),
    seq_trim = as.character(trim),
    stringsAsFactors = FALSE
  )

  # Compare
  df_compare <- full_join(df_orig, df_trim, by = c("file", "ID")) %>%
    mutate(
      present_in_filtered = !is.na(seq_orig),
      present_in_trimmed = !is.na(seq_trim),
      len_orig = ifelse(!is.na(seq_orig), nchar(seq_orig), NA),
      len_trim = ifelse(!is.na(seq_trim), nchar(seq_trim), NA),
      bases_removed = ifelse(!is.na(len_orig) & !is.na(len_trim), len_orig - len_trim, NA),

      forward_present_before = ifelse(
        !is.na(seq_orig),
        grepl(forward_primer, seq_orig, ignore.case = TRUE),
        NA
      ),

      reverse_present_before = ifelse(
        !is.na(seq_orig),
        grepl(reverse_primer, seq_orig, ignore.case = TRUE),
        NA
      ),

      changed = ifelse(!is.na(seq_orig) & !is.na(seq_trim), seq_orig != seq_trim, NA),

      trim_type = case_when(
        !present_in_filtered ~ "missing_in_filtered",
        !present_in_trimmed ~ "missing_in_trimmed",
        changed & forward_present_before & reverse_present_before ~ "both",
        changed & forward_present_before & !reverse_present_before ~ "forward_only",
        changed & !forward_present_before & reverse_present_before ~ "reverse_only",
        !changed ~ "none",
        changed ~ "changed_other",
        TRUE ~ "unknown"
      )
    )

  all_results[[f]] <- df_compare
}

# Combine all files
all_results_df <- bind_rows(all_results)

# Per-sample summary
summary_by_file <- all_results_df %>%
  group_by(file) %>%
  summarise(
    n_reads_filtered = sum(present_in_filtered, na.rm = TRUE),
    n_reads_trimmed = sum(present_in_trimmed, na.rm = TRUE),
    same_read_count = n_reads_filtered == n_reads_trimmed,
    n_changed = sum(changed, na.rm = TRUE),
    n_unchanged = sum(trim_type == "none", na.rm = TRUE),
    n_forward_only = sum(trim_type == "forward_only", na.rm = TRUE),
    n_reverse_only = sum(trim_type == "reverse_only", na.rm = TRUE),
    n_both = sum(trim_type == "both", na.rm = TRUE),
    n_changed_other = sum(trim_type == "changed_other", na.rm = TRUE),
    n_missing_in_trimmed = sum(trim_type == "missing_in_trimmed", na.rm = TRUE),
    max_bases_removed = suppressWarnings(max(bases_removed, na.rm = TRUE)),
    min_bases_removed = suppressWarnings(min(bases_removed, na.rm = TRUE)),
    .groups = "drop"
  )

# Fix Inf / -Inf if no valid values
summary_by_file$max_bases_removed[is.infinite(summary_by_file$max_bases_removed)] <- NA
summary_by_file$min_bases_removed[is.infinite(summary_by_file$min_bases_removed)] <- NA

# Only reads where something changed
changed_reads <- all_results_df %>%
  filter(changed == TRUE)

# Save outputs
write.csv(all_results_df,
          "trimmed_vs_original_all_reads.csv",
          row.names = FALSE)

write.csv(summary_by_file,
          "trimmed_vs_original_summary_by_file.csv",
          row.names = FALSE)

write.csv(changed_reads,
          "trimmed_vs_original_changed_reads_only.csv",
          row.names = FALSE)

# Print useful summaries
cat("\n===== OVERALL SUMMARY =====\n")
print(table(all_results_df$trim_type, useNA = "ifany"))

cat("\n===== PER-FILE SUMMARY =====\n")
print(summary_by_file)

cat("\n===== FIRST CHANGED READS =====\n")
print(
  changed_reads %>%
    select(file, ID, len_orig, len_trim, bases_removed,
           forward_present_before, reverse_present_before, trim_type) %>%
    head(20)
)

cat("\nFiles written:\n")
cat(" - trimmed_vs_original_all_reads.csv\n")
cat(" - trimmed_vs_original_summary_by_file.csv\n")
cat(" - trimmed_vs_original_changed_reads_only.csv\n")

############################
# STEP 3. CHECK MERGED RANGE / TRIMMED LENGTHS
############################

seq_lengths = NULL
for(f in dir(path = '05_TrimmedData/', full.names = TRUE))
{
  print(f)
  fasta_file = read.fasta(f)
  seq_lengths = c(seq_lengths, as.numeric(lengths(fasta_file)))
}

hist(seq_lengths) # Use this to adjust merge length range if needed
quantile(seq_lengths)

pdf("hist_post_adjustment_all_samples.pdf")
hist(seq_lengths)
dev.off()

sum(seq_lengths < merge_minlen) / sum(seq_lengths) * 100
sum(seq_lengths > merge_maxlen) / sum(seq_lengths) * 100

############################
# STEP 4. JOIN ALL DEREPLICATED FILES
############################

system("cat 06_DereplicatedData/*.fasta > 07_JoinedData/AllSamples.fasta") # Linux and MacOS
#shell("type 06_DereplicatedData\\*.fasta > 07_JoinedData\\AllSamples.fasta") # Windows

FNA = readLines("07_JoinedData/AllSamples.fasta")
FNA[grep(pattern = ">", x = FNA, invert = TRUE)] = toupper(FNA[grep(pattern = ">", x = FNA, invert = TRUE)])
write(x = FNA, file = "07_JoinedData/AllSamples2.fasta")

############################
# STEP 5. DEREPLICATE ACROSS ALL SAMPLES
############################

system("./usearch11.0.667 -fastx_uniques 07_JoinedData/AllSamples2.fasta -fastaout 08_UniqueSequences/AllSamples_uniques.fasta -sizein -sizeout -strand both")

############################
# STEP 6. GENERATE DENOISED SEQUENCES WITH UNOISE3
############################

system("./usearch11.0.667 -unoise3 08_UniqueSequences/AllSamples_uniques.fasta -zotus 09_DenoisedSequences/AllSamples_denoised.fasta")

############################
# STEP 7. CHIMERA REMOVAL
############################

system("./usearch11.0.667 -uchime2_ref 09_DenoisedSequences/AllSamples_denoised.fasta -db DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna -strand plus -mode high_confidence -notmatched 10_UchimeReference/AllSamples_unoise_nc.fasta")

############################
# STEP 8. OTU / ASV TABLE GENERATION
############################

system("./usearch11.0.667 -otutab 07_JoinedData/AllSamples2.fasta -otus 10_UchimeReference/AllSamples_unoise_nc.fasta -id 0.97 -otutabout 11_OtuTable/AllSamples_unoise_otu_table1.txt")

############################
# STEP 9. BLCA PREPARATION
############################

# Add clustalo and muscle to system path (required by BLCA)

# FOR MacOS:
system("cp DataNeeded/clustalo/clustal-omega-1.2.3-macosx DataNeeded/clustalo/clustalo ; chmod +x DataNeeded/clustalo/clustalo")
system("cp DataNeeded/muscle/muscle3.8.31_i86darwin64 DataNeeded/muscle/muscle ; chmod +x DataNeeded/muscle/muscle")

# FOR WINDOWS: run these in terminal if needed
# copy DataNeeded\\clustalo\\clustal-omega-1.2.2-win64 DataNeeded\\clustalo\\clustalo
# attrib +x DataNeeded\\clustalo\\clustalo
# copy DataNeeded\\muscle\\muscle3.8.31_i86darwin64 DataNeeded\\muscle\\muscle
# attrib +x DataNeeded\\muscle\\muscle

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "DataNeeded/clustalo", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "DataNeeded/muscle", sep = ":"))

# Install Biopython if needed
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

############################
# STEP 10. RUN BLCA
############################

# If you are using another local implementation of BLCA,
# replace this section accordingly.

system("git clone https://github.com/qunfengdong/BLCA.git")
system("makeblastdb -in DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna -parse_seqids -blastdb_version 5 -dbtype nucl")

system('python3 BLCA/2.blca_main.py -r DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_taxonomy.csv -q DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna -i 10_UchimeReference/AllSamples_unoise_nc.fasta -o 12_TaxAssignmentGTDB_BLCA_GTDB226/AllSamples_unoise_BLCA_GTDB226_out.1.txt')

############################
# STEP 11. CLEAN BLCA OUTPUT
############################

m = 1
for (each_line in readLines(file('12_TaxAssignmentGTDB_BLCA_GTDB226/AllSamples_unoise_BLCA_GTDB226_out.1.txt', open = "r")) ){
  each_line_split = strsplit(each_line, '\t')
  OTU_ID = each_line_split[[1]][1]
  taxonomy = each_line_split[[1]][2]
  taxonomy_split = strsplit(taxonomy, ';')
  taxonomy_no_rank = ''
  n = 1

  for (taxon in taxonomy_split[[1]]){
    if (n %% 2 == 1){
      taxon_split = strsplit(taxon, ':')
      if (length(taxon_split[[1]]) == 2)
      {taxon_no_rank = taxon_split[[1]][2]}
      else
      {taxon_no_rank = taxon_split[[1]][1]}
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon_no_rank, sep = ");")
    } else {
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon, sep = "(")
    }
    n = n + 1
  }

  taxonomy_no_rank = paste(taxonomy_no_rank, ')', sep = "")
  taxonomy_no_rank = substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))
  if (taxonomy_no_rank == "Unclassified)"){taxonomy_no_rank = "Unclassified"}

  taxonomy_no_rank_with_OTU = paste(OTU_ID, taxonomy_no_rank, sep = "\t")

  if (m == 1)
  {cat(taxonomy_no_rank_with_OTU, file = '12_TaxAssignmentGTDB_BLCA_GTDB226/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt', sep = "\n", append = FALSE)}
  else
  {cat(taxonomy_no_rank_with_OTU, file = '12_TaxAssignmentGTDB_BLCA_GTDB226/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt', sep = "\n", append = TRUE)}
  m = m + 1
}

############################
# STEP 12. MERGE OTU TABLE AND TAXONOMY
############################

OTU = read.delim("11_OtuTable/AllSamples_unoise_otu_table1.txt", header = TRUE)
TAX = read.delim("12_TaxAssignmentGTDB_BLCA_GTDB226/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt", header = FALSE)
names(TAX) = c("X.OTU.ID", "Taxonomy")

OTU_TAX = merge(OTU, TAX, by = "X.OTU.ID")

write.table(
  OTU_TAX,
  "13_FinalOtuTableGTDB_BLCA_GTDB226/AllSamples_unoise_otu_table_BLCA.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("\nPipeline completed successfully.\n")
