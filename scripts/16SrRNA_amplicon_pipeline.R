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

############################################################
# 1. USER-EDITABLE SETTINGS
############################################################

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
merge_minlen   <- 210
merge_maxlen   <- 310

############################
# INPUT / OUTPUT SETTINGS
############################

raw_data_dir <- "01_RawData_B26"
preprocessed_dir <- "02_PreProcessedData"
merged_dir <- "03_MergedData"
filtered_dir <- "04_FilteredData"
trimmed_dir <- "05_TrimmedData"
derep_dir <- "06_DereplicatedData"
joined_dir <- "07_JoinedData"
unique_dir <- "08_UniqueSequences"
denoised_dir <- "09_DenoisedSequences"
uchime_dir <- "10_UchimeReference"
otu_dir <- "11_OtuTable"
tax_dir <- "12_TaxAssignmentGTDB_BLCA_GTDB226"
final_dir <- "13_FinalOtuTableGTDB_BLCA_GTDB226"

############################################################
# 2. SET UP USEARCH
############################################################

# MacOS
system("ln -s 00_DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667")

# Linux / Katana
# system('ln -s 00_DataNeeded/usearch/usearch11.0.667_i86linux64 usearch11.0.667 ; chmod +x usearch11.0.667')

# Test USEARCH
system("./usearch11.0.667")

############################################################
# 3. CREATE OUTPUT DIRECTORIES
############################################################

dir.create(preprocessed_dir, showWarnings = FALSE)
dir.create(merged_dir, showWarnings = FALSE)
dir.create(filtered_dir, showWarnings = FALSE)
dir.create(trimmed_dir, showWarnings = FALSE)
dir.create(derep_dir, showWarnings = FALSE)

############################################################
# 4. PROCESS EACH SAMPLE
############################################################
# Steps per sample:
# - identify R1 / R2 pair
# - trim low-quality regions
# - merge paired reads
# - filter merged reads
# - trim primers
# - dereplicate reads
############################################################

files <- dir(path = raw_data_dir, pattern = "_R1_")
files

for (file1 in files) {

  file2 <- gsub(x = file1, pattern = "_R1_", replacement = "_R2_")
  filename <- gsub(
    x = file1,
    pattern = "_L001_R1_0012.fastq",
    replacement = "",
    perl = TRUE
  )

  ##########################################################
  # 4.1 Quality trimming with Trimmomatic
  ##########################################################
  command <- paste(
    "java -jar 00_DataNeeded/trimmomatic-0.38/trimmomatic-0.38.jar PE ",
    raw_data_dir, "/", file1, " ",
    raw_data_dir, "/", file2, " ",
    preprocessed_dir, "/", filename, "_pF.fastq ",
    preprocessed_dir, "/", filename, "_upF.fastq ",
    preprocessed_dir, "/", filename, "_pR.fastq ",
    preprocessed_dir, "/", filename, "_upR.fastq ",
    "HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:100",
    sep = ""
  )
  system(command)

  ##########################################################
  # 4.2 Merge paired reads with USEARCH
  ##########################################################
  command <- paste(
    "./usearch11.0.667 -fastq_mergepairs ", preprocessed_dir, "/", filename, "_pF.fastq ",
    "-reverse ", preprocessed_dir, "/", filename, "_pR.fastq ",
    "-fastqout ", merged_dir, "/", filename, ".fastq ",
    "-relabel @ ",
    "-minhsp 12 ",
    "-fastq_minovlen 8 ",
    "-fastq_maxdiffs ", merge_maxdiffs, " ",
    "-fastq_pctid ", merge_pctid, " ",
    "-fastq_minmergelen ", merge_minlen, " ",
    "-fastq_maxmergelen ", merge_maxlen, " ",
    "-sample ", filename,
    sep = ""
  )
  system(command)

  ##########################################################
  # 4.3 Quality filtering of merged reads
  ##########################################################
  command <- paste(
    "./usearch11.0.667 -fastq_filter ", merged_dir, "/", filename, ".fastq ",
    "-fastaout ", filtered_dir, "/", filename, ".fasta ",
    "-fastq_maxns 1 -fastq_maxee 1",
    sep = ""
  )
  system(command)

  ##########################################################
  # 4.4 Primer trimming
  ##########################################################
  # Reads are checked for:
  # 1. forward primer at the start of the sequence
  # 2. reverse primer at the end of the sequence
  #
  # If a primer is found, it is removed.
  # If not found, the original sequence is kept.
  ##########################################################

  FQ <- read.fasta(
    file = paste0(filtered_dir, "/", filename, ".fasta"),
    as.string = TRUE
  )

  IDs <- names(FQ)
  SEQs <- as.character(FQ)

  # Check/remove forward primer
  for_primer_present <- grepl(
    pattern = forward_primer,
    x = SEQs,
    ignore.case = TRUE
  )

  SEQs <- ifelse(
    for_primer_present,
    gsub(
      pattern = forward_primer,
      replacement = "",
      SEQs,
      ignore.case = TRUE
    ),
    SEQs
  )

  # Check/remove reverse primer
  rev_primer_present <- grepl(
    pattern = reverse_primer,
    x = SEQs,
    ignore.case = TRUE
  )

  SEQs <- ifelse(
    rev_primer_present,
    gsub(
      pattern = reverse_primer,
      replacement = "",
      SEQs,
      ignore.case = TRUE
    ),
    SEQs
  )

  OUT <- file(
    description = paste0(trimmed_dir, "/", filename, ".fasta"),
    open = "w"
  )

  for (i in 1:length(IDs)) {
    write(x = paste0(">", IDs[i], "\n", SEQs[i]), file = OUT)
  }
  close(OUT)

  ##########################################################
  # 4.5 Dereplication per sample
  ##########################################################
  command <- paste(
    "./usearch11.0.667 -fastx_uniques ",
    trimmed_dir, "/", filename, ".fasta ",
    "-fastaout ", derep_dir, "/", filename, ".fasta ",
    "-sizeout",
    sep = ""
  )
  system(command)
}

############################################################
# 5. COMPARE FILTERED VS TRIMMED READS
############################################################
# This section checks:
# - whether read counts stayed the same after primer trimming
# - which reads changed
# - whether trimming affected the forward primer, reverse primer, both, or neither
#
# Output files:
# - trimmed_vs_original_all_reads.csv
# - trimmed_vs_original_summary_by_file.csv
# - trimmed_vs_original_changed_reads_only.csv
############################################################

files <- dir(filtered_dir, pattern = "\\.fasta$")
all_results <- list()

for (f in files) {

  orig <- read.fasta(
    file = file.path(filtered_dir, f),
    as.string = TRUE
  )

  trim <- read.fasta(
    file = file.path(trimmed_dir, f),
    as.string = TRUE
  )

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

      seq_after_forward = ifelse(
        !is.na(seq_orig),
        ifelse(
          grepl(forward_primer, seq_orig, ignore.case = TRUE),
          gsub(forward_primer, "", seq_orig, ignore.case = TRUE),
          seq_orig
        ),
        NA
      ),

      reverse_present_after_forward = ifelse(
        !is.na(seq_after_forward),
        grepl(reverse_primer, seq_after_forward, ignore.case = TRUE),
        NA
      ),

      changed = ifelse(!is.na(seq_orig) & !is.na(seq_trim), seq_orig != seq_trim, NA),

      trim_type = case_when(
        !present_in_filtered ~ "missing_in_filtered",
        !present_in_trimmed ~ "missing_in_trimmed",
        changed & forward_present_before & reverse_present_after_forward ~ "both",
        changed & forward_present_before & !reverse_present_after_forward ~ "forward_only",
        changed & !forward_present_before & reverse_present_after_forward ~ "reverse_only",
        !changed ~ "none",
        changed ~ "changed_other",
        TRUE ~ "unknown"
      )
    )

  all_results[[f]] <- df_compare
}

all_results_df <- bind_rows(all_results)

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

summary_by_file$max_bases_removed[is.infinite(summary_by_file$max_bases_removed)] <- NA
summary_by_file$min_bases_removed[is.infinite(summary_by_file$min_bases_removed)] <- NA

changed_reads <- all_results_df %>%
  filter(changed == TRUE)

write.csv(
  all_results_df,
  "trimmed_vs_original_all_reads.csv",
  row.names = FALSE
)

write.csv(
  summary_by_file,
  "trimmed_vs_original_summary_by_file.csv",
  row.names = FALSE
)

write.csv(
  changed_reads,
  "trimmed_vs_original_changed_reads_only.csv",
  row.names = FALSE
)

cat("\n===== OVERALL SUMMARY =====\n")
print(table(all_results_df$trim_type, useNA = "ifany"))

cat("\n===== PER-FILE SUMMARY =====\n")
print(summary_by_file)

cat("\n===== FIRST CHANGED READS =====\n")
print(
  changed_reads %>%
    select(
      file, ID, len_orig, len_trim, bases_removed,
      forward_present_before, reverse_present_after_forward, trim_type
    ) %>%
    head(20)
)

cat("\nFiles written:\n")
cat(" - trimmed_vs_original_all_reads.csv\n")
cat(" - trimmed_vs_original_summary_by_file.csv\n")
cat(" - trimmed_vs_original_changed_reads_only.csv\n")

############################################################
# 6. CHECK MERGED/TRIMMED SEQUENCE LENGTH DISTRIBUTION
############################################################
# This section helps assess whether the selected merge-length
# settings are appropriate for the region under study.
############################################################

seq_lengths <- NULL

for (f in dir(path = trimmed_dir, full.names = TRUE)) {
  print(f)
  fasta_file <- read.fasta(f)
  seq_lengths <- c(seq_lengths, as.numeric(lengths(fasta_file)))
}

hist(seq_lengths)
quantile(seq_lengths)

pdf("hist_post_adjustment_ALLsamples.pdf")
hist(seq_lengths)
dev.off()

# Example interpretation
sum(seq_lengths < 200) / sum(seq_lengths) * 100
sum(seq_lengths > 500) / sum(seq_lengths) * 100

############################################################
# 7. JOIN ALL DEREPLICATED FASTA FILES
############################################################

dir.create(joined_dir, showWarnings = FALSE)

system(paste(
  "cat", derep_dir, "/*.fasta >", joined_dir, "/AllSamples.fasta"
)) # Linux and MacOS

# shell("type 06_DereplicatedData\\*.fasta > 07_JoinedData\\AllSamples.fasta") # Windows

FNA <- readLines(paste0(joined_dir, "/AllSamples.fasta"))
FNA[grep(pattern = ">", x = FNA, invert = TRUE)] <- toupper(
  FNA[grep(pattern = ">", x = FNA, invert = TRUE)]
)
write(x = FNA, file = paste0(joined_dir, "/AllSamples2.fasta"))

############################################################
# 8. DEREPLICATE ACROSS ALL SAMPLES
############################################################

dir.create(unique_dir, showWarnings = FALSE)

system(paste(
  "./usearch11.0.667 -fastx_uniques ",
  joined_dir, "/AllSamples2.fasta ",
  "-fastaout ", unique_dir, "/AllSamples_uniques.fasta ",
  "-sizein -sizeout -strand both",
  sep = ""
))

############################################################
# 9. DENOISE WITH UNOISE
############################################################

dir.create(denoised_dir, showWarnings = FALSE)

system(paste(
  "./usearch11.0.667 -unoise3 ",
  unique_dir, "/AllSamples_uniques.fasta ",
  "-zotus ", denoised_dir, "/AllSamples_denoised.fasta",
  sep = ""
))

############################################################
# 10. CHIMERA REMOVAL
############################################################

dir.create(uchime_dir, showWarnings = FALSE)

system(paste(
  "./usearch11.0.667 -uchime2_ref ",
  denoised_dir, "/AllSamples_denoised.fasta ",
  "-db 00_DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna ",
  "-strand plus -mode high_confidence ",
  "-notmatched ", uchime_dir, "/AllSamples_unoise_nc.fasta",
  sep = ""
))

############################################################
# 11. GENERATE OTU TABLE
############################################################

dir.create(otu_dir, showWarnings = FALSE)

system(paste(
  "./usearch11.0.667 -otutab ",
  joined_dir, "/AllSamples2.fasta ",
  "-otus ", uchime_dir, "/AllSamples_unoise_nc.fasta ",
  "-id 0.97 ",
  "-otutabout ", otu_dir, "/AllSamples_unoise_otu_table1.txt",
  sep = ""
))

############################################################
# 12. PREPARE BLCA DEPENDENCIES
############################################################
# BLCA requires clustalo, muscle, BLAST, Python, and Biopython.
############################################################

# MacOS
system("cp 00_DataNeeded/clustalo/clustal-omega-1.2.3-macosx 00_DataNeeded/clustalo/clustalo ; chmod +x 00_DataNeeded/clustalo/clustalo")
system("cp 00_DataNeeded/muscle/muscle3.8.31_i86darwin64 00_DataNeeded/muscle/muscle ; chmod +x 00_DataNeeded/muscle/muscle")

# Windows notes:
# Run manually in terminal if needed:
# copy 00_DataNeeded\\clustalo\\clustal-omega-1.2.2-win64 00_DataNeeded\\clustalo\\clustalo
# attrib +x 00_DataNeeded\\clustalo\\clustalo
# copy 00_DataNeeded\\muscle\\muscle3.8.31_i86darwin64 00_DataNeeded\\muscle\\muscle
# attrib +x 00_DataNeeded\\muscle\\muscle

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/clustalo", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/muscle", sep = ":"))

############################################################
# 13. RUN BLCA AGAINST GTDB
############################################################

dir.create(tax_dir, showWarnings = FALSE)

# If needed, install Biopython
system(paste(system("which python3", intern = TRUE), "-m pip install biopython --user"))

# Windows alternative:
# py -m pip install biopython

system("git clone https://github.com/qunfengdong/BLCA.git")
system("makeblastdb -in 00_DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna -parse_seqids -blastdb_version 5 -dbtype nucl")

system(paste(
  "python3 BLCA/2.blca_main.py ",
  "-r 00_DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_taxonomy_oldF.csv ",
  "-q 00_DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna ",
  "-i ", uchime_dir, "/AllSamples_unoise_nc.fasta ",
  "-o ", tax_dir, "/AllSamples_unoise_BLCA_GTDB226_out.1.txt",
  sep = ""
))

############################################################
# 14. CLEAN BLCA TAXONOMY OUTPUT
############################################################
# This section removes taxonomic rank prefixes from BLCA output
# and keeps confidence values in parentheses for readability.
############################################################

m <- 1

for (each_line in readLines(file(paste0(tax_dir, "/AllSamples_unoise_BLCA_GTDB226_out.1.txt"), open = "r"))) {

  each_line_split <- strsplit(each_line, "\t")
  OTU_ID <- each_line_split[[1]][1]
  taxonomy <- each_line_split[[1]][2]
  taxonomy_split <- strsplit(taxonomy, ";")
  taxonomy_no_rank <- ""
  n <- 1

  for (taxon in taxonomy_split[[1]]) {
    if (n %% 2 == 1) {
      taxon_split <- strsplit(taxon, ":")
      if (length(taxon_split[[1]]) == 2) {
        taxon_no_rank <- taxon_split[[1]][2]
      } else {
        taxon_no_rank <- taxon_split[[1]][1]
      }
      taxonomy_no_rank <- paste(taxonomy_no_rank, taxon_no_rank, sep = ");")
    } else {
      taxonomy_no_rank <- paste(taxonomy_no_rank, taxon, sep = "(")
    }
    n <- n + 1
  }

  taxonomy_no_rank <- paste(taxonomy_no_rank, ")", sep = "")
  taxonomy_no_rank <- substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))

  if (taxonomy_no_rank == "Unclassified)") {
    taxonomy_no_rank <- "Unclassified"
  }

  taxonomy_no_rank_with_OTU <- paste(OTU_ID, taxonomy_no_rank, sep = "\t")

  if (m == 1) {
    cat(
      taxonomy_no_rank_with_OTU,
      file = paste0(tax_dir, "/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt"),
      sep = "\n",
      append = FALSE
    )
  } else {
    cat(
      taxonomy_no_rank_with_OTU,
      file = paste0(tax_dir, "/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt"),
      sep = "\n",
      append = TRUE
    )
  }

  m <- m + 1
}

############################################################
# 15. MERGE OTU TABLE AND TAXONOMY
############################################################

dir.create(final_dir, showWarnings = FALSE)

OTU <- read.delim(
  paste0(otu_dir, "/AllSamples_unoise_otu_table1.txt"),
  header = TRUE
)

# Note:
# The taxonomy output may still contain confidence values in parentheses
# (e.g. 100.0000). If you manually edit/filter that file, document it clearly.
TAX <- read.delim(
  paste0(tax_dir, "/AllSamples_unoise_nc.fasta.blca_GTDB226.2.txt"),
  header = FALSE
)

names(TAX) <- c("X.OTU.ID", "Taxonomy")

OTU_TAX <- merge(OTU, TAX, by = "X.OTU.ID")

write.table(
  OTU_TAX,
  paste0(final_dir, "/AllSamples_unoise_otu_table_BLCA_filtered2.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# rm(list = ls(all = TRUE))
