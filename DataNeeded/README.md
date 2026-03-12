# External Software Required

This folder contains the external software and reference databases required to run the 16S amplicon pipeline.

These files are not included in the GitHub repository because they may be large or have licensing restrictions.

## Required software

### USEARCH
Used for:
- merging paired reads
- quality filtering
- dereplication
- denoising (UNOISE3)
- chimera removal
- OTU/ASV table generation

Download from:
https://www.drive5.com/usearch/

Example location:
DataNeeded/usearch/usearch11.0.667_i86osx64

---

### Trimmomatic
Used for quality trimming of raw reads.

Download from:
http://www.usadellab.org/cms/?page=trimmomatic

Example location:
DataNeeded/trimmomatic-0.38/trimmomatic-0.38.jar

---

### GTDB reference database

This pipeline was developed using a **filtered version of GTDB release 226** in which **all sequences shorter than 1000 bp were removed** prior to taxonomic assignment.

Required files:

DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_over1000bp.fna  
DataNeeded/GTDB_r226_for_BLCA/GTDB_r226_taxonomy.csv

GTDB can be downloaded from:
https://gtdb.ecogenomic.org/

Users should update these files when new GTDB releases become available.

---

### Clustal Omega

Required for BLCA taxonomy assignment.

Download from:
http://www.clustal.org/omega/

Example location:
DataNeeded/clustalo/

---

### MUSCLE

Required for BLCA taxonomy assignment.

Download from:
https://www.drive5.com/muscle/

Example location:
DataNeeded/muscle/

---

### BLAST+

Required to create the reference database used by BLCA.

Download from:
https://blast.ncbi.nlm.nih.gov/

Make sure the command `makeblastdb` is available in your system path.

---

### Python dependencies

BLCA requires Python 3 and Biopython.

Install Biopython with:

python3 -m pip install biopython --user
