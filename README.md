# 16SrRNA_amplicon_pipeline
R pipeline for processing 16SrRNA amplicon sequencing data including GTDB-based taxonomic assignment using a BLCA approach.
---

## Overview

This workflow processes Illumina paired-end 16S amplicon data and includes the following steps:

1. Quality trimming (Trimmomatic)  
2. Paired-end read merging (USEARCH)  
3. Quality filtering  
4. Primer trimming  
5. Dereplication  
6. Denoising (UNOISE3)  
7. Chimera removal  
8. OTU table generation  
9. Taxonomic assignment using a BLCA-based approach with GTDB reference taxonomy  

---

## Workflow

```
Raw FASTQ
   ↓
Quality trimming (Trimmomatic)
   ↓
Merge paired reads (USEARCH)
   ↓
Quality filtering
   ↓
Primer trimming
   ↓
Dereplication
   ↓
UNOISE denoising
   ↓
Chimera removal
   ↓
OTU table generation
   ↓
Taxonomic assignment (BLCA + GTDB)
```

---

## Primer regions

By default the pipeline targets the **V4 region of the 16S rRNA gene** using the primer pair:

- **515F**
- **806R**

The pipeline can also process **V3–V4 amplicons** using:

- **314F**
- **785R**

To run the pipeline with V3–V4 amplicons, modify the primer patterns in the script.

---

## Taxonomic assignment (BLCA)

Taxonomic classification follows the **Bayesian Lowest Common Ancestor (BLCA)** method described in:

**Gao X., Lin H., Revanna K., Dong Q. (2017)**  
*A Bayesian taxonomic classification method for 16S rRNA gene sequences with improved species-level accuracy.*  
BMC Bioinformatics 18:247.

The BLCA approach is **implemented within this pipeline following the description provided in the publication**.

---

## Reference database

Taxonomy was assigned to the **ASV sequences against a filtered version of the Genome Taxonomy Database (GTDB) version 226**, in which **all sequences shorter than 1000 bp were removed** prior to taxonomic assignment.

Users should note that **new GTDB releases are periodically published**, and the database used in this pipeline should be updated accordingly when newer versions become available.

The GTDB database can be downloaded from:

https://gtdb.ecogenomic.org/

---

## Required software

External tools must be downloaded separately and placed in the `DataNeeded` directory:

- USEARCH  
- Trimmomatic  
- BLAST+  
- Clustal Omega  
- MUSCLE  
- GTDB reference database  

See `DataNeeded/README.md` for details.

---

## Repository structure

```
scripts/
    16S_amplicon_pipeline.R

DataNeeded/
    external tools and databases

example_data/
    optional example FASTQ files
```

---

## Citation

If you use this pipeline, please cite:

Gao X., Lin H., Revanna K., Dong Q. (2017).  
*A Bayesian taxonomic classification method for 16S rRNA gene sequences with improved species-level accuracy.*  
BMC Bioinformatics 18:247.

Additionally, please cite the **Genome Taxonomy Database (GTDB)** when using the reference taxonomy database.
