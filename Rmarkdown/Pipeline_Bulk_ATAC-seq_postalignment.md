Bulk ATAC-seq Postalignment Pipeline
================
Victor Enrique Goitea
2024-03-11

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Description](#description)
- [Output](#output)
- [Output files reorganization
  (optional)](#output-files-reorganization-optional)
  - [Folder Structure:](#folder-structure)
  - [Creating a Multiqc report:](#creating-a-multiqc-report)

# Overview

This script is designed to process BAM files generated from sequencing
data. It performs various quality control (QC) analyses and data
preprocessing steps to prepare the data for downstream analysis. The
script utilizes several bioinformatics tools and workflows.

# Prerequisites

- SLURM scheduler.  
- Modules: SAMtools, Picard tools, Bedtools, Qualimap, Preseq,
  Phantompeakqualtools, Homer, Deeptools, MACS2.  
- Reference genomes (provided in /work/References).

# Usage

1.  Modify SLURM Parameters (Optional): Open the script
    (**pe_postalign_ATACseq_multigenome.sh**) and modify SLURM
    parameters at the beginning of the file, such as account, output
    file, email notifications, nodes, memory, CPU cores, and runtime.
    Alternatively, you can modify these parameters on-the-fly when
    executing the script.

2.  On UCloud, start a **Terminal Ubuntu** run:

    - Enable **Slurm cluster**
    - To process several samples consider requesting nodes \> 1
    - Set the modules path to **FGM \> Utilities \> App \> easybuild**

    ![](../Img/terminal_slurm.png)

    - Include the References folder **FGM \> References \> References**

    ![](../Img/terminal_folders.png)

    - Include your Scripts folder and the folder with the bam files.

    - **Notes:**

      - Make sure the scripts have executing permission. If not run:
        `chmod 700 script.sh`
      - Match the job CPUs to the amounts requested in the script.
      - If you modify the memory parameter in the script, specify 5-10%
        less than the memory available in the terminal run.
      - Although it is not necessary to enable **tmux**, it is a good
        practise to always do it.
      - The script also uses some files in `References/Multiqc/` to make
        some metric files compatible with Multiqc.

3.  **Run the Script:** Submit the script to the SLURM cluster:

        sbatch -J <job_name> path_to/Scripts_folder/pe_postalign_ATACseq_multigenome.sh -g <mm10|mm39|hg38> <input-bam-file> 

    **Required Arguments**

    - **-g:** specify the genome to use (mm10, mm39, or hg38).
    - Replace **input-bam-file** with the path to your bam file.

    For several samples you can use a for loop:

        for i in *.bam; do sbatch -J <job_name> path_to/Scripts_folder/pe_postalign_ATACseq_multigenome.sh -g <mm10|mm39|hg38> $i; sleep 1; done

4.  **Monitor Job:** You can monitor the job using the SLURM commands,
    such as **squeue**, **scontrol show job <job-id>**, and check the
    log files generated.

# Description

This script performs the following main tasks:

1.  **Extensive quality control:**
    - **Samtools:** stats, idxstats, flagstat.
    - **Picard:** CollectAlignmentSummaryMetrics,
      QualityScoreDistribution, MeanQualityByCycle,
      CollectInsertSizeMetrics, EstimateLibraryComplexity,
      Markduplicates (flag/not remove).
    - **Qualimap:** bamqc.
    - **Library Complexity Estimation:** Preseq, Picard and PBC metrics.
    - **Deeptools:** plotCoverage, estimateReadFiltering,
      plotFingerprint (in 10 million mates subsample)
    - **Phantompeakqualtools:** cross-correlation, NSC and RSC
      coefficients.
    - **Agregate profile around genes:** genes were scaled to a region
      body length of 1 kb and -3/+3 kb flanking regiions were included
      (Deeptools).
2.  **Fixing Coodinates:** for example if chromosome nomenclature is
    “1”, convert to “chr1”. It keeps only canonical chromosomes
    including sex and mitochondrial chromosomes.  
3.  **Filtering and Sorting:** filters and sorts the processed BAM files
    to remove low-quality reads and ensure proper alignment:
    - Files with duplicates marked but not removed, i.e. **samtools -q
      30 -F 780**
    - Files with duplicates removed and only properly paired reads are
      kept, i.e. **samtools -q 30 -F 1804 -f 2**
4.  **Reads splitting:** Bam files are split into NFR bin (38-100 bp)
    and MN bin (180-247 bp) using **alignmentSieve** (Deeptools) with
    **–ATACshift** to produce the shift commonly done for ATAC-seq,
    i.e. equivalent to –shift 4 -5 5 -4.
5.  **Indexing:** indexes the final BAM files.
6.  **Bigwig coverage tracks:** generates coverage tracks with bin size
    1 and 10, the former is RPKM normalized by **bamCoverrage**
    (Deeptools). Blacklist has been taken into consideration to generate
    the tracks.
7.  **Tag Directories:** generates a sample tag directory using Homer
    for downstream analysis.
8.  **PEAK call:** macs2 callpeak in NFR bam files and also calculates
    FRiP scores.
9.  **GENCODE GTF Version Information:**
    - **Human (GRCh38):** v44
    - **Mouse (GRCm38 - mm10):** M25
    - **Mouse (GRCm39 - mm39):** M33
10. **Blacklists:** version 2 of Encode blacklists. Blacklist regions
    separated less than 1000 bp were merged.

# Output

The script generates various output files containing QC metrics,
alignment statistics, and processed BAM files ready for downstream
analysis. The script store the output files in a directory of basename
<input-filename> (Main output directory).

# Output files reorganization (optional)

After running the script for all samples, each sample will have its own
folder with the basename <input-filename> (Main output directory). If
you prefer to organize your files by category, you can execute the
provided script **reorganize_files_atacseq.sh** in the terminal.

## Folder Structure:

    ├── FASTQ_Raw
    │   ├── *fastq.gz
    ├── BAM_Raw
    │   ├── *bam
    │   ├── *bai
    ├── BAM_Markdown
    │   ├── *dupmark.bam
    │   ├── *dupmark.bai
    ├── BAM_NODUP
    │   ├── *nodup.bam
    │   ├── *nodup.bai
    │   ├── *nodup.bai
    │   ├── BAM_NFR
    │   │   ├── *nodup.NFR.bam
    │   │   ├── *nodup.NFR.bai
    │   ├── BAM_MN
    │   │   ├── *nodup.MN.bam
    │   │   ├── *nodup.MN.bai
    ├── BAM_NMSRT (NMSRT = name sorted)
    │   │   ├── *nmsrt.nodup.NFR.bam
    ├── BW_COVERAGE
    │   ├── *.nodup.bs1.bw
    │   ├── *.nodup.bs10.bw
    │   ├── BW_NFR
    │   │   ├── *.NFR.bs1.bw
    │   │   ├── *.NFR.bs10.bw
    │   ├── BW_MN
    │   │   ├── *.MN.bs1.bw
    │   │   ├── *.MN.bs10.bw
    ├── TAGDIRECTORIES
    ├── MACS2_PEAKS
    │   ├── *peaks_NFR
    │   │   ├── *.macs2.callpeak.log
    │   │   ├── *.NFR_peaks.narrowPeak
    ├── Metrics
    │   ├── QC_FASTQC
    │   │   ├── *fastqc.html
    │   │   ├── *fastqc.zip
    │   ├── QC_fastq_Screen
    │   │   ├── *screen.html
    │   │   ├── *screen.txt
    │   ├── QC_FASTQ_AdapterRemoval
    │   │   ├── *settings
    │   │   ├── *pair1.truncated.gz
    │   │   ├── *pair2.truncated.gz
    │   │   ├── *discarded.gz
    │   ├── QC_SAMSTATS
    │   │   ├── *samstat.qc.txt
    │   ├── QC_SAMFLAG
    │   │   ├── *flagstat.qc.txt
    │   ├── QC_IDXSTAT
    │   │   ├── *idxstat.qc.txt
    │   ├── QC_PRESEQ
    │   │   ├── *libraryComplexity.preseq.qc.tab
    │   ├── QC_PBC
    │   │   ├── *qc.pbc_mqc.tsv
    │   ├── QC_PICARD
    │   │   ├── *filt.dup.qc.txt
    │   │   ├── *libraryComplexity.picardMetrics.qc.tab
    │   │   ├── *InsertSize.picardMetrics.qc.tab
    │   │   ├── *InsertSize.picardMetrics.qc.pdf
    │   │   ├── *QByCycle.picardMetrics.qc.txt
    │   │   ├── *QByCycle.picardMetrics.qc.pdf
    │   │   ├── *QScoreDist.picardMetrics.qc.pdf
    │   │   ├── *QScoreDist.picardMetrics.qc.txt
    │   │   ├── *alignmentSummary.picardMetrics.qc.txt
    │   ├── QC_BAMQC
    │   │   ├── *qc.qualimap.bamqc
    │   ├── QC_Deeptools
    │   │   ├── *qc.fingerprints.txt
    │   │   ├── *qc.fingerprints.png
    │   │   ├── *qc.coveragePlot.png
    │   │   ├── *qc.estimatereadfiltering.txt
    │   │   ├── *qc.matrix.scaled.gz
    │   │   ├── *qc.plotprofile.txt
    │   │   ├── *qc.plotprofile.png
    │   ├── QC_PHANTOMPEAKS
    │   │   ├── *phantom.qc.spp.out
    │   │   ├── *phantom.qc.top.spp.out
    │   │   ├── *phantom.qc.pdf
    │   │   ├── *phantom.qc.spp.Rdata
    │   ├── QC_MACS
    │   │   ├── *.macs.qc.FRiP_mqc.tsv

## Creating a Multiqc report:

If you wish to create a report of the collected metrics, run the
following in a ubuntu-terminal job with modules:

    # load multiQC
    module load MultiQC
    # Run multiqc in the directory with all the analysis folders:
    multiqc -c /work/Refereneces/Multiqc/multiqc_config_preseq_human.yaml ./  

**Note:** the yaml config file (-c) is optional and it is design to
adjust the genome coverage scale in a plot from Preseq. In case of a
study in mouse, there is a mouse version of the yaml file in the same
directory.

**Notes:**  
- Ensure that the necessary modules are available on your cluster.  
- The script includes Slurm directives to specify resource requirements.
Review and customize the script based on your specific requirements.  
- For additional information on individual tools and parameters, refer
to their official documentation.
