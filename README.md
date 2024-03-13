# Pipeline for bulk ATAC-seq in UCloud using Slurm

## Overview of the Pipeline

This pipeline facilitates the analysis of bulk ATAC-seq data in UCloud HPC using Ubuntu-Terminal with Slurm workload manager: 
1. **pe_fastq_qc.sh:** it is designed to perform quality control (QC) on
paired-end FASTQ files.
2. **pe_align_BWA_sampe_multigenome.sh:** it is designed to align paired-end FASTQ files using BWA sampe. 
3. **pe_postalign_ATACseq_multigenome.sh:** it is designed to process BAM files. It performs mapping QC, and data
preprocessing steps to prepare the data for downstream analysis, for example filtering-split of reads into nucleosomal-free (NFR) and mononucleosomal (MN) among others.
4. **reorganize_files_atacseq.sh:** it is designed to organize files by category.

**Supported genomes:** hg38, mm10, mm39.

## Access to the guides

1. [**Quality Control of fastq files**]()
Runs FASTQC, Fastq_Screen and AdapterRemoval2.

2. [**Alignment of ATAC-seq data using BWA (optionally Salmon)**]()

3. [**Postalignment of BAM files**]()
Read counting is done by FeatureCounts.

4. [**Organize files by category**]()

## General Usage

1.  Clone Repository and copy the script to your Scripts folder
   
<!-- -->
    git clone <repository-url> 
    cd <repository-directory> 

2.  Modify SLURM Parameters (Optional): Open a script
    (**script.sh**) and modify SLURM parameters at the beginning of
    the file, such as account, output file, email notifications, nodes,
    memory, CPU cores, and runtime. Alternatively, you can modify these
    parameters on-the-fly when executing the script.

3.  On UCloud, start a **Terminal Ubuntu** run:

    - Enable **Slurm cluster**
    - To process several samples consider requesting nodes \> 1
    - Set the modules path to **FGM \> Utilities \> App \> easybuild**

![](./Img/terminal_slurm.png)

- Include the References folder **FGM \> References \> References**

![](./Img/terminal_folders.png)

- Include your Scripts folder and the folder with the fastq.gz/bam files.

- **Notes:**
  - Make sure the scripts have executing permission. If not run: `chmod 700 script.sh`
  - Match the job CPUs to the amounts requested in the script.
  - If you modify the memory parameter in the script, specify 5-10% less
    than the memory available in the terminal run.
  - Although it is not necessary to enable **tmux**, it is a good
    practise to always do it.
  - The configuration file of Fastq_Sreen is also located in the
    /References folder.

4.  **Run the Script:** Submit the script to the SLURM cluster:

        sbatch -J <job_name> path_to/Scripts_folder/script.sh <input--file> 

    Replace **input-file** with the full path to your file. 

    For several samples you can use a for loop:

        for i in *<file-pattern>; do sbatch -J <job_name> path_to/Scripts_folder/script.sh $i; sleep 1; done

5.  **Monitor Job:** You can monitor the job using the SLURM commands,
    such as squeue, scontrol show job <job-id>, and check the log files
    generated.

**Notes:** Find test data in UCloud at the FGM project (Utilities/Example_data/bulkRNA/Fastq)

## Contributing
If you find mistakes, errors, or just have suggestions, please create an issue here on GitHub.