#!/bin/bash -l

#SBATCH -A mandrup    # Names the account for tracking usage, will stay as your lab name
#SBATCH --output slurm-%x.%A.%a.log   # Names an output file for what would go to STDOUT, %x %A %a represent jobname jobid jobarrayindex  
#SBATCH --mail-user victorg@bmb.sdu.dk   # Names an address to be emailed when the job finishes
#SBATCH --mail-type END,FAIL,ARRAY_TASKS  # Specifies when the job finishes (either correctly or failed)
#SBATCH --job-name this_job   # Gives the job a name, so you can find its log file & see it in the queue status, etc
#SBATCH --nodes 1         # How many nodes to be used, ie one compute nodes
#SBATCH --mem 90G        # The job can use up to __GB of ram. It is mutually exclusive with --mem-per-cpu.
#SBATCH -p CLOUD       # Names to use the serial partition
#SBATCH --cpus-per-task 16    # How many cores on that node
##SBATCH --mem-per-cpu 2500M   # he job can use up to 2.5 GB (non-integers are not allowed) of ram per cpu, i.e. 80 GB ram. NOT IN USE
#SBATCH -t 1-20:6:30       # Means to use up to days-hours:minutes:seconds of run time before it will get killed

# run in u1-standard-16 node

# Start runtime
START=$(date +%s)
echo -e "\nStarting processing"

#we set OMP_NUM_THREADS to the number of available cores
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
SAM_THREADS=8 # Use only 8 threads for samtools sort 

# If using #SBATCH --cpus-per-task then
#we set MEM_PER_THREADS to the max memory per CPU
#echo "Running on $SLURM_MEM_PER_CPU M max RAM per CPU"
#export MEM=$((${SLURM_CPUS_PER_TASK}*${SLURM_MEM_PER_CPU}))
#echo "Running on ${MEM} M max RAM"


# else using #SBATCH --mem then
echo "Running on $SLURM_MEM_PER_NODE M max RAM"

# One should only specify 80-90% of available memory 
MEM=$((${SLURM_MEM_PER_NODE}/1024)) # Memory in GB
MEM=$((9*${MEM}/10)) # Leave 19 GB ram aside of the available RAM
MEM_PER_THREAD=$((${MEM}/${SLURM_CPUS_PER_TASK}))
MEM_SAMTOOLS=$((${MEM}/${SAM_THREADS})) # Memory per thread for samtools sort

echo "${OMP_NUM_THREADS} CPUs per task"
echo "${MEM} G total mem"
echo "${MEM_PER_THREAD} G mem per thread"
echo "Samtools sorting conditions: ${SAM_THREADS} threads and ${MEM_SAMTOOLS} G Mem Per Thread"

## Read optional arguments    
function usage {
  echo -e "\n Usage:$(basename $0) -g <genome>  <input_files>"
  echo "Options:"
  echo " -g <genome>   - Specify the genome (mm10, mm39, hg38)"
  echo " -h            - Display this help message"
  exit 1
}

while getopts g:h opt; do
    case "${opt}" in
      g) GENOME="${OPTARG}"
      ;;
      h)          
        usage     
        exit 0              
      ;;
      \?) ## Invalid option
      echo "Invalid option: -${OPTARG}" >&2
      usage
      exit 1
      ;;
      :) ## Argument required
      echo "Option -${OPTARG} requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done
shift "$((OPTIND-1))"  #This tells getopts to move on to the next argument.

if [ -z "${GENOME}" ]; then
  echo "Genome option (-g) is required."
  usage
  exit 1
fi

# PATH to the reference genomes
if [ "${GENOME}" == "mm10" ]; then
  REFERENCE=/work/References/Mouse/mm10/mm10.fa
  BLACKLIST=/work/References/Blacklist/Blacklist_Mouse/mm10.blacklist.v2_merge1000.bed
  GENE=/work/References/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf

elif [ "${GENOME}" == "mm39" ]; then
  REFERENCE="/work/References/Mouse/mm39/mm39.fa"
  BLACKLIST=/work/References/Blacklist/Blacklist_Mouse/mm39.blacklist.v2_merge1000.bed
  GENE=/work/References/GENCODE/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf

elif [ "${GENOME}" == "hg38" ]; then
  REFERENCE=/work/References/Human/hg38_analysisSet/hg38.analysisSet.fa
  BLACKLIST=/work/References/Blacklist/Blacklist_Human/hg38.blacklist.v2_merge1000.bed
  GENE=/work/References/GENCODE/Gencode_human/release_44/gencode.v44.annotation.gtf

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

echo "Using reference: ${REFERENCE}"

# Input Files
if [ $# -eq 0 ]; then
  echo "No input files provided."
  usage
fi

RAW_BAM_FILE=${1?Missing input bam file}
RAW_BAM_FILE=$(readlink -f "${RAW_BAM_FILE}")

echo "Input file: ${RAW_BAM_FILE}"

# ===============================================
# Create Main Output Directory
# ================================================

NAME=$(basename "${RAW_BAM_FILE}" .bam)

if [ -f "${NAME}" ]; then
  >&2 echo "Error: Output location (${NAME}) already exists as a file"
  exit 1
fi

if [ -d "${NAME}" ]; then
  echo "Warning: Output location (${NAME}) is already a directory, reusing, could overwrite"
  # If you don't want to reuse, you could make it exit 1 here to kill the script if
  # the folder already exists
else
  mkdir "${NAME}"
fi

cd "${NAME}"

# Aditional input files
HEADER_PBC=/work/References/Multiqc/pbc_libcomp_header.txt
HEADER_SPP_CORR=/work/References/Multiqc/spp_correlation_header.txt
HEADER_SPP_METRICS=/work/References/Multiqc/spp_rsc_nsc_coefficients_header.txt
HEADER_MACS=/work/References/Multiqc/frip_score_header.txt

# Output Files

BAM_FILE_SAMSTATS="$(basename "${RAW_BAM_FILE}" bam)samstat.qc.tab" # QC file for multiQC
BAM_FILE_IDXSTATS="$(basename "${RAW_BAM_FILE}" bam)idxstat.qc.tab" # QC file for multiQC
BAM_FILE_MAPSTATS="$(basename "${RAW_BAM_FILE}" bam)flagstat.qc.txt" # QC file

ALIGMENTSUMMARY_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)alignmentSummary.picardMetrics.qc.tab" # QC file
QSCOREDIST_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QScoreDist.picardMetrics.qc.tab" # QC file
CHART_QSCOREDIST_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QScoreDist.picardMetrics.qc.pdf" # QC graph
QBYCYCLE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QByCycle.picardMetrics.qc.tab" # QC file
CHART_QBYCYCLE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QByCycle.picardMetrics.qc.pdf" # QC graph

# Only for Pair-End
INSERTSIZE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)InsertSize.picardMetrics.qc.tab" # QC file
CHART_INSERTSIZE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)InsertSize.picardMetrics.qc.pdf" # QC graph
LIBRARYCOMPLEXITY_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)libraryComplexity.picardMetrics.qc.tab" # QC file

# Library complexity
PBC_FILE_QC="$(basename "${RAW_BAM_FILE}" bam)qc.pbc_mqc.tsv" # QC file
# Library complexity description
# TotalReadPairs(MT)	DistinctReadPairs(M1)	OneReadPair	TwoReadPairs(M2)	NRF=Distinct/Total	PBC1=OnePair/Distinct	PBC2=OnePair/TwoPair

LIBRARYCOMPLEXITY_PRESEQ="$(basename "${RAW_BAM_FILE}" bam)libraryComplexity.preseq.qc.tab" # QC file

# Temporary (Tmp) files: transitory files in the proccesses of fixing mates in pair-end data and filtering contigs (remove of chrMT, chrUn, etc.)
TMP_QFILT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.bam"
TMP_NMSRT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.nmsrt.bam"
TMP_QFILT_BED_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.bed"
SAM_FILE_TMP="$(basename "${RAW_BAM_FILE}" bam)filt_contigs.sam"
BAM_FILE_CHR_FILETERED_TMP="$(basename "${RAW_BAM_FILE}" bam)filt_contigs.bam"

# Output Files of Deeptools
FINGERPRINTS="$(basename "${RAW_BAM_FILE}" bam)qc.fingerprints.txt" # QC file
CHART_FINGERPRINTS="$(basename "${RAW_BAM_FILE}" bam)qc.fingerprints.png" # QC graph
CHART_COVERAGE="$(basename "${RAW_BAM_FILE}" bam)qc.coveragePlot.png" # QC graph
BIGWIG_COVERAGE_PER_BASE="$(basename "${RAW_BAM_FILE}" bam)nodup.bs1.bw" # Bigwig from BAM file without duplicates. Bin size 1
BIGWIG_COVERAGE_BIN="$(basename "${RAW_BAM_FILE}" bam)nodup.bs10.bw" # Bigwig from BAM file without duplicates. Bin size 10. Normalization RPKM
ESTIMATE_READ_FILTERING="$(basename "${RAW_BAM_FILE}" bam)qc.estimatereadfiltering.txt" # QC file
MATRIX_SCALED_GENES="$(basename "${RAW_BAM_FILE}" bam)qc.matrix.scaled.gz"
PLOTPROFILE_DATA="$(basename "${RAW_BAM_FILE}" bam)qc.plotprofile.txt" # QC file
PLOTPROFILE_PLOT="$(basename "${RAW_BAM_FILE}" bam)qc.plotprofile.png" # QC file

# Output Files for Markduplicates
MARKDUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)raw.dupmark.bam"	# Intermediate file later removed
NMSRT_MARKDUP_QFILT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.dupmark.nmsrt.bam" # Intermediate file later removed	
FINAL_MARKDUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)dupmark.bam"	# Bam file with duplicates Flagged for allow a downstream program to remove/ignore duplicates or keep them
FINAL_MARKDUP_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)dupmark.bai"	# Index file
DUP_FILE_QC="$(basename "${RAW_BAM_FILE}" bam)filt.dup.qc.txt"	# QC file


# Output Files of Filtering
NODUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.bam"
NODUP_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.bai"

# Output Files of Phantompeakqualtools
PHANTOMPEAKS_METRICS="$(basename "${RAW_BAM_FILE}" bam)phantom.qc.spp.out" # QC file
PHANTOMPEAKS_METRICS_TOP="$(basename "${RAW_BAM_FILE}" bam)phantom.qc.top.spp.out" # QC file
PHANTOMPEAKS_PLOT="$(basename "${RAW_BAM_FILE}" bam)phantom.qc.pdf"	# QC graph
PHANTOMPEAKS_RDATA="$(basename "${RAW_BAM_FILE}" bam)phantom.qc.spp.Rdata" # QC Rdata

# Split Files ATAC
TMP_BAM_NFR="$(basename "${RAW_BAM_FILE}" bam)nodup.NFR.tmp.bam" # nuclesome-free
TMP_BAM_MN="$(basename "${RAW_BAM_FILE}" bam)nodup.MN.tmp.bam" # mono-nuclesome

BAM_NFR="$(basename "${RAW_BAM_FILE}" bam)nodup.NFR.bam" # nuclesome-free
BAM_MN="$(basename "${RAW_BAM_FILE}" bam)nodup.MN.bam" # mono-nuclesome

BAM_NFR_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.NFR.bai"
BAM_MN_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.MN.bai"

NMSRT_BAM_NFR="$(basename "${RAW_BAM_FILE}" bam)nmsrt.nodup.NFR.bam"


BIGWIG_COVERAGE_PER_BASE_NFR="$(basename "${RAW_BAM_FILE}" bam)NFR.bs1.bw" # nuclesome-free. Bin size 1 
BIGWIG_COVERAGE_BIN_NFR="$(basename "${RAW_BAM_FILE}" bam)NFR.bs10.bw"  # nuclesome-free. Bin size 10. Normalization RPKM 

BIGWIG_COVERAGE_PER_BASE_MN="$(basename "${RAW_BAM_FILE}" bam)MN.bs1.bw" # mono-nuclesome. Bin size 1 
BIGWIG_COVERAGE_BIN_MN="$(basename "${RAW_BAM_FILE}" bam)MN.bs10.bw"  # mono-nuclesome. Bin size 10. Normalization RPKM 

# ===============================================
# Flagstat
# Remove low MAPQ reads
# Compute Library Complexity
# ================================================
#Sort by name
#convert to bedPE and obtain fragment coordinates
#sort by position and strand
#Obtain unique count statistics

# Commands

module load SAMtools

echo "Samtools QC..."
START_SUBPROCESS=$(date +%s)

samtools stats ${RAW_BAM_FILE} > ${BAM_FILE_SAMSTATS}
samtools idxstats ${RAW_BAM_FILE} > ${BAM_FILE_IDXSTATS}
samtools flagstat ${RAW_BAM_FILE} > ${BAM_FILE_MAPSTATS}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload SAMtools
module load Qualimap

echo "qualimap..." 
START_SUBPROCESS=$(date +%s)

qualimap bamqc -bam ${RAW_BAM_FILE} -c -outdir ./${NAME}.qc.qualimap.bamqc -nt "${OMP_NUM_THREADS}" --java-mem-size="${MEM}G"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload Qualimap
module load SAMtools R picard BEDTools

echo "Picard QC..."
START_SUBPROCESS=$(date +%s)

java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics -R ${REFERENCE} -I ${RAW_BAM_FILE} -O ${ALIGMENTSUMMARY_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.alignmentSummary.error.log"
java -jar $EBROOTPICARD/picard.jar QualityScoreDistribution --CHART ${CHART_QSCOREDIST_PICARDMETRICS} -I ${RAW_BAM_FILE} -O ${QSCOREDIST_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.QScoreDist.error.log"
java -jar $EBROOTPICARD/picard.jar MeanQualityByCycle -I ${RAW_BAM_FILE} -O ${QBYCYCLE_PICARDMETRICS} -CHART ${CHART_QBYCYCLE_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.QByCycle.error.log"

#Only Pair-End
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -H ${CHART_INSERTSIZE_PICARDMETRICS} -M 0.5 -I ${RAW_BAM_FILE} -O ${INSERTSIZE_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.InsertSizeMetrics.error.log"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Quality filtering..."
START_SUBPROCESS=$(date +%s)

MAPQ_THRESH=30
samtools view -q ${MAPQ_THRESH} -@ "${OMP_NUM_THREADS}" -Sbh ${RAW_BAM_FILE} -o ${TMP_QFILT_BAM_FILE} 

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Library complexity.."
START_SUBPROCESS=$(date +%s)
# Get PBC bottlenecking metrics
  # mt: number of reads (TotalReadPairs)
  # m0: number of all genomic locations where reads mapped (DistinctReadPairs)
  # m1: number of genomic locations where only one read maps uniquely (OneReadPair)
  # m2: number of genomic locations where 2 reads map uniquely (TwoReadPairs)
  # NRF: Non-Redundant Fraction
  # PBC1: PCR Bottlenecking Coefficient 1
  # PBC2: PCR Bottlenecking Coefficient 2  

#echo -e "TotalReadPairs \tDistinctReadPairs \tOneReadPair \tTwoReadPairs \tNRF=Distinct/Total \tPBC1=OnePair/Distinct \tPBC2=OnePair/TwoPair" > ${PBC_FILE_QC}
samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -n ${TMP_QFILT_BAM_FILE} -o ${TMP_NMSRT_BAM_FILE} # Will produce name sorted BAM

# Convert bam to bedpe, exclude chrM, sort and count number of occurrences of each unique entry. Finally awk calculte the metrics.
bedtools bamtobed -i ${TMP_NMSRT_BAM_FILE} -bedpe | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | \
sort | uniq -c | awk -v NAME="${NAME}" 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",NAME,mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'| cat ${HEADER_PBC} - > "${PBC_FILE_QC}" 

# Only Pair-End
#java -jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity -I ${TMP_QFILT_BAM_FILE} -O ${LIBRARYCOMPLEXITY_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.EstimateLibraryComplexity.error.log"

# Create a bed file for Preseq run more stable than from bam
(bedtools bamtobed -bedpe -i ${TMP_NMSRT_BAM_FILE} 2>/dev/null) | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 - > "${TMP_QFILT_BED_FILE}"

# Preseq
module load preseq 

preseq lc_extrap -pe -o ${LIBRARYCOMPLEXITY_PRESEQ} ${TMP_QFILT_BED_FILE} # computation in Bam files is less stable

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# ===================================
# Mark duplicates
# Insert Size
# ====================================
echo "Fixing Coodinates and sorting: filter chrM, random and chrUn..."
START_SUBPROCESS=$(date +%s)

samtools view -@ "${OMP_NUM_THREADS}" -Sh ${TMP_QFILT_BAM_FILE} | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | sed '/chrM/d;/random/d;/chrUn/d' > ${SAM_FILE_TMP}
samtools view -@ "${OMP_NUM_THREADS}" -Shb ${SAM_FILE_TMP} | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${BAM_FILE_CHR_FILETERED_TMP}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Marking duplicates..."
START_SUBPROCESS=$(date +%s)

java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ${BAM_FILE_CHR_FILETERED_TMP} -O ${MARKDUP_BAM_FILE} -M ${DUP_FILE_QC} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.MarkDuplicates.error.log"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Filtering..."
echo "Filtering 780 (keep duplicates and all paired reads, Remove unmapped), Fixming mate, Sorting and indexing..."
START_SUBPROCESS=$(date +%s)
# ===========================================================================================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Obtain name sorted BAM file
# Obtain position sorted BAM
# Index final position sorted BAM
# Create final name sorted BAM
# ===================================================================================================
# Filter out read and mate unmapped, not primary alignment (multimapped), read fails platform/vendor quality checks 
samtools view -F 780 -Shb ${MARKDUP_BAM_FILE} -@ "${OMP_NUM_THREADS}" | samtools sort -n -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${NMSRT_MARKDUP_QFILT_BAM_FILE}
# create mate-cigar strings in SAM/BAM files
samtools fixmate -r -@ "${OMP_NUM_THREADS}" ${NMSRT_MARKDUP_QFILT_BAM_FILE} - | samtools view -@ "${OMP_NUM_THREADS}" -F 780 -Shb - | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${FINAL_MARKDUP_BAM_FILE}
samtools index ${FINAL_MARKDUP_BAM_FILE} ${FINAL_MARKDUP_INDEX_FILE}

rm ${TMP_NMSRT_BAM_FILE} ${TMP_QFILT_BED_FILE} ${MARKDUP_BAM_FILE} ${TMP_QFILT_BAM_FILE} ${SAM_FILE_TMP} ${BAM_FILE_CHR_FILETERED_TMP} ${NMSRT_MARKDUP_QFILT_BAM_FILE}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Filtering -F 1804: remove duplicates, remove unmapped, mate unmapped not primary aligment and reads failing platform/vendor quiality checks..."
echo "Filtering -f 2: Keep only properly paired reads (read mapped in proper pair)..."
START_SUBPROCESS=$(date +%s)
# ===========================================================================================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
# Obtain position sorted BAM
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ===================================================================================================
samtools view -@ "${OMP_NUM_THREADS}" -F 1804 -f 2 -Sbh ${FINAL_MARKDUP_BAM_FILE} | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${NODUP_BAM_FILE} # Will produce coordinate sorted BAM / Final NODUP BAM file

echo "Indexing..."
samtools index ${NODUP_BAM_FILE} ${NODUP_INDEX_FILE}		# Index Final NODUP BAM file

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."


module load deepTools
echo "Filtering-Split into nucleosomal-free bin 38-100 bp and nucleosomal bin 180-247 bp bin..."
START_SUBPROCESS=$(date +%s)

alignmentSieve -b ${NODUP_BAM_FILE} -o ${TMP_BAM_NFR} -p "${OMP_NUM_THREADS}" --ATACshift --minFragmentLength 38 --maxFragmentLength 100
alignmentSieve -b ${NODUP_BAM_FILE} -o ${TMP_BAM_MN} -p "${OMP_NUM_THREADS}" --ATACshift --minFragmentLength 180 --maxFragmentLength 247

echo "Filtering-Split into nucleosomal-free bin 38-100 bp and nucleosomal bin 180-247 bp bin...Sorting"
samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -o ${BAM_NFR} ${TMP_BAM_NFR}
rm ${TMP_BAM_NFR}
samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -o ${BAM_MN} ${TMP_BAM_MN}
rm ${TMP_BAM_MN}

echo "Filtering-Split into nucleosomal-free bin 38-100 bp and nucleosomal bin 180-247 bp bin...Indexing"
samtools index ${BAM_NFR} ${BAM_NFR_INDEX_FILE}
samtools index ${BAM_MN} ${BAM_MN_INDEX_FILE}

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Name Sorting..."
START_SUBPROCESS=$(date +%s)

samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -n ${BAM_NFR} -o ${NMSRT_BAM_NFR}	# Create final name sorted NODUP_BAM

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Deeptools QC..."
START_SUBPROCESS=$(date +%s)

plotCoverage -b ${FINAL_MARKDUP_BAM_FILE} -o ${CHART_COVERAGE} --ignoreDuplicates -T "CoveragePlot"
estimateReadFiltering -b ${RAW_BAM_FILE} -o ${ESTIMATE_READ_FILTERING} -p "${OMP_NUM_THREADS}" -bl ${BLACKLIST}

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Deeptools QC: plotFingerprint on 10 million mates/pairs subsampled..."
START_SUBPROCESS=$(date +%s)

READS=10000000
# Temporary Files   
TMP_DOWN_BAM="${NAME}.10M.tmp.bam" # 10 million mates/pairs (fragments), i.e. x2 # reads
TMP_DOWN_BAM_INDEX="${NAME}.10M.tmp.bai"

# Calculate FRACTION for each target BAM file
FRACTION=$(samtools idxstats "${FINAL_MARKDUP_BAM_FILE}" | cut -f3 | awk -v ct=${READS} 'BEGIN {total=0} {total += $1} END {print ct/total}')
echo "Fraction for target ${FINAL_MARKDUP_BAM_FILE} = ${FRACTION}"
    
if (( $(echo "${FRACTION} 1" | awk '{print ($1 > $2)}') )); then
  echo "Warning: File (${FINAL_MARKDUP_BAM_FILE}) is low depth (less than 10 M fragments)"
  else     
  echo "Samtools downsampling target bam files..."    
  samtools view -Sbh -@ "${OMP_NUM_THREADS}" -s ${FRACTION} ${FINAL_MARKDUP_BAM_FILE} | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${TMP_DOWN_BAM}
  samtools index ${TMP_DOWN_BAM} ${TMP_DOWN_BAM_INDEX}    
  plotFingerprint -b ${TMP_DOWN_BAM} -plot ${CHART_FINGERPRINTS} --outRawCounts ${FINGERPRINTS} -p "${OMP_NUM_THREADS}" --skipZeros --ignoreDuplicates -T "Fingerprints"
fi  

rm ${TMP_DOWN_BAM} ${TMP_DOWN_BAM_INDEX}

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload deepTools SAMtools R picard BEDTools preseq
module load R/4.1.2-foss-2021b Phantompeakqualtools

echo "Phantompeakqualtools QC..."
START_SUBPROCESS=$(date +%s)

run_spp.R -c=${NODUP_BAM_FILE}	-savd=${PHANTOMPEAKS_RDATA} -savp=${PHANTOMPEAKS_PLOT} -out=${PHANTOMPEAKS_METRICS}

# keep only the top value of COL3: estFragLen: comma separated strand cross-correlation peak(s)
sed -r 's/,[^\t]+//g' ${PHANTOMPEAKS_METRICS} > "${PHANTOMPEAKS_METRICS_TOP}"

# Prepare files for reporting cross-correlation
cp "${HEADER_SPP_CORR}" "${NAME}.phantom.qc.spp.corr_mqc.tsv"
Rscript --max-ppsize=500000 -e "load('${PHANTOMPEAKS_RDATA}'); write.table(crosscorr\$cross.correlation, file=\"${NAME}.phantom.qc.spp.corr_mqc.tsv\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

# Prepare files for reporting NSC and RSC coefficients
SPP_NSC_RSC="${NAME}.phantom.qc.spp.coeff_mqc.tsv"
awk -v OFS='\t' '{print "${NAME}", $9,$10}' ${PHANTOMPEAKS_METRICS_TOP} | cat ${HEADER_SPP_METRICS} - > "${SPP_NSC_RSC}"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload R/4.1.2-foss-2021b Phantompeakqualtools
module load deepTools

echo "Deeptools Coverages Profiles..."
START_SUBPROCESS=$(date +%s)

# bin size 10 bp
# The smooth length defines a window, larger than the binSize, to average the number of reads. For example, if the –binSize is set to 20 and the –smoothLength is set to 60, then, for each bin, the average of the bin and its left and right neighbors is considered.
bamCoverage -b ${FINAL_MARKDUP_BAM_FILE} -o ${BIGWIG_COVERAGE_PER_BASE} -bs 1 -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}" --ignoreDuplicates --samFlagExclude 1804
bamCoverage -b ${FINAL_MARKDUP_BAM_FILE} -o ${BIGWIG_COVERAGE_BIN} -bs 10 -e --normalizeUsing RPKM -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}" --ignoreDuplicates --samFlagExclude 1804

# Coverage track of split files  NFR /MN
bamCoverage -b ${BAM_NFR} -o ${BIGWIG_COVERAGE_PER_BASE_NFR} -bs 1 -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}"
bamCoverage -b ${BAM_NFR} -o ${BIGWIG_COVERAGE_BIN_NFR} -bs 10 -e --normalizeUsing RPKM -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}"

bamCoverage -b ${BAM_MN} -o ${BIGWIG_COVERAGE_PER_BASE_MN} -bs 1 -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}"
bamCoverage -b ${BAM_MN} -o ${BIGWIG_COVERAGE_BIN_MN} -bs 10 -e --normalizeUsing RPKM -bl ${BLACKLIST} -p "${OMP_NUM_THREADS}"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Deeptools plotProfile..."
START_SUBPROCESS=$(date +%s)

computeMatrix scale-regions -R ${GENE} -S ${BIGWIG_COVERAGE_BIN} -o ${MATRIX_SCALED_GENES} --regionBodyLength 1000 -b 3000 -a 3000 -bl ${BLACKLIST} --missingDataAsZero  --skipZeros --smartLabels -p "${OMP_NUM_THREADS}" -q
plotProfile -m ${MATRIX_SCALED_GENES} -out ${PLOTPROFILE_PLOT} --outFileNameData ${PLOTPROFILE_DATA} --plotTitle "Read Distribution Profile"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload deepTools
module load SAMtools Homer

echo "HOMER: Making TagDirectory..."
START_SUBPROCESS=$(date +%s)

makeTagDirectory ./${NAME}_TagDirectory ${BAM_NFR} -single -tbp 1 -genome ${REFERENCE} -checkGC

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload SAMtools Homer
module load SAMtools MACS2

echo "MACS2: Calling peaks..."
START_SUBPROCESS=$(date +%s)

if [ "${GENOME}" == "mm10" ] || [ "${GENOME}" == "mm39" ]; then
  GENOME_SZ="mm"

elif [ "${GENOME}" == "hg38" ]; then
  GENOME_SZ="hs"

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

# Specify the path to the MACS2 output directory
MACS2_OUTPUT_DIR="${NAME}_peaks_NFR"

macs2 callpeak -g ${GENOME_SZ} -f BAMPE --keep-dup all -t ${BAM_NFR} -n "${NAME}.NFR" --outdir ${MACS2_OUTPUT_DIR} --verbose 2 2> "${NAME}.macs2.callpeak.log"

module unload MACS2
module load SAMtools BEDTools

# Extract the total number of peaks from the MACS2 output
TOTAL_PEAKS=$(cat "${MACS2_OUTPUT_DIR}/${NAME}.NFR_peaks.narrowPeak" | wc -l)

### FRiP Score ###
# Calculate reads in peaks using intersectBed
READS_IN_PEAKS=$(bedtools intersect -a "${BAM_NFR}" -b "${MACS2_OUTPUT_DIR}/${NAME}.NFR_peaks.narrowPeak" | samtools view -c)

# Calculate total mapped reads using samtools flagstat
samtools flagstat "${BAM_NFR}" > "${BAM_NFR}.flagstat"
MAPPED_READS=$(grep 'mapped (' "${BAM_NFR}.flagstat" | grep -v "primary" | awk '{print $1}')

# Calculate FRiP score
FRIP_SCORE=$(awk -v a="${READS_IN_PEAKS}" -v b="${MAPPED_READS}" 'BEGIN {print a/b}')

# Write Total peaks and FRiP score to output file for multiqc
MACS_FRIP="${NAME}.macs.qc.FRiP_mqc.tsv"
TMP_MACS_METRICS="${NAME}.macs.qc.txt"

echo "${NAME},${TOTAL_PEAKS},${FRIP_SCORE}"
echo -e "${NAME}\t${TOTAL_PEAKS}\t${FRIP_SCORE}" > "${TMP_MACS_METRICS}"
awk -v OFS='\t' '{print $1,$2,$3}' ${TMP_MACS_METRICS} | cat ${HEADER_MACS} - > "${MACS_FRIP}"

rm "${BAM_NFR}.flagstat" "${TMP_MACS_METRICS}"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

# Finalize
END=$(date +%s)
RUNTIME=$((END-START))
H=$((RUNTIME / 3600 ))  # Calculate hours
M=$(( (RUNTIME / 60 ) % 60 ))  # Calculate minutes
S=$(( RUNTIME % 60 ))  # Calculate seconds
echo -e "\tProcessing completed. Total run time: ${H} hours, ${M} minutes, and ${S} seconds."

echo "Done!"
