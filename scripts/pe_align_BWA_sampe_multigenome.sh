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
#SBATCH -t 20:30:00       # Means to use up to hours:minutes:seconds of run time before it will get killed

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
  echo -e "\n Usage:$(basename $0) -g <genome> <input_files> -t"
  echo "Options:"
  echo " -g <genome>        - Specify the genome (mm10, mm39, hg38)"
  echo " -t                 - Specify to trim (optional)"
  echo " -h                 - Display this help message"
  exit 1
}

while getopts g:t:h opt; do
    case "${opt}" in
      g) GENOME="${OPTARG}"
      ;;
      t) INCLUDE_TRIMMING=true
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

## Input files
if [ $# -eq 0 ]; then
  echo "No input files provided."
  usage
fi

INPUT1=${1?Missing input fastq-r1.gz file}
INPUT2=$(basename "${INPUT1}" _R1_001.fastq.gz)_R2_001.fastq.gz
PREFIX=$(basename "${INPUT1}" _R1_001.fastq.gz)

echo "Input file: ${INPUT1} ${INPUT2}"

# Load programs
module load AdapterRemoval BWA SAMtools picard

# Reference Genomes
if [ -z "${GENOME}" ]; then
  echo "Genome option (-g) is required."
  usage
  exit 1
fi

# PATH to the reference genomes
if [ "${GENOME}" == "mm10" ]; then
  REFERENCE="/work/References/Mouse/mm10/bwa0.7.17/no_haplotypes/mm10.fa"

elif [ "${GENOME}" == "mm39" ]; then
  REFERENCE="/work/References/Mouse/mm39/bwa_0.7.17-r1188/mm39.fa"

elif [ "${GENOME}" == "hg38" ]; then
  REFERENCE="/work/References/Human/hg38_analysisSet/bwa_0.7.17-r1188/hg38.analysisSet.fa"

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

echo "Using reference: ${REFERENCE}"

if [ "${INCLUDE_TRIMMING}" = true ]; then
  echo "fastq files will be trimmed by AdapterRemovalv2"
  # Temporay files
  TRIMMED1=$(basename "${INPUT1}" _R1_001.fastq.gz)_trimmed_R1_001.fastq.gz
  TRIMMED2=$(basename "${INPUT1}" _R1_001.fastq.gz)_trimmed_R2_001.fastq.gz
  # Commands
  echo "AdapterRemoval..."
  AdapterRemoval --threads "${OMP_NUM_THREADS}" --file1 "${INPUT1}" --file2 "${INPUT2}" --basename "${PREFIX}" --gzip --output1 ${TRIMMED1} --output2 ${TRIMMED2}
else
  echo "fastq files will aligned without trimming"
  TRIMMED1="${INPUT1}"
  TRIMMED2="${INPUT2}"
fi

# Alignning
START_SUBPROCESS=$(date +%s)

echo "Aligning and coordinate sorting..."
echo "Running BWA..."

bwa sampe ${REFERENCE} \
  <(bwa aln -t "${OMP_NUM_THREADS}" ${REFERENCE} "${TRIMMED1}") \
  <(bwa aln -t "${OMP_NUM_THREADS}" ${REFERENCE} "${TRIMMED2}") \
  "${TRIMMED1}" \
  "${TRIMMED2}" | \
  samtools view -@ "${OMP_NUM_THREADS}" -Shu - |  \
  samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -O bam -o "${PREFIX}.tmp.bam"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Add or correct the read group information
echo "Picard AddOrReplaceReadGroups..."
START_SUBPROCESS=$(date +%s)

ID=$(zcat "${INPUT1}" | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) ## read group identifier 
LB=$(echo "${INPUT1}" | cut -d "_" -f1,2)                                ## library ID
PL="illumina"                                                           ## platform (e.g. illumina, solid)
PU=$(echo "${ID}"."${LB}")                                                            ##Platform Unit
SM=$(echo "${INPUT1}" | cut -d"_" -f1)                                          ##sample ID

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I "${PREFIX}.tmp.bam" -O "${PREFIX}.bam" --RGID "${ID}" --RGLB "${LB}" --RGPL "${PL}" --RGPU "${PU}" --RGSM "${SM}" --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.AddOrReplaceReadGroups.error.log"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Generate a index of the bam file
echo "Indexing..."
START_SUBPROCESS=$(date +%s)

samtools index -b "${PREFIX}.bam" "${PREFIX}.bai"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Deleting temporary files
if [ "${INCLUDE_TRIMMING}" = true ]; then
  rm "${PREFIX}.tmp.bam" "${TRIMMED1}" "${TRIMMED2}" "${PREFIX}.singleton.truncated.gz" "${PREFIX}.discarded.gz" "${PREFIX}.settings"
else
  rm "${PREFIX}.tmp.bam"
fi

module unload AdapterRemoval BWA SAMtools picard

# Finalize
END=$(date +%s)
RUNTIME=$((END-START))
H=$((RUNTIME / 3600 ))  # Calculate hours
M=$(( (RUNTIME / 60 ) % 60 ))  # Calculate minutes
S=$(( RUNTIME % 60 ))  # Calculate seconds
echo -e "\tProcessing completed. Total run time: ${H} hours, ${M} minutes, and ${S} seconds."

echo "... Done!!!"
