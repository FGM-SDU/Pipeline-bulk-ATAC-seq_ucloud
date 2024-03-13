#!/bin/bash 


# ===============================================
# Create Directories
# ================================================

Metrics=Metrics && mkdir -p "${Metrics}"
BAM_NODUP=BAM_NODUP && mkdir -p "${BAM_NODUP}"
BAM_MarkDup=BAM_MarkDup && mkdir -p "${BAM_MarkDup}"
BAM_NMSRT=BAM_NMSRT && mkdir -p "${BAM_NMSRT}"
BAM_Raw=BAM_Raw && mkdir -p "${BAM_Raw}"
FASTQ_Raw=FASTQ_Raw && mkdir -p "${FASTQ_Raw}"
TAGDIRECTORIES=TAGDIRECTORIES && mkdir -p "${TAGDIRECTORIES}"
BW_COVERAGE=BW_COVERAGE && mkdir -p "${BW_COVERAGE}"
MACS2_PEAKS=MACS2_PEAKS && mkdir -p "${MACS2_PEAKS}"

mv *_S*/*.qc.* "${Metrics}"
mv *_S*/*picard*.error.log* "${Metrics}"
mv *_S*/*.bw "${BW_COVERAGE}"	# Bamcoverage File
mv *_S*/*nmsrt*.ba* "${BAM_NMSRT}"	# bam and bai files sorted by name
mv *_S*/*nodup.ba* "${BAM_NODUP}"	# bam and bai files without duplicates
mv *_S*/*nodup.NFR.ba* "${BAM_NODUP}"	# bam and bai files of nucleosomal free regions
mv *_S*/*nodup.MN.ba* "${BAM_NODUP}"	# bam and bai files of mononucleosomal regions

mv *_S*/*_peaks_NFR/* "${MACS2_PEAKS}"	# folder with peaks called by macs2 
mv *_S*/*callpeak.log* "${MACS2_PEAKS}"	# folder with peaks called by macs2 

mv *_S*/*dupmark.ba* "${BAM_MarkDup}"	# bam and bai files with duplicates FLAGed but not removed
mv ./*.ba* "${BAM_Raw}"	# sorted bam and bai files obtained after alignment and used as input in Post-Align proccessing
mv *_S*/*TagDirectory "${TAGDIRECTORIES}"	# Move TagDirectory directories
mv ./*fastq.gz "${FASTQ_Raw}"	# Move Raw fastq files

mv *_S*/QC_AdapterRemoval/* "${Metrics}"
mv *_S*/QC_FASTQC/* "${Metrics}"
mv *_S*/QC_fastq_Screen/* "${Metrics}"

mv *.discarded.gz "${FASTQ_Raw}" # discarded FASTQ reads by AdapterRemoval
mv *.truncated.gz "${FASTQ_Raw}" # trimmed FASTQ reads #  # .singleton. containing reads where one mate was discarded by AdapterRemoval
mv *.settings "${FASTQ_Raw}"

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

cd ${BAM_NODUP}

BAM_NFR=BAM_NFR && mkdir -p "${BAM_NFR}"
BAM_MN=BAM_MN && mkdir -p "${BAM_MN}"

mv *nodup.NFR.ba* "${BAM_NFR}"
mv *peaks* "${BAM_NFR}"	
mv *nodup.MN.ba* "${BAM_MN}"

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

cd ../${BW_COVERAGE}

BW_NFR=BW_NFR && mkdir -p "${BW_NFR}"
BW_MN=BW_MN && mkdir -p "${BW_MN}"

mv *.NFR.bs*.bw "${BW_NFR}"
mv *.MN.bs*.bw "${BW_MN}"

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

cd ../${Metrics}

QC_AdapterRemoval=QC_AdapterRemoval && mkdir -p "${QC_AdapterRemoval}"
QC_FASTQC=QC_FASTQC && mkdir -p "${QC_FASTQC}"
QC_fastq_Screen=QC_fastq_Screen && mkdir -p "${QC_fastq_Screen}"

mv *.discarded.gz "${QC_AdapterRemoval}" # discarded FASTQ reads by AdapterRemoval
mv *.truncated.gz "${QC_AdapterRemoval}" # trimmed FASTQ reads #  # .singleton. containing reads where one mate was discarded by AdapterRemoval
mv *.settings "${QC_AdapterRemoval}" # settings and summary statistics by AdapterRemoval
mv *fastqc* "${QC_FASTQC}" # QC file "fastqc" on fastq file before alignment
mv *screen* "${QC_fastq_Screen}" # QC file "fastq_screen" on fastq file before alignment


QC_PBC=QC_PBC && mkdir -p "${QC_PBC}"
QC_SAMFLAG=QC_SAMFLAG && mkdir -p "${QC_SAMFLAG}"
QC_SAMSTATS=QC_SAMSTATS && mkdir -p "${QC_SAMSTATS}"
QC_IDXSTAT=QC_IDXSTAT && mkdir -p "${QC_IDXSTAT}"
QC_PRESEQ=QC_PRESEQ && mkdir -p "${QC_PRESEQ}"
QC_PICARD=QC_PICARD && mkdir -p "${QC_PICARD}"
QC_Deeptools=QC_Deeptools && mkdir -p "${QC_Deeptools}"
QC_BAMQC=QC_BAMQC && mkdir -p "${QC_BAMQC}"
QC_PHANTOMPEAKS=QC_PHANTOMPEAKS && mkdir -p "${QC_PHANTOMPEAKS}"
QC_MACS=QC_MACS && mkdir -p "${QC_MACS}"

mv *flagstat* "${QC_SAMFLAG}" # QC file "samtools flagstat" on bam file input after alignment "INITIAL STATS"
mv *.samstat.* "${QC_SAMSTATS}"	# QC file "samtools samstats"
mv *idxstat* "${QC_IDXSTAT}"	# QC file "samtools idxstats"
mv *preseq* "${QC_PRESEQ}" 		# QC file "preseq lc_extrap"
mv *.filt.dup.qc.txt* "${QC_PICARD}"	# QC file picard Markduplicates
mv *picard* "${QC_PICARD}"	# QC files picard metris: CollectAlignmentSummaryMetrics; QualityScoreDistribution; MeanQualityByCycle; **CollectInsertSizeMetrics; **EstimateLibraryComplexity (**)Only PE data
mv *pbc* "${QC_PBC}"	# QC file library complexity 
							# "TotalReadPairs \tDistinctReadPairs \tOneReadPairTwoReadPairs \tNRF=Distinct/Total \tPBC1=OnePair/Distinct \tPBC2=OnePair/TwoPair"
mv *coveragePlot* "${QC_Deeptools}"	# QC files deeptools: plotFingerprint; plotCoverage; estimateReadFiltering.py
mv *fingerprints* "${QC_Deeptools}" # QC files deeptools: plotFingerprint; plotCoverage; estimateReadFiltering.py
mv *estimatereadfiltering* "${QC_Deeptools}" # QC files deeptools: plotFingerprint; plotCoverage; estimateReadFiltering
mv *bamqc "${QC_BAMQC}" # QC files bamqc from qualimap
mv *phantom* "${QC_PHANTOMPEAKS}" # QC files phantompeakqualtools
mv *macs* "${QC_MACS}" # QC file with #peaks and FRiP scores 

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*
chmod 770 ./*/*/*/*
chmod 770 ./*/*/*/*/*

echo "Done!"
