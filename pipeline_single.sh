## 2018.02.08
## ATAC-seq FASTQ processing pipeline for single-end reads
## "./pipeline_single.sh ${file_prefix}" will process ${file_prefix}.fastq.gz.  

file_prefix=$1

fastq=${file_prefix}.fastq.gz

trim_prog_path="/path/Trimmomatic-0.32/"
picard_prog_path="/path/picard-tools-1.88/"
adapter="/path/adapter.fa"
FASTA="/path/mm9.fa"
bowt_mm9_path="/path/mm9_bowtie2/"
tmpdir="/path/tmp"

## 1. Trim adapter (Trimmomatic)
echo "1. Trim Adapter (by Trimmomatic)"
java -jar ${trim_prog_path}trimmomatic-0.32.jar SE -trimlog ${file_prefix}_trimmomatic.log ${fastq} ${file_prefix}_trimmed ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
gzip -f ${file_prefix}_trimmomatic.log

##2. Map to mm9 (bowtie2)
echo "2. Map to mm9 (bowtie2)"
fastq="${file_prefix}_trimmed"
bamfile="${file_prefix}_trimmed"

bowtie2 -x ${bowt_mm9_path}mm9 -U $fastq | samtools view -bhS - > ${bamfile}.bam

## Count # of reads per chromosome (chr1-19,X,Y,M)
samtools sort ${bamfile}.bam ${bamfile}.sorted
samtools index ${bamfile}.sorted.bam
samtools idxstats ${bamfile}.sorted.bam > ${bamfile}.sorted.idxstats.txt

rm ${bamfile}.bam

##3. Mark Duplicate (picard)
echo "3. Mark Duplicates (by Picard)"
java -Xmx4g -Djava.io.tmpdir=${tmpdir} -jar ${picard_prog_path}MarkDuplicates.jar INPUT=${bamfile}.sorted.bam OUTPUT=${bamfile}.rmdup.bam METRICS_FILE=${bamfile}.rmdup.bam.picard_metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

samtools sort ${bamfile}.rmdup.bam ${bamfile}.rmdup.sorted
samtools index ${bamfile}.rmdup.sorted.bam
samtools idxstats ${bamfile}.rmdup.sorted.bam > ${bamfile}.rmdup.sorted.bam.idxstats.txt

rm ${bamfile}.rmdup.bam

##4. Filter by mapping quality
echo "4. Filter by Mapping Quality (MQ30)"
bamfile2=${bamfile}.rmdup

samtools view -bhq 30 ${bamfile2}.sorted.bam > ${bamfile2}.MQ30.bam

samtools sort ${bamfile2}.MQ30.bam ${bamfile2}.MQ30.sorted
samtools index ${bamfile2}.MQ30.sorted.bam
samtools idxstats ${bamfile2}.MQ30.sorted.bam > ${bamfile2}.MQ30.sorted.bam.idxstats.txt

rm ${bamfile2}.MQ30.bam

##5. Separate GENOME & chrM
echo "5. Separate GENOME & chrM"
samtools view -bh ${bamfile2}.MQ30.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${bamfile2}.MQ30.sorted.genome.bam
samtools view -bh ${bamfile2}.MQ30.sorted.bam chrM > ${bamfile2}.MQ30.sorted.chrM.bam

samtools sort ${bamfile2}.MQ30.sorted.genome.bam ${bamfile2}.MQ30.genome.sorted
samtools index ${bamfile2}.MQ30.genome.sorted.bam
rm ${bamfile2}.MQ30.sorted.genome.bam

samtools sort ${bamfile2}.MQ30.sorted.chrM.bam ${bamfile2}.MQ30.chrM.sorted
samtools index ${bamfile2}.MQ30.chrM.sorted.bam
rm ${bamfile2}.MQ30.sorted.chrM.bam
