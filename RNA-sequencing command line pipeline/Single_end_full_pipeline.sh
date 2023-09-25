#!/bin/bash

#Directory from where the script is executed
parent_directory=$(dirname $0)

project_name=$1 #name of project (plain input)
organism_name=$2 #short name of organism, needed to build reference when mapping with bowtie2 (plain input)
fastq_gz_dir=$3 #directory with all of the fastq files (path to directory)
fasta_reference=$4 #reference fasta (path to .fna file)
seq_adapter=$5 #Seqencing adapter file (path to .fa file)
gtf_ref=$6 #GTF/GFF3 reference file (path to .gtf/.gff file)

#Check to see if the user input is a directory
if [[ ! -d $fastq_gz_dir ]]; then
  echo "Enter valid directory with fastq files. Exiting..."
  exit 
fi

if [[ -d $parent_directory/$project_name ]]; then
  echo "The directory $parent_directory/$project_name exists, skipping creation of directory..."
else 
  echo "Directory does not exist... creating $parent_directory/$project_name"
  #Creates the directory structure
  mkdir -p $parent_directory/$project_name/{01_fastq,02_trimmed_fastq,03_mapped_reads/{SAM,sBAM},04_gene_counts,Quality_Control/{01_fastq,02_trimmed_fastq,03_mapped_reads},ref_seqs/$organism_name}
fi 

#Pathnames stored as variables to be called on 
QC1=$parent_directory/$project_name/Quality_Control/01_fastq

tfastq=$parent_directory/$project_name/02_trimmed_fastq

QC2=$parent_directory/$project_name/Quality_Control/02_trimmed_fastq

indexed=$parent_directory/$project_name/ref_seqs/$organism_name

mappedSAM=$parent_directory/$project_name/03_mapped_reads/SAM

sortedBAM=$parent_directory/$project_name/03_mapped_reads/sBAM

QC3=$parent_directory/$project_name/Quality_Control/03_mapped_reads

counts=$parent_directory/$project_name/04_gene_counts

ref_seqs=$parent_directory/$project_name/ref_seqs

#number of sample files
fil_num=$(find $fastq_gz_dir -name "*.fastq.gz" | wc -l)

echo "Processing $fil_num files in $fastq_gz_dir"

#########################################################################
####RUNS INITIAL QUALITY CONTROL ON RAW FILES 

#Creates blank .txt files for stdout and stderr
> $QC1/raw_fastqc_stdout.txt
> $QC1/raw_fastqc_stderr.txt

#Runs quality control on all of the raw files 
for fastq_files in $fastq_gz_dir/*.fastq.gz
do
  echo "Running quality control with FastQC on $fastq_files..."

  #run the fastQC analysis for the raw read
  fastqc --extract $fastq_files \
  --outdir $QC1 \
  1>>$QC1/raw_fastqc_stdout.txt \
  2>>$QC1/raw_fastqc_stderr.txt

done 

#########################################################################
###########################Paired End Trimmomatic########################
> $tfastq/trimmomatic_stdout.txt
> $tfastq/trimmomatic_stderr.txt

mate1=""
mate2=""

for fastq_files in $fastq_gz_dir/*.fastq.gz
do

  file_name="$(basename $fastq_files .fastq.gz)"

  cut_name="${file_name: -6}"

  if [[ "$cut_name" == *"R1"* ]]; then
    mate1="$fastq_files"
    basename_mate1="$(basename $mate1)"
    echo "Assigning $basename_mate1 to mate1 in pair"
  fi

  if [[ "$cut_name" == *"R2"* ]]; then
    mate2="$fastq_files"
    basename_mate2="$(basename $mate2)"
    echo "Assigning $basename_mate2 to mate2 in pair"
  fi

  if [ "$mate1" != "" ] && [ "$mate2" != "" ]; then

    mate1_basename="$(basename $mate1)"
    mate2_basename="$(basename $mate2)"

    paired_mate1="$(basename $mate1 .fastq.gz)_paired_trimmed.fastq.gz"
    unpaired_mate1="$(basename $mate1 .fastq.gz)_unpaired_trimmed.fastq.gz"
    paired_mate2="$(basename $mate2 .fastq.gz)_paired_trimmed.fastq.gz"
    unpaired_mate2="$(basename $mate2 .fastq.gz)_unpaired_trimmed.fastq.gz"

    echo "Trimming $mate1_basename into $paired_mate1 and $mate2_basename into $paired_mate2 ..."

    echo "Samples: $mate1_basename and $mate2_basename" | cat >> $tfastq/trimmomatic_stdout.txt
    echo "Samples: $mate1_basename and $mate2_basename" | cat >> $tfastq/trimmomatic_stderr.txt

    trimmomatic PE -threads 4 \
    $mate1 $mate2 \
    $tfastq/$paired_mate1 $tfastq/$unpaired_mate1 $tfastq/$paired_mate2 $tfastq/$unpaired_mate2 \
    ILLUMINACLIP:$seq_adapter:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 \
    1>>$tfastq/trimmomatic_stdout.txt \
    2>>$tfastq/trimmomatic_stderr.txt

    mate1=""
    mate2=""

    echo "
    " | cat >> $tfastq/trimmomatic_stdout.txt

    echo "
    " | cat >> $tfastq/trimmomatic_stderr.txt
  fi
done

#########################################################################
#Runs quality control on all of the trimmed files

> $QC2/trimmed_fastqc_stdout.txt
> $QC2/trimmed_fastqc_stderr.txt

#Fastqc analysis on trimmed files 
for trimmed_files in $tfastq/*.fastq.gz
do
  echo "Running quality control with FastQC on $trimmed_files..."

  fastqc --extract $trimmed_files \
  -outdir $QC2 \
  1>>$QC2/trimmed_fastqc_stdout.txt \
  2>>$QC2/trimmed_fastqc_stderr.txt

done

###########################################################################
######################Build Bowtie2 reference genome#######################
cd $indexed

echo "Building bowtie2 index from $fasta_reference..."

#Builds bowtie2 index within $indexed
bowtie2-build $fasta_reference $organism_name

> $mappedSAM/sam_mapping_stdout.txt
> $mappedSAM/sam_mapping_stderr.txt


########################Paired End Bowtie2 Processing######################
cd $parent_directory

echo "Processing $fil_num files in $fastq_gz_dir"

mate1=""
mate2=""

for files in $tfastq/*_paired_trimmed.fastq.gz
do

  file_name="$(basename $files _paired_trimmed.fastq.gz)"

  cut_name="${file_name: -6}"

  if [[ "$cut_name" == *"R1"* ]]; then
    mate1="$files"
    basename_mate1="$(basename $mate1)"
    echo "Assigning $basename_mate1 to mate1 in pair"
  fi

  if [[ "$cut_name" == *"R2"* ]]; then
    mate2="$files"
    basename_mate2="$(basename $mate2)"
    echo "Assigning $basename_mate2 to mate2 in pair"
  fi

  if [ "$mate1" != "" ] && [ "$mate2" != "" ]; then

    cut_name_file="$(basename $mate1 _R1_001_paired_trimmed.fastq.gz)_mapped.sam"
    basename_mate1="$(basename $mate1)"
    basename_mate2="$(basename $mate2)"

    echo "Converting $basename_mate1 and $basename_mate2 into $cut_name_file SAM file ..."

    echo "Sample: $cut_name_file" | cat >> $mappedSAM/sam_mapping_stdout.txt

    bowtie2 -p 8 -x $organism_name \
    -1 $mate1 -2 $mate2 -S $mappedSAM/$cut_name_file \
    1>>$mappedSAM/sam_mapping_stderr.txt \
    2>>$mappedSAM/sam_mapping_stdout.txt

    echo "
    " | cat >> $mappedSAM/sam_mapping_stdout.txt
    
    mate1=""
    mate2=""
  fi
done 

###########################################################################
#Converts read alignment file (SAM) to Binary Alignment Map (BAM):

#Converts sam files to bam files 
> $sortedBAM/bam_sorting_stdout.txt
> $sortedBAM/bam_sorting_stderr.txt

#Converts sam files to bam files 
for sam_files in $mappedSAM/*.sam
do

  sortedBAM_file="$(basename $sam_files _mapped.sam)_sorted.bam"

  echo "Sample: $sortedBAM_file" | cat >> $sortedBAM/bam_sorting_stdout.txt
  echo "Sample: $sortedBAM_file" | cat >> $sortedBAM/bam_sorting_stderr.txt

  echo "Converting $sam_files to $sortedBAM_file ..."
  samtools view -b $sam_files > $sortedBAM/$sortedBAM_file 

  echo "Sorting $sortedBAM_file ..."
  samtools sort $sortedBAM/$sortedBAM_file -o $sortedBAM/$sortedBAM_file \
  1>>$sortedBAM/bam_sorting_stdout.txt \
  2>>$sortedBAM/bam_sorting_stderr.txt

  echo "
  " | cat >> $sortedBAM/bam_sorting_stdout.txt

  echo "
  " | cat >> $sortedBAM/bam_sorting_stderr.txt

done
###########################################################################
#Create summary statistics for aligned reads:

mkdir -p $QC3/RSeQC/bam_stat

> $QC3/RSeQC/bam_stat/std_err.txt

#Creates summary statistics for the bam files 
for bam_files in $sortedBAM/*.bam
do

  echo "Creating summary statistics for $bam_files..."

  #name each summary statistic in the file
  file_name="$(basename $bam_files .bam).txt"
  error_name="$(basename $bam_files)"

  echo "Sample: $error_name" | cat >> $QC3/RSeQC/bam_stat/std_err.txt

  #write out the summary statistcs for each file
  bam_stat.py -i $bam_files > $QC3/RSeQC/bam_stat/$file_name \
  2>> $QC3/RSeQC/bam_stat/std_err.txt

  echo "
  " | cat >> $QC3/RSeQC/bam_stat/std_err.txt

done

###########################################################################
#Quality control of read alignment:

mkdir -p $QC3/QualiMap

> $QC3/QualiMap/std_err.txt

for bam_files in $sortedBAM/*.bam
do
  echo "Running quality control with BamQC on $bam_files ..."

  #Naming the summary statistics 
  name="$(basename $bam_files .bam)"
  error_name="$(basename $bam_files)"

  echo "Sample: $error_name" | cat >> $QC3/QualiMap/std_err.txt

  qualimap bamqc \
  -outdir $QC3/QualiMap/$name \
  -bam $bam_files \
  -gff $gtf_ref \
  2>>$QC3/QualiMap/std_err.txt

  echo "
  " | cat >> $QC3/QualiMap/std_err.txt

done 

###########################################################################
#Quantify gene counts from reads aligned to reference genome

> $counts/counting_stderr.txt

for bam_files in $sortedBAM/*.bam
do
  echo "Getting gene counts for $bam_files"

  bam_count_file="$(basename $bam_files _sorted.bam).txt"

  echo "Sample: $bam_count_file" | cat >> $counts/counting_stderr.txt

  htseq-count \
  --format=bam \
  --stranded=yes \
  --type=transcript \
  --idattr=gene_id \
  --additional-attr=gene_name \
  $bam_files \
  $gtf_ref > \
  $counts/$bam_count_file \
  2>>$counts/counting_stderr.txt 

  echo "
  " | cat >> $counts/counting_stderr.txt

done

###########################################################################
#Complete quality control:
multiqc -outdir $parent_directory/$project_name/multiqc $parent_directory/$project_name