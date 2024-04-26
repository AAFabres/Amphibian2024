#! /bin/bash

##This scrip is a modified verision of a default script provided by Dr. Tonia Schwartz for the Functional Genomics Class of Spring 2024
######### Script for analyzing data without cleaning and trimming ############
## Objective: Download raw read from NCBI and run the full RNAseq pipeline without trimming
## This script is ment to run on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	suggested paramenters are below to submit this script.
## 		queue: medium
##		core: 8
##		time limit (HH:MM:SS): 20:00:00 
##		Memory: 16gb
##		
###############################################

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1
module load multiqc
#NO TRIMMING
#module load trimmomatic/0.39
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread

# Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
MyID=aubclsc0324

WD=/scratch/$MyID/RNAseqFrog
OP=/scratch/$MyID/RNAseqFrog/NoTrim
DD=$WD/RawData
RDQ=RawDataQuality
#adapters=TruSeq3-PE-2.fa		#specific adapter file for illumina in Trimmomatic
CD=$OP/CleanData
PCQ=PostCleanQualityNoTrim
REFD=$WD/XTropicalisRefGenome		#folder with mapped reference genome
MAPD=$OP/Map_HiSat2
COUNTSD=/$OP/Counts_StringTie
RESULTSD=/home/$MyID/PracticeRNAseq_FullNoTrim/Counts_H_S_2024
REF=UCB_Xtro_10.0			#name of the genome in NCBI



##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there.
mkdir -p ${WD}
mkdir -p ${DD}
mkdir -p ${OP}
## move to the Data Directory
cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
		## for fasterq-dump
		## or for fastq-dump 
			## 	-F 	Defline contains only original sequence name.
			## 	-I 	Append read id after spot id as 'accession.spot.readid' on defline.
		## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are from Bioproject PRJDB12187. An experiment on Buergeria otai, 3 samples heat stress (DRR316901, DRR316903 & DRR316904), 3 samples control (DRR316902, DRR316905 & DRR316906)
## https://www.ncbi.nlm.nih.gov/bioproject/PRJDB12187

vdb-config --interactive
fastq-dump -F --split-files DRR316901

fasterq-dump --split-files DRR316901
fasterq-dump --split-files DRR316902
fasterq-dump --split-files DRR316903
fasterq-dump --split-files DRR316904
fasterq-dump --split-files DRR316905
fasterq-dump --split-files DRR316906

##### Extra ####
## If you are downloading data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then you can compare the values in this file with the ones provided by the company.
md5sum ./* > md5sum.txt

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
mkdir -p ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

##### MultiQC to summarized the fastqc results!
cd ${WD}/${RDQ}
multiqc ${WD}/${RDQ}

#######  Tarball the directory containing the FASTQC and MultiQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.

#######
################****************** Step 2  Cleaning the data with Trimmomatic ###################################
#######
######Fabres
######Option1: No Trimming

## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
#mkdir -p ${CD}
#mkdir -p ${WD}/${PCQ}


## Move to Raw Data Directory
#cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
#ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
#scp /home/${MyID}/class_shared/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
#while read i
#do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
#	       java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
#					PE -threads 6 -phred33 \
#        	"$i"_1.fastq "$i"_2.fastq  \
#       	 ${CD}/"$i"_1_paired.fastq ${CD}/"$i"_1_unpaired.fastq  ${CD}/"$i"_2_paired.fastq ${CD}/"$i"_2_unpaired.fastq \
#       	 ILLUMINACLIP:TruSeq3-PE-2.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.

	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

#fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
#fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}


#done<list			# This is the end of the loop

################## Run MultiQC to summarize the fastqc results
### move to the directory with the cleaned data
#cd ${WD}/${PCQ}
#multiqc ${WD}/${PCQ}

########################  Now compress your results files from the Quality Assessment by FastQC 
#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
#tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*


#######
############################*********** Step 3 Mapping and Counting ************###########################
#######


## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
##Identify where your reference genome is
#scp /home/${MyID}/XTropicalis/${REF}.fna .
#scp /home/${MyID}/XTropicalis/${REF}.gff .

###  Identify exons and splice sites on the reference genome
gffread ${REF}.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 or .gff to .gft formate for HiSat2 to use.
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fasta XTropicalis_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
## Option1: NO TRIMMING so Data will be our RAWDATA in foler ${DD}
cd ${DD}  #### This is where our Raw paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: DRR316901_1_paired.fastq DRR316901_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: DRR316901

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${DD}/list  . 

## process the samples in the list, one by one using a while loop
while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "${REFD}"/XTropicales_index       \
    -1 "${DD}"/"$i"_1.fastq  -2 "${DD}"/"$i"_2.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: DRR316901.sam; Output: DRR316901_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  

	######################  Step 3b  Counting  		################

	mkdir -p "${COUNTSD}"/"$i"
	stringtie -p 6 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"_sorted.bam

done<list


#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt ${RESULTSD}


### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix that can be used for other programs like DESeq2. 
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv ${RESULTSD}
## move these results files to your personal computer for downstream statistical analyses in R studio.
