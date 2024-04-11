#! /bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to run the full RNAseq pipeline
## For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##  You may need to increase these for bigger datasets
## 		queue: medium
##		core: 6
##		time limit (HH:MM:SS): 18:00:00 
##		Memory: 12gb
##		
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
#module load sra
module load fastqc/0.10.1
module load multiqc
module load trimgalore/0.6.7
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread
#module load gffcompare

#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x


##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsc0324         ## Example: MyID=aubtss


WD=/scratch/$MyID/RNAseqFrog            ## Example:/scratch/$MyID/PracticeRNAseq  
OP=/scratch/$MyID/RNAseqFrog/StrictGalore
DD=$WD/RawData
#RDQ=RawDataQuality
#adapters=TruSeq3-PE-2.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
				## In the future, for your data, you will likely need to edit this for other projects based on how your libraries 
				## were made to search for the correct adapters for your project
CD=$OP/CleanData            				## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files are located from the last script
PCQ=PostCleanQualityStrictGalore
REFD=$WD/XTropicalisRefGenome          ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed reference genome for the garter snake
MAPD=$OP/Map_HiSat2           			## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/$OP/Counts_StringTie       ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/PracticeRNAseq_FullStrictGalore/Counts_H_S_2024      ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S
REF=UCB_Xtro_10.0                  ## This is what the "easy name" will be for the genome reference



##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there.
#mkdir -p ${WD}
#mkdir -p ${DD}
## move to the Data Directory
#cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
		## for fasterq-dump
		## or for fastq-dump 
			## 	-F 	Defline contains only original sequence name.
			## 	-I 	Append read id after spot id as 'accession.spot.readid' on defline.
		## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are from Bioproject PRJNA437447. An experiment on Daphnia pulex, 5 samples on ad lib feed, 5 samples on caloric restriction diet
## https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=5206312
## For class only do the 6 that you are assigned, delete the other 4 from this list

#vdb-config --interactive
#fastq-dump -F --split-files DRR316901

#fasterq-dump DRR316901
#fasterq-dump DRR316902
#fasterq-dump DRR316903
#fasterq-dump DRR316904
#fasterq-dump DRR316905
#fasterq-dump DRR316906



##### Extra ####
## If you are downloading data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then you can compare the values in this file with the ones provided by the company.
#md5sum ./* > md5sum.txt

##### Extra ####
## If you data comes with multiple R1 and R2 files per individual. You can contatenate them together using "cat" before running FASTQC
## see examples below for one file. You will probably want to use a loop to process through all the files.
#cat SRR6819014*_R1_*.fastq.gz > SRR6819014_All_R1.fastq.gz
#cat SRR6819014*_R2_*.fastq.gz > SRR6819014_All_R2.fastq.gz


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
#mkdir ${WD}/${RDQ}
#fastqc *.fastq --outdir=${WD}/${RDQ}

##### MultiQC to summarized the fastqc results!
#cd ${WD}/${RDQ}
#multiqc ${WD}/${RDQ}

#######  Tarball the directory containing the FASTQC and MultiQC results so we can easily bring it back to our computer to evaluate.
#tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.


#######
################****************** Step 2  Cleaning the data with Trim Galore! ###################################
#######
#### Option strict


## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir -p ${CD}
mkdir -p ${OP}/${PCQ}


## Move to Raw Data Directory
cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
#cp /home/${MyID}/class_shared/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
	      
trim_galore --paired --q 30 --illumina --three_prime_clip_R2 10 --length 36 "$i"_1.fastq "$i"_2.fastq -o ${CD}/

# --paired is to tell the program that your using paired data
# --q this is the Qthred score that will filter data below Qphred scroe 20, which is the default non-strict model
# --ilumina only filter for 'AGATCGGAAGAGC' this adapter	#Our data is illumina
# not putting an adapter option it give me that the adapter is 'CTGTCTCTTATA'
# -o where the output files should go

        
                ## Trim read for quality when quality drops below Q20
                ## output file 1: DRR316901_1_trimmed.fq
		            ## output file 2: DRR316901_2_trimmed.fq 
		 
	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

fastqc ${CD}/"$i"_1_trimmed.fq --outdir=${OP}/${PCQ}
fastqc ${CD}/"$i"_2_trimmed.fq --outdir=${OP}/${PCQ}


done<list			# This is the end of the loop

################## Run MultiQC to summarize the fastqc results
### move to the directory with the cleaned data
cd ${OP}/${PCQ}
multiqc ${OP}/${PCQ}

########################  Now compress your results files from the Quality Assessment by FastQC 
#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${PCQ}.tar.gz ${OP}/${PCQ}/*


#######
############################*********** Step 3 Mapping and Counting ************###########################
#######


## Make the directories and all subdirectories defined by the variables above
#mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
#cd $REFD
### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
#scp /home/${MyID}/XTropicalis/${REF}.fna .
#scp /home/${MyID}/XTropicalis/${REF}.gff .

###  Identify exons and splice sites on the reference genome
#gffread ${REF}.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
#hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
#hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
#hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna XTropicales_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd ${CD}  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${CD}/list  . 

## process the samples in the list, one by one using a while loop
while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x ${REFD}/XTropicales_index       \
    -1 ${CD}/"$i"_1_trimmed.fq  -2 ${CD}/"$i"_2_trimmed.fq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
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

	mkdir "${COUNTSD}"/"$i"
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
