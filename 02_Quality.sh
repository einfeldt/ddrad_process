####################################################################################################
#STEP 2: Quality filter and trim
####################################################################################################
##Trims Illumina adapters, cutsite&barcode read through, R1 cutsite readthrough (R2 readthrough too variable to trim), 
##Trims leading/trailing low quality bases, low quality bases via sliding window
##Removes short reads

cd ~/Desktop/Process_libraries
echo 'Started Step 2 - quality filter and trimming:'
date

#Trim sequences containing Illumina adapters with trimmomatic v0.33
#Order of trimming steps occurs in order listed in command line
#16 threads
#Lists of adapters must be in working directory 
#Trim Illumina adapters from R1 and R2 sequences using simple clipping sequence method (note: does not use palindrome method, as results in collapsing of R1 and R2 sequences)
#1#ILLUMINACLIP seedmismatches=2 palindromecliphreshold=30 simpleclipthreshold=7 (recommended 7-15, Score=0.6*Nmatchingbases-(Q/10)*NMismatches
#2#Trim cutsite/barcodes and cutsite/Yadapters sequentially using simple method. Lists: Nbpcutsitebarcodes.fa for N=8-13 (GGNC/barcode & AATT/Yadapter)
#2#length=13 simpleclipthreshold=6 (2 mismatch allowed)
#2#length=12 simpleclipthreshold=5 (2 mismatch allowed)
#2#length=11 simpleclipthreshold=5 (1 mismatch allowed)
#2#length=10 simpleclipthreshold=5 (1 mismatch allowed)
#2#length=9  simpleclipthreshold=4 (1 mismatch allowed)
#2#length=8  simpleclipthreshold=4 (1 mismatch allowed)
#3#Trim cutsite remnant (AATT) for R1 reads (can't trim GNCC from R2, as GNCC is not same as cutsite recognition GGNCC, results in trimming of actual sequence)
#3#length=4  (0 seed mismatches allowed) simpleclipthreshold=2
#4#Final quality trimming
#4#Removes leading low quality or N bases (below quality 30) (LEADING:30)
#4#Removes trailing low quality or N bases (below quality 30) (TRAILING:30)
#4#Scans the read with a 4-base wide sliding window, cutting when the average quality per base drops below 22 (SLIDINGWINDOW:4:22)
#4#Quality 22 based on Del Fabbro et al 2013 - Trimmomatic sliding window trim maximizes % mapping reads at Q=22
#4#Drops reads below 30 bases long (MINLEN:30)

#Regular trims on reads from 138 individuals, rename to be in format POP01.F.fq POP01.R.fq

for file in ./Individuals/?????_R1.fastq.gz; do 
	filename=$(basename $file)
	individual=${filename:0:5}
	TrimmomaticPE -threads 8 -phred33 ./Individuals/$individual"_R1.fastq.gz" ./Individuals/$individual"_R2.fastq.gz" ./Individuals/$individual".F.fastq.gz" ./Individuals/$individual"_R1unpaired.fastq.gz" ./Individuals/$individual".R.fastq.gz" ./Individuals/$individual"_R2unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE-ALL.fa:2:25:7 ILLUMINACLIP:13bpcutsitebarcodes.fa:2:30:6 ILLUMINACLIP:12bpcutsitebarcodes.fa:2:30:5 ILLUMINACLIP:11bpcutsitebarcodes.fa:2:30:5 ILLUMINACLIP:10bpcutsitebarcodes.fa:2:30:5 ILLUMINACLIP:9bpcutsitebarcodes.fa:1:30:4 ILLUMINACLIP:8bpcutsitebarcodes.fa:1:30:4 ILLUMINACLIP:R1_Cutsite.fa:0:30:2 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:28 MINLEN:30
done

#Cleanup
rm ./Individuals/*unpaired.fastq.gz
rm ./Individuals/*R?.fastq.gz

#Rename files to adhere to convention of POP_##.R1.fastq POP_##.R2.fastq (add underscore between POP and INDIVIDUAL)
cd ./Individuals
for file in * ; do
    mv ./"$file" "${file:0:3}_${file:3}"
done

####################################################################################################################
##Merge overlapping paired end reads with PEAR
####################################################################################################################
##Threads j=16

cd ~/Desktop/Process_libraries/02_stacks
mkdir ./01_merged
mkdir ./02_forward
mkdir ./03_reverse

for file in ~/Desktop/Process_libraries/Individuals/*.F.fastq.gz; do
	filename=$(basename $file)
	individual=${filename:0:6}
	pear -j 8 -f ../Individuals/"$individual".F.fastq.gz -r ../Individuals/"$individual".R.fastq.gz -o ./01_merged/"$individual"
	gzip ./01_merged/"$individual".assembled.fastq
	gzip ./01_merged/"$individual".unassembled.forward.fastq
	mv ./01_merged/"$individual".unassembled.forward.fastq.gz ./02_forward
	#PEAR output gives reverse complement of unassembled R sequences 
	#Reverse complement R sequences to maintain alignment of bases next to restriction site at start of sequence read for trimming to uniform size
	cat ./01_merged/$individual".unassembled.reverse.fastq" | fastx_reverse_complement -Q33 -o ./03_reverse/$individual".unassembled.reverse.fastq"
	gzip ./03_reverse/$individual".unassembled.reverse.fastq"
	rm ./01_merged/$individual".unassembled.reverse.fastq"
	rm ./01_merged/$individual".discarded.fastq"
	rm ../Individuals/$individual".F.fastq.gz"
	rm ../Individuals/$individual".R.fastq.gz"
done

##Remove .unassembled.forward., unassembled.reverse., .assembled from filename (this matters for catenating each branch together downstream)
cd ~/Desktop/Process_libraries/02_stacks/01_merged
for file in ~/Desktop/Process_libraries/02_stacks/01_merged/*.assembled.fastq.gz; do
    filename=$(basename "$file")
    mv "$file" "${filename%%.*}.fastq.gz"
done

cd ~/Desktop/Process_libraries/02_stacks/02_forward
for file in ~/Desktop/Process_libraries/02_stacks/02_forward/*.unassembled.forward.fastq.gz; do
    filename=$(basename "$file")
    mv "$file" "${filename%%.*}.fastq.gz"
done

cd ~/Desktop/Process_libraries/02_stacks/03_reverse
for file in ~/Desktop/Process_libraries/02_stacks/03_reverse/*.unassembled.reverse.fastq.gz; do
    filename=$(basename "$file")
    mv "$file" "${filename%%.*}.fastq.gz"
done

####################################################################################################################
##Trim files to standard length for stacks input
####################################################################################################################
##Trim all reads >85 bp, discard all reads <85bp
#85bp maximizes information=nseq*seqlength for each of F, R, merged
#Keeps ~65% of information (for trimmed/filtered 100000 subsample, originalinfo=85bp*(139432seqs*0.75halfoverlapping)=8888790bp; finalinfo=85bp*68319seqs=5807115bp)
#Performed before filtering to keep R reads in alignment with F and merged

mkdir ~/Desktop/Process_libraries/02_stacks/04_85bpreads
cd ~/Desktop/Process_libraries/02_stacks/04_85bpreads

for file in ~/Desktop/Process_libraries/02_stacks/01_merged/*.fastq.gz; do
	filename=$(basename "$file")
	individual=${filename:0:6}
	TrimmomaticSE -threads 8 -phred33 $file $individual"_merged".fastq.gz MINLEN:85 CROP:85
	rm ~/Desktop/Process_libraries/02_stacks/01_merged/$individual".fastq.gz"
done

for file in ~/Desktop/Process_libraries/02_stacks/02_forward/*.fastq.gz; do
	filename=$(basename "$file")
	individual=${filename:0:6}
	TrimmomaticSE -threads 8 -phred33 $file $individual"_forward".fastq.gz MINLEN:85 CROP:85
	rm ~/Desktop/Process_libraries/02_stacks/02_forward/$individual".fastq.gz"
done

for file in ~/Desktop/Process_libraries/02_stacks/03_reverse/*.fastq.gz; do
	filename=$(basename "$file")
	individual=${filename:0:6}
	TrimmomaticSE -threads 8 -phred33 $file $individual"_rev_reverse".fastq.gz MINLEN:85 CROP:85
	#Reverse complement R sequences after trimming to be in same orientation as merged and forward reads
	gunzip $individual"_rev_reverse".fastq.gz
	cat $individual"_rev_reverse".fastq | fastx_reverse_complement -Q33 -z -o $individual"_reverse.fastq.gz"
	rm $individual"_rev_reverse".fastq
	rm ~/Desktop/Process_libraries/02_stacks/03_reverse/$individual".fastq.gz"
done

##Catenate F, R, merged reads into single file for each individual 
mkdir ~/Desktop/Process_libraries/02_stacks/05_allreads

for file in ~/Desktop/Process_libraries/02_stacks/04_85bpreads/*_merged.fastq.gz; do
	filename=$(basename $file)
	individual=${filename:0:6}
	cat $individual"_merged.fastq.gz" $individual"_forward.fastq.gz" $individual"_reverse.fastq.gz" > ~/Desktop/Process_libraries/02_stacks/05_allreads/$individual".fastq.gz"
	rm ~/Desktop/Process_libraries/02_stacks/04_85bpreads/$individual"_merged.fastq.gz"
	rm ~/Desktop/Process_libraries/02_stacks/04_85bpreads/$individual"_forward.fastq.gz"
	rm ~/Desktop/Process_libraries/02_stacks/04_85bpreads/$individual"_reverse.fastq.gz"
done

echo 'Finished Step 2 - quality filter and trimming:'
date
