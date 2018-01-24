#!/bin/bash

##Before running, check species barcode ID lists

##Start of outer loop that runs all steps below for each library in ~/Desktop/Process_libraries/Input
for libraryfile in ~/Desktop/Process_libraries/Input/*R1.fastq.gz; do
	libraryfilename=$(basename $libraryfile)
	LIBRARY=${libraryfilename%_R1.fastq.gz}
	

	####################################################################################################
	#STEP 1: Demultiplex paired end reads
	####################################################################################################
	##Filters all reads by R2 barcodes, demultiplexes by F barcode, removes barcodes and cutsite remnants
	##All inputs must be renamed to format Libraryname_R1.fastq.gz Libraryname_R2.fastq.gz before running

	cd ~/Desktop/Process_libraries/01_Demultiplex
	echo 'Started demultiplexing by R2 Y-adapters:'
	date

	##Demultiplex by all R2 Y-adapters with 0 mismatches (Input1=R2, Input2=R1)
	mkdir Yadapters_demultiplexed
	fastq-multx ./R2_barcode_sequences/All_Yadapters.txt ../Input/$LIBRARY"_R2.fastq.gz" ../Input/$LIBRARY"_R1.fastq.gz" -o ./Yadapters_demultiplexed/r2.%.fastq.gz -o ./Yadapters_demultiplexed/r1.%.fastq.gz -b -m 0 -d 1

	##Move sequences that do not match Y-adapter barcodes 
	mkdir Yadapters_unmatched
	mv ./Yadapters_demultiplexed/r?.unmatched.fastq.gz ./Yadapters_unmatched

	##Trim 4bp cut-site remnant from demultiplexed R2 sequences & first 2 bases of reads (due to non-random per-base sequence content)
	for file in ./Yadapters_demultiplexed/r2.Yadapter*; do 
		filename=$(basename $file)
		fileroot=${filename:0:14}
		gunzip $file	
		cat ~/Desktop/Process_libraries/01_Demultiplex/Yadapters_demultiplexed/$fileroot.fastq | fastx_trimmer -f 6 -Q33 -z -o ~/Desktop/Process_libraries/01_Demultiplex/Yadapters_demultiplexed/$fileroot.fastq.gz	
		rm ~/Desktop/Process_libraries/01_Demultiplex/Yadapters_demultiplexed/$fileroot.fastq
	done

	##Put all demultiplexed reads into single file (for each of R1 and R2)
	cat ./Yadapters_demultiplexed/r1.Yadapter* >> ./Yadapters_demultiplexed/R1_Yadapters_demultiplexed.fastq.gz
	cat ./Yadapters_demultiplexed/r2.Yadapter* >> ./Yadapters_demultiplexed/R2_Yadapters_demultiplexed.fastq.gz

	##Echo the date and time at finish
	echo 'Finished demultiplexing R2 Y-adapters and trimming 4bp cut-site remnant. Started demultiplexing R1:'
	date

	##Demultiplex by R1 barcodes
	mkdir ../Individuals_unnamed
	fastq-multx ./R1_barcode_sequences/All_R1barcodes.txt ./Yadapters_demultiplexed/R1_Yadapters_demultiplexed.fastq.gz ./Yadapters_demultiplexed/R2_Yadapters_demultiplexed.fastq.gz -o ../Individuals_unnamed/r1.%.fastq.gz -o ../Individuals_unnamed/r2.%.fastq.gz -b -m 0 -d 1
	rm ../Individuals_unnamed/r?.unmatched.fastq.gz

	##Trim 4bp cut-site remnant from demultiplexed R1 sequences & first base of reads (due to non-random per-base sequence content)
	cd ~/Desktop/Process_libraries/Individuals_unnamed
	for file in ./r1.barcode_*; do 
		filename=$(basename $file)
		fileroot=${filename:0:13}
		gunzip $file	
		cat ./$fileroot.fastq | fastx_trimmer -f 6 -Q33 -z -o ./$fileroot.fastq.gz
		rm ./$fileroot.fastq
	done

	##Cleanup
	rm -r ~/Desktop/Process_libraries/01_Demultiplex/Yadapters_demultiplexed
	rm -r ~/Desktop/Process_libraries/01_Demultiplex/Yadapters_unmatched

	##Rename sequences to match individual ID
	COUNTER=0
	for file in `ls r?.barcode_??.fastq.gz | sort`; do
	    COUNTER=$(($COUNTER+1))
		mv $file `sed -n ${COUNTER}p ../01_Demultiplex/$LIBRARY"_IDs.txt"`
	done

	##Remove Hediste sequences if they exist
	cd ~/Desktop/Process_libraries/Individuals_unnamed
	find . -name '*_C1.fastq.gz' -delete
	find . -name '*_C2.fastq.gz' -delete
	
	##Remove sequences with ghost barcodes if they exist 
	##Demultiplexing with barcodes not used in library: false positive match rate of 0.0015% for 9bp barcodes, 0.008% for 4bp barcodes 
	##False matches increase with decreasing barcode complexity
	##Therefore, false matches likely due to errors in sequencing and/or custom oligo impurities; less likely due to contamination (false matches not in use at time of library prep)
	cd ~/Desktop/Process_libraries/Individuals_unnamed
	find . -name '*_B1.fastq.gz' -delete
	find . -name '*_B2.fastq.gz' -delete

	##Move processed sequences from /Individuals_unnamed to /Individuals
	mkdir ~/Desktop/Process_libraries/Individuals
	mv ~/Desktop/Process_libraries/Individuals_unnamed/*.fastq.gz ~/Desktop/Process_libraries/Individuals
	#rm -r ~/Desktop/Process_libraries/Individuals_unnamed

	echo 'Finished demultiplexing:'
	date

done
