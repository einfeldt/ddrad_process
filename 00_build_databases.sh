#!/bin/bash

#############################################################################################################################
##Build databases for contaminant and parasite filtering using KRAKEN
#############################################################################################################################

##Download genomes of interest from NCBI FTP site, organise by database for filtering
##Before build, create folder for each database, directory structure: $DBNAME> library> $DBNAME> seqs.fa.gz
##Requires files: genomeaccession_keys.txt (corresponding to accession in filenames) genomeaccession_values.txt (corresponding to kraken:taxid|####) in ~/Desktop/Process_libraries/03_filters
	
##Build databases

##Builds databases for filters:
	##Algae (contaminants)
	##Trematodes (parasites)
	##Nematodes (parasites)
	##Myxozoa (parasites)
	##Microsporidia (parasites)
	##Cryptococci (parasites)
	##Oomycetes (parasites)
	##Apicomplexans (parasites)

#algae apicomplexans cryptococci microsporidians myxozoans nematodes oomycetes already processed

cd ~/Desktop/Process_libraries/03_filters
for DBNAME in trematodes; do
	echo "Building database for" "$DBNAME"
	##Download taxonomy and unzip genome files
	kraken-build --download-taxonomy --db $DBNAME

	for file in ~/Desktop/Process_libraries/03_filters/$DBNAME/library/$DBNAME/*.fna.gz; do
		filename=${file%%.fna.gz}
		gunzip "$filename.fna.gz"
		mv "$filename.fna" "$filename.fa"
	done

	##For every file in database, replace ">" with ">kraken:taxid|#### " value corresponding to accession key in filename
	##Fixes problem arising from Kraken's dependency on GIs in .fa files, but NCBI phasing out GIs
	declare -a accnkeys
	mapfile -t accnkeys < "~/Desktop/Process_libraries/03_filters/genomeaccession_keys.txt"
	declare -a accnvals
	mapfile -t accnvals < "~/Desktop/Process_libraries/03_filters/genomeaccession_vals.txt"	
	##Gets number of genomes in key and value arrays
	num_genomekeys=$(wc -l < "~/Desktop/Process_libraries/03_filters/genomeaccession_keys.txt")   	
	##Replacement loop for each sequence
	for seqfile in ~/Desktop/Process_libraries/03_filters/$DBNAME/library/$DBNAME/*.fa; do
		##Pulls accession from filename accession.*.fa
		accession=$(basename ${seqfile%%.*.fa})
		for i in $(seq 0 $num_genomekeys); do
			##Checks if $accession = accessionkey for position $i in array $accnkeys
			if [ "$accession" = "${accnkeys[$i]}" ]
				then
					##Replaces ">" with ">kraken:taxid|#### "where ">kraken:taxid|#### " corresponds to "$accession"
					sed -i "s/>/${accnvals[$i]}/g" $seqfile
			fi
		done
	done

	##Add sequences to database library
	cd ~/Desktop/Process_libraries/03_filters
	for seqs in ~/Desktop/Process_libraries/03_filters/$DBNAME/library/$DBNAME/*.fa; do
		kraken-build --add-to-library $seqs --db $DBNAME
	done		

	##Build database (takes ~36 hours for human genome reduced to 4GB) and clean initial files
	kraken-build --build --threads 8 --minimizer-len 13 --max-db-size 6 --work-on-disk --jellyfish-hash-size 640M --db $DBNAME
	kraken-build --db $DBNAME --clean

	echo "Finished building database for" "$DBNAME"
done

##After clean, directory structure: DBNAME> db.idx db.kdb taxonomy> names.dmp nodes.dmp

