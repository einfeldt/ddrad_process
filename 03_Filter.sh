#!/bin/bash

#############################################################################################################################
##Filter sequences for parasites and contaminants using KRAKEN 
#############################################################################################################################

echo 'Started Step 3 - filtering contaminants and parasites:'
date

mkdir ~/Desktop/Process_libraries/04_output
mkdir ~/Desktop/Process_libraries/04_output/contaminants
mkdir ~/Desktop/Process_libraries/04_output/parasites
mkdir ~/Desktop/Process_libraries/04_output/sequences

##input: POP_ID.fastq.gz

## Filter databases
#contaminants
	#g general (minikraken)
	#h human (also filters 18s)
	#l algae
#parasites
	#a apicomplexans
	#c cryptococci
	#i microsporidians
	#m myxozoans
	#n nematodes
	#o oomycetes
	#t trematodes

## Copy input to working directory, unzip
for file in ~/Desktop/Process_libraries/02_stacks/05_allreads/*.fastq.gz; do
	cp $file ~/Desktop/Process_libraries/03_filters
done

for file in ~/Desktop/Process_libraries/03_filters/*.fastq.gz; do
	gunzip $file
done

## Filter each individual for general (contaminants)
cd ~/Desktop/Process_libraries/03_filters

dbname='general'
mkdir ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/???_??.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	db=minikraken_20141208
	kraken --db $db --preload --threads 8 --fastq-input "$individual".fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_g.fastq \
		> ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $db ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-g for human (contaminants)
dbname='human'
mkdir ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_g.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	db=kraken_human
	kraken --db $db --preload --threads 8 --fastq-input "$individual"_g.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_gh.fastq \
		> ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $db ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-gh for algae (contaminants)
dbname='algae'
mkdir ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_gh.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_gh.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghl.fastq \
		> ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/contaminants/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghl for apicomplexans (parasites)
dbname='apicomplexans'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghl.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghl.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghla.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghla for cryptococci (parasites)
dbname='cryptococci'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghla.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghla.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghlac.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done


## Filter each individual-ghlac for microsporidians (parasites)
dbname='microsporidians'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghlac.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghlac.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghlaci.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghlaci for myxozoans (parasites)
dbname='myxozoans'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghlaci.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghlaci.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghlacim.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghlacim for nematodes (parasites)
dbname='nematodes'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghlacim.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghlacim.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghlacimn.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghlacimn for oomycetes (parasites)
dbname='oomycetes'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghlacimn.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghlacimn.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/03_filters/"$individual"_ghlacimno.fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

## Filter each individual-ghlacimno for trematodes (parasites)
dbname='trematodes'
mkdir ~/Desktop/Process_libraries/04_output/parasites/"$dbname"
for file in ~/Desktop/Process_libraries/03_filters/*_ghlacimno.fastq; do
	filename=$(basename $file)
	individual=${filename:0:6}
	kraken --db $dbname --preload --threads 8 --fastq-input "$individual"_ghlacimno.fastq \
		--classified-out ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".fastq \
		--unclassified-out ~/Desktop/Process_libraries/04_output/sequences/"$individual".fastq \
		> ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken
	kraken-translate --db $dbname ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".kraken > ~/Desktop/Process_libraries/04_output/parasites/"$dbname"/"$individual"_"$dbname".labels
	rm $file
done

for file in ~/Desktop/Process_libraries/04_output/sequences/*.fastq; do
	gzip $file
done

echo 'Finished Step 3 - filtering contaminants and parasites:'
date
