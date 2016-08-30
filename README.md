# pollen-DNA-metabarcoding
Evaluating the Relative Performance of DNA Metabarcoding Sequence Classifiers
Rodney T. Richardson, Johan Bengtsson-Palme and Reed M. Johnson

## Table of Contents
1. Introduction
2. Retrieve data and filter by Ohio plant species list
2.1 Obtain data from NCBI
2.2 Preliminary formatting and quality curation
2.3 Filter by geography
3. Randomly sample testing sequences from training sequences
4. Format reference sequences and train classifiers
4.1 Format training sequences
4.2 Train RDP classifier
4.3 Train UTAX classifier
4.4 Re-format references for RTAX classifier
5. Format testing sequences
5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages
5.2 Replace Gi numbers in fasta with lineages
5.3 Format sequences to facilitate later analyses
6. Classify testing sequences
6.1 Classify sequences using RDP
6.2 Classify sequences using RTAX
6.3 Classify sequences using UTAX
7. Format classifier output and perform accuracy and sensitivity analysis
7.1 Format RDP output
7.2 Format RTAX output
7.3 Format UTAX output
7.4 Perform accuracy and sensitivity analysis
8. Calculate DB coverage
9. Isolate sequences with no genus level reference representation

## 1. Introduction
This document contains the commands used for evaluation of classifier performance. These commands are provided such that readers can work through our analyses independently and apply the approach to their own research endeavors. It should be noted however that the syntax, commands and software used here may not be entirely transferrable for future applications given differences in computational architecture, software updates, etc. Further, these commands are given without guidance in terms of directory organization, which we leave at the discretion of the reader. Additionally, due to the ease of transferring commands from one analysis to another, we do not provide commands for all the analyses performed in the paper. Lastly, while these commands were performed on files representing five loci, for ease of organization we generally give only the commands used for the rbcL locus as an example.

2. Retrieve data and filter by Ohio plant species list

## 2.1 Obtain data from NCBI
### search the following in NCBI
•	ITS2 AND vascular plants [ORGN]
•	rbcL AND vascular plants [ORGN]
•	matK AND vascular plants [ORGN]
•	trnH AND psbA AND vascular plants [ORGN]
•	trnL AND vascular plants [ORGN]
2.2 Preliminary formatting and quality curation
### remove newline instances from NCBI formatted fasta files
`awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNextLine.fasta` 
### filter by minimum length threshold
`awk '!/^>/ { next } { getline seq } length(seq) >= 400 { print $0 "\n" seq }' GenBank_rbcL_030416_rmNextLine.fasta > GenBank_rbcL_030416_400.fasta`
### filter by maximum length threshold
`awk '!/^>/ { next } { getline seq } length(seq) <= 1500 { print $0 "\n" seq }' GenBank_rbcL_030416_500.fasta > GenBank_rbcL_030416_400_1500bp.fasta`
### Length trimming thresholds used
•	200 – 800 bp for ITS2
•	500 – 1000 bp for matK
•	400 – 1500 bp for rbcL
•	150 – 1100 bp for trnL
•	200 – 800 bp for trnH
### remove sequences with more than two consecutive uncalled base pairs
`sed -i '/NN/I,+1 d' GenBank_rbcL_030416_400_1500bp.fasta`
## 2.3 Filter by geography 
### get Gi numbers from NCBI entries for species found in Ohio (OhioSpecies.txt file is provided in Appendix S2). It should be noted that this approach to filtering the sequences is conservative as only one sequence representative per species is retained. Furthermore, one should be aware of that there could be naming inconsistencies between the list of species used for filtering and, e.g., the NCBI naming convention. Thus, this step should only be carried out if relevant for the questions investigated.
`grep -f OhioSpecies.txt GenBank_rbcL_030416_400_1500bp.fasta | awk '!seen[$2$3]++' | perl -pe 's/^>gi\|(\d+)\|.*/$1/' > OH_GB_030416_rbcL.txt`
### make blast database with reference sequences
`~/PATH/TO/ncbi-blast-2.2.29+/bin/makeblastdb -in GenBank_rbcL_030416_400_1500bp.fasta -dbtype nucl -parse_seqids`
### use blastdbcmd to retrieve sequences matching the Gi numbers of plants known to occur in Ohio
`~/PATH/TO/ncbi-blast-2.2.29+/bin/blastdbcmd -db GenBank_rbcL_030416_400_1500bp.fasta -entry_batch OH_GB_030416_rbcL.txt -out OH_Spcs_rbcL.fa`
### remove everything from sequence headers except “>Gi” to meet formatting requirements of downstream applications
`cat OH_Spcs_rbcL.fa | perl -pe 's/^>gi\|(\d+)\|.*/>$1/g' > OH_Spcs_rbcL_400_1500bp_CLEAN.fasta`
### remove next line instances which were added during blast database usage
`awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNL.fasta`
## 3. Randomly sample testing sequences from training sequences
### sample the headers from 15% of the sequences to serve as the testing set
`awk 'BEGIN {srand()} !/^$/ { if (rand() <= .15) print $0}' < OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TestSet.txt`
### use reverse grep to get headers of sequences which were not sampled for the testing set
`grep -v -f OH_GB_030416_rbcL_TestSet.txt OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TrainSet.txt`
### make a fasta file with the header and sequence for the testing set entries
`grep -A 1 -f OH_GB_030416_rbcL_TestSet.txt OH_Spcs_rbcL_400_1500bp_CLEAN_2.fasta > OH_rbcL_TestSet.fasta`
### make a fasta file with the header and sequence for the training set entries
`grep -A 1 -f OH_GB_030416_rbcL_TrainSet.txt OH_Spcs_rbcL_200_800bp_CLEAN_2.fasta > OH_rbcL_TrainSet.fasta`
### delete grep “--” artifact 
`sed -i '/--/d' OH_rbcL_*Set.fasta`
## 4. Format reference sequences and train classifiers
## 4.1 Format training sequences
### get Gi numbers from training sequences
`grep '>' OH_rbcL_TrainSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_ohio.gis`
### use the following two Perl scripts provided in Sickel at al. (2015) and the NCBI Taxonomy module to produce required taxonomy and fasta files for later use with RDP and UTAX. See https://github.com/iimog/meta-barcoding-dual-indexing for the Perl scripts and instructions on installing the taxonomy module
`perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_ohio.gis --out rbcL_ohio.tax --species rbcL_ohio.species.taxids --genus rbcL_ohio.genus.taxids`
`perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/tax2rdp_utax.pl rbcL_ohio.tax OH_rbcL_TrainSet.fasta OH_rbcL_trainRDP`
### remove duplicate sequences
`java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-dupseq --infile OH_rbcL_trainRDP.rdp.fa --outfile OH_rbcL_trainRDP.rmDS.rdp.fa --duplicates --min_seq_length 50`
### remove partial sequences
`java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-partialseq OH_rbcL_trainRDP.rdp.fa OH_rbcL_trainRDP.rmDS.rdp.fa OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --alignment-mode overlap --min_gaps 50 --knn 20`
### remove entries with incomplete taxonomy data
`sed -i '/undef/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa`
`sed -i '/ sp/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa`
`sed -i '/incertae/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa`
`sed -i '/incerti/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa`
`sed -i '/unknown/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa`
## 4.2 Train RDP classifier
`java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar train --out_dir rbcL_trained --seq OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --tax_file OH_rbcL_trainRDP.rdp.tax`
### move rbcL.properties, produced during training, from data directory to rdp_trained directory
`mv ./rbcL.properties ./rbcL_trained/rbcL.properties`
## 4.3 Train UTAX classifier
### use training sequences curated for RDP, OH_rbcL_trainRDP.FinalTrainSet.rdp.fa from previous page, as downstream classifier training data for UTAX and as reference data for RTAX
### reformat fasta headers from RDP format to UTAX format
`cat OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)\tRoot;/>$1;tax=/g' | perl -pe 's/(.)__/$1:/g' | perl -pe 's/_\d+;/,/g' | perl -pe 's/,$/;/g' > rbcL.utax.fa`
### make UTAX training database
`~/PATH/TO/Usearch/usearch_v8.1 -makeudb_utax rbcL.utax.fa -output rbcL.udb`
## 4.4 Re-format references for RTAX classifier
### reformat the taxonomic lineage headers of the RDP training fasta file to make a taxonomy file which meets the requirements of the RTAX classifier 
`grep '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)(\t)Root;/$1$2/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__/ /g' | perl -pe 's/_\d+//g' > rbcL_training.rtax.tax`
### reformat RDP training fasta file for use with RTAX
`grep -A 1 '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+).*/>$1/g' > rbcL_training.rtax.fasta`
## 5. Format testing sequences
## 5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages
### get file of Gi numbers for each test sequence entry
`grep '>' OH_rbcL_TestSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_test_ohio.gis`
### use Gi numbers along with the Perl script from Sickel et al. (2015) and the NCBI taxonomy module to get the taxonomic lineage of each entry
`perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_test_ohio.gis --out rbcL_test_ohio.tax --species rbcL_test_ohio.species.taxids --genus rbcL_test_ohio.genus.taxids`
## 5.2 Replace Gi numbers in fasta with lineages
### use the following bash script to replace fasta Gi headers with the taxonomic lineage for each entry
`#!/bin/bash`
`cat rbcL_test_ohio.gis | while read Gi`
`do`
`	grep -w "$Gi" rbcL_test_ohio.tax | perl -pe 's/(\d+)\tRoot;(.*)/>$1;$2/g' >> rbcL_test_RFMT.fasta`
`	grep -A 1 -w "$Gi" OH_rbcL_TestSet.fasta | grep -v '>' >> rbcL_test_RFMT.fasta`
`done`
## 5.3 Format sequences to facilitate later analyses 
### remove taxonomically incomplete entries
`sed -i '/undef/,+1 d' rbcL_test_RFMT.fasta`
`sed -i '/ sp/,+1 d' rbcL_test_RFMT.fasta`
`sed -i '/incertae/,+1 d' rbcL_test_RFMT.fasta`
`sed -i '/incerti/,+1 d' rbcL_test_RFMT.fasta`
`sed -i '/unknown/,+1 d' rbcL_test_RFMT.fasta`
## 6. Classify testing sequences
## 6.1 Classify sequences using RDP 
### replace tabs in header with “-“ so classifier doesn't delete the characters after the space.
`grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/ 	/-/g' > rbcL_test_rdpFMT.fasta`
`java -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar classify -t ~/PATH/TO/rbcL.properties --outputFile rbcL.rdp.output rbcL_test_rdpFMT.fasta`
## 6.2 Classify sequences using RTAX
### replace “;” with “–“ in header so that whole lineage is retained in RTAX output
`grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/;/-/g' > rbcL_test_rtaxFMT.fasta`
### classify sequences
`~/PATH/TO/rtax -r rbcL_training.rtax.fasta -t rbcL_training.rtax.tax -a rbcL_rtaxFMT.fasta -o rbcL.rtax.output`
## 6.3 Classify sequences using UTAX 
`~/local/src/Usearch/usearch_v8.1 -utax rbcL_test_RFMT.fasta -db ~/PATH/TO/rbcL.udb -rdpout rbcL.utax.output -strand plus`
## 7. Format classifier output and perform accuracy and sensitivity analysis
## 7.1 Format RDP output
### replace '\t-\t' with '\t' and replace tabs with commas 
`cat rbcL.rdp.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_clean.rdp.output`
### replace RDP assignments with less than 90% bootstrap support with “NULL” using the following bash script (this facilitates later analysis in R)
`#!/bin/bash`
`### this for loop can be applied to all five loci if all testing sequences have been classified`
`for i in rbcL`
`do`
`cat ${i}_clean.rdp.output | while read LINE`
`do` 
`ACTUAL_LINEAGE_FINAL=$(echo $LINE | cut -d',' -f1)`
`KINGDOM=$(echo $LINE | cut -d',' -f5,6,7)`
`PHYLUM=$(echo $LINE | cut -d',' -f8,9,10)`
`CLASS=$(echo $LINE | cut -d',' -f11,12,13)`
`ORDER=$(echo $LINE | cut -d',' -f14,15,16)`
`FAMILY=$(echo $LINE | cut -d',' -f17,18,19)`
`GENUS=$(echo $LINE | cut -d',' -f20,21,22)`
`SPECIES=$(echo $LINE | cut -d',' -f23,24,25)`

`K_VAL=$(echo $KINGDOM | cut -d',' -f3)`
`P_VAL=$(echo $PHYLUM | cut -d',' -f3)`
`C_VAL=$(echo $CLASS | cut -d',' -f3)`
`O_VAL=$(echo $ORDER | cut -d',' -f3)`
`F_VAL=$(echo $FAMILY | cut -d',' -f3)`
`G_VAL=$(echo $GENUS | cut -d',' -f3)`
`S_VAL=$(echo $SPECIES | cut -d',' -f3)`
`echo {$K,$P,$C,$O,$G,$S}_VAL >> QC_CHECK_${i}`
`if [ $K_VAL \< 0.90 ]`
`then`
`K_FINAL=NULL`
`else`
`K_FINAL=$(echo $LINE | cut -d',' -f5,6,7)`
`fi`
`if [ $P_VAL \< 0.90 ]`
`then`
`P_FINAL=NULL`
`else`
`P_FINAL=$(echo $LINE | cut -d',' -f8,9,10)`
`firbcL_0.5DB_OUT`
`if [ $C_VAL \< 0.90 ]`
`then`
`C_FINAL=NULL`
`else`
`C_FINAL=$(echo $LINE | cut -d',' -f11,12,13)`
`fi`
`if [ $O_VAL \< 0.90 ]`
`thenmatK_Test_OUTPUT`
`O_FINAL=NULL`
`else`
`O_FINAL=$(echo $LINE | cut -d',' -f14,15,16)`
`fi`
`if [ $F_VAL \< 0.90 ]`
`then`
`F_FINAL=NULL`
`else`
`F_FINAL=$(echo $LINE | cut -d',' -f17,18,19)`
`fi`
`if [ $G_VAL \< 0.90 ]`
`then`
`G_FINAL=NULL`
`else`
`G_FINAL=$(echo $LINE | cut -d',' -f20,21,22)`
`fi`
`if [ $S_VAL \< 0.90 ]`
`then`
`S_FINAL=NULL`
`else`
`S_FINAL=$(echo $LINE | cut -d',' -f23,24,25)`
`fi`
`echo {$ACTUAL_LINEAGE,$K,$P,$C,$O,$F,$G,$S}_FINAL >> ${i}_RDP_OUTPUT`
`done`
`done`

### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RDP assignments (the predicted identity)
`cat rbcL_RDP_OUTPUT | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi' > rbcL_RDP_FINAL`

## 7.2 Format RTAX output
### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RTAX assignments (the predicted identity)
`cat rbcL.rtax.output | perl -pe 's/(\d+)-(.__)/$1,$2/g' | perl -pe 's/\t/,/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' | perl -pe 's/, /,/g' | perl -pe 's/-//g' > rbcL_RTAX_FINAL`

### import data into excel and replace empty cells with “NULL”, save as .csv and import to R
## 7.3 Format UTAX output
### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the `UTAX assignments (the predicted identity)
cat rbcL.utax.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_UTAX_output`

### use bash shell from RDP output reformatting (section 7.1) to remove assignments below 0.90 confidence and replace empty cells with “NULL”

### further reformatting for later analysis
`cat rbcL_UTAX_OUTPUT | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g' > rbcL_UTAX_FINAL`
## 7.4 Perform accuracy and sensitivity analysis
### isolate classified sequences to determine sensitivity in R (the example given is for RTAX at the genus level but this would have to be changed for other ranks and classifiers)
`genus_assigned <- as.data.frame(rbcL_RTAX_FINAL[which(rbcL_RTAX_FINAL $V23 != 'NULL'),c(7,23)])`
### compare testing sequence taxonomy with classifier predicted taxonomy to isolate mis-classifications in R
`genus_misIDs <- as.data.frame(rbcL_RTAX_FINAL[which(as.character(rbcL_RTAX_FINAL$V7) != as.character(rbcL_RTAX_FINAL$V15) & as.character(rbcL_RTAX_FINAL$V15) != 'NULL'),c(7,15)])`
### visually analyze the resulting data frames to ensure that the correct columns were compared and that there are no obviously confounding artifacts in the data
## 8. Calculate DB coverage
### Obtain lists of genera from Ohio and both the training sets and the testing sets
`cat OhioSpecies.txt | cut -d' ' -f1 | sort | uniq > rbcL_Ohio_genera`
`cat rbcL_test_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi’ | sort | uniq  > rbcL_test_genera`
`cat rbcL_train_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi' | sort | uniq > rbcL_train_genera`
### import these to excel and save as .csv
### use %in% function in R to find overlap of genera lists
`rbcL_test_train <- as.data.frame(rbcL_test[which(rbcL_test[,1] %in% rbcL_train[,1]),1])`

## 9. Isolate sequences with no genus level reference representation
### find genera present in the testing set but not in the training set
`grep -v -f rbcL_train_genera rbcL_test_genera > rbcL_in_test_not_train`
