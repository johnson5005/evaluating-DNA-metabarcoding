Appendix S1</span>

Evaluating the Relative Performance of DNA Metabarcoding Sequence Classifiers</span>

Rodney T. Richardson, Johan Bengtsson-Palme and Reed M. Johnson</span>

Table of Contents</span>

[1\. Introduction](#h.gjdgxs)</span>[        ](#h.gjdgxs)</span>

[2\. Retrieve data and filter by Ohio plant species list](#h.30j0zll)</span>[        ](#h.30j0zll)</span>

[2.1 Obtain data from NCBI](#h.1fob9te)</span>[        ](#h.1fob9te)</span>

[2.2 Preliminary formatting and quality curation](#h.3znysh7)</span>[        ](#h.3znysh7)</span>

[2.3 Filter by geography](#h.2et92p0)</span>[        ](#h.2et92p0)</span>

[3\. Randomly sample testing sequences from training sequences](#h.tyjcwt)</span>[        ](#h.tyjcwt)</span>

[4\. Format reference sequences and train classifiers](#h.3dy6vkm)</span>[        ](#h.3dy6vkm)</span>

[4.1 Format training sequences](#h.1t3h5sf)</span>[        ](#h.1t3h5sf)</span>

[4.2 Train RDP classifier](#h.4d34og8)</span>[        ](#h.4d34og8)</span>

[4.3 Train UTAX classifier](#h.2s8eyo1)</span>[        ](#h.2s8eyo1)</span>

[4.4 Re-format references for RTAX classifier](#h.17dp8vu)</span>[        ](#h.17dp8vu)</span>

[5\. Format testing sequences](#h.3rdcrjn)</span>[        ](#h.3rdcrjn)</span>

[5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages](#h.26in1rg)</span>[        ](#h.26in1rg)</span>

[5.2 Replace Gi numbers in fasta with lineages](#h.lnxbz9)</span>[        ](#h.lnxbz9)</span>

[5.3 Format sequences to facilitate later analyses](#h.35nkun2)</span>[        ](#h.35nkun2)</span>

[6\. Classify testing sequences](#h.1ksv4uv)</span>[        ](#h.1ksv4uv)</span>

[6.1 Classify sequences using RDP](#h.44sinio)</span>[        ](#h.44sinio)</span>

[6.2 Classify sequences using RTAX](#h.z337ya)</span>[        ](#h.z337ya)</span>

[6.3 Classify sequences using UTAX](#h.3j2qqm3)</span>[        ](#h.3j2qqm3)</span>

[7\. Format classifier output and perform accuracy and sensitivity analysis](#h.1y810tw)</span>[        ](#h.1y810tw)</span>

[7.1 Format RDP output](#h.4i7ojhp)</span>[        ](#h.4i7ojhp)</span>

[7.2 Format RTAX output](#h.1ci93xb)</span>[        ](#h.1ci93xb)</span>

[7.3 Format UTAX output](#h.2bn6wsx)</span>[        ](#h.2bn6wsx)</span>

[7.4 Perform accuracy and sensitivity analysis](#h.3as4poj)</span>[        ](#h.3as4poj)</span>

[8\. Calculate DB coverage](#h.1pxezwc)</span>[        ](#h.1pxezwc)</span>

[9\. Isolate sequences with no genus level reference representation](#h.49x2ik5)</span>[        ](#h.49x2ik5)</span>

[](#_Toc453838845)</span>

[](#_Toc453838845)</span>

# 1\. Introduction</span>

This document contains the commands used for evaluation of classifier performance. These commands are provided such that readers can work through our analyses independently and apply the approach to their own research endeavors. It should be noted however that the syntax, commands and software used here may not be entirely transferrable for future applications given differences in computational architecture, software updates, etc. Further, these commands are given without guidance in terms of directory organization, which we leave at the discretion of the reader. Additionally, due to the ease of transferring commands from one analysis to another, we do not provide commands for all the analyses performed in the paper. Lastly, while these commands were performed on files representing five loci, for ease of organization we generally give only the commands used for the</span> rbcL</span> locus as an example.</span>

# 2\. Retrieve data and filter by Ohio plant species list</span>

## 2.1 Obtain data from NCBI</span>

 search the following in NCBI</span>

*   ITS2 AND vascular plants [ORGN]</span>
*   rbcL AND vascular plants [ORGN]</span>
*   matK AND vascular plants [ORGN]</span>
*   trnH AND psbA AND vascular plants [ORGN]</span>
*   trnL AND vascular plants [ORGN]</span>

## 2.2 Preliminary formatting and quality curation</span>

 remove newline instances from NCBI formatted fasta files</span>

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNextLine.fasta</span>

 filter by minimum length threshold</span>

awk '!/^>/ { next } { getline seq } length(seq) >= 400 { print $0 "\n" seq }' GenBank_rbcL_030416_rmNextLine.fasta > GenBank_rbcL_030416_400.fasta</span>

 filter by maximum length threshold</span>

awk '!/^>/ { next } { getline seq } length(seq) <= 1500 { print $0 "\n" seq }' GenBank_rbcL_030416_500.fasta > GenBank_rbcL_030416_400_1500bp.fasta</span>

 Length trimming thresholds used</span>

*   200 – 800 bp for ITS2</span>
*   500 – 1000 bp for matK</span>
*   400 – 1500 bp for rbcL</span>
*   150 – 1100 bp for trnL</span>
*   200 – 800 bp for trnH</span>

 remove sequences with more than two consecutive uncalled base pairs</span>

sed -i '/NN/I,+1 d' GenBank_rbcL_030416_400_1500bp.fasta</span>

## 2.3 Filter by geography</span>

 get Gi numbers from NCBI entries for species found in Ohio (OhioSpecies.txt file is provided in Appendix S2). It should be noted that this approach to filtering the sequences is conservative as only one sequence representative per species is retained. Furthermore, one should be aware of that there could be naming inconsistencies between the list of species used for filtering and, e.g., the NCBI naming convention. Thus, this step should only be carried out if relevant for the questions investigated.</span>

grep -f OhioSpecies.txt GenBank_rbcL_030416_400_1500bp.fasta | awk '!seen[$2$3]++' | perl -pe 's/^>gi\|(\d+)\|.*/$1/' > OH_GB_030416_rbcL.txt</span>

 make blast database with reference sequences</span>

~/PATH/TO/ncbi-blast-2.2.29+/bin/makeblastdb -in GenBank_rbcL_030416_400_1500bp.fasta -dbtype nucl -parse_seqids</span>

 use blastdbcmd to retrieve sequences matching the Gi numbers of plants known to occur in Ohio</span>

~/PATH/TO/ncbi-blast-2.2.29+/bin/blastdbcmd -db GenBank_rbcL_030416_400_1500bp.fasta -entry_batch OH_GB_030416_rbcL.txt -out OH_Spcs_rbcL.fa</span>

 remove everything from sequence headers except “>Gi” to meet formatting requirements of downstream applications</span>

cat OH_Spcs_rbcL.fa | perl -pe 's/^>gi\|(\d+)\|.*/>$1/g' > OH_Spcs_rbcL_400_1500bp_CLEAN.fasta</span>

 remove next line instances which were added during blast database usage</span>

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNL.fasta</span>

# 3\. Randomly sample testing sequences from training sequences</span>

 sample the headers from 15% of the sequences to serve as the testing set</span>

awk 'BEGIN {srand()} !/^$/ { if (rand() <= .15) print $0}' < OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TestSet.txt</span>

 use reverse grep to get headers of sequences which were not sampled for the testing set</span>

grep -v -f OH_GB_030416_rbcL_TestSet.txt OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TrainSet.txt</span>

</span>

</span>

 make a fasta file with the header and sequence for the testing set entries</span>

grep -A 1 -f OH_GB_030416_rbcL_TestSet.txt OH_Spcs_rbcL_400_1500bp_CLEAN_2.fasta > OH_rbcL_TestSet.fasta</span>

 make a fasta file with the header and sequence for the training set entries</span>

grep -A 1 -f OH_GB_030416_rbcL_TrainSet.txt OH_Spcs_rbcL_200_800bp_CLEAN_2.fasta > OH_rbcL_TrainSet.fasta</span>

 delete grep “--” artifact</span>

sed -i '/--/d' OH_rbcL_*Set.fasta</span>

# 4\. Format reference sequences and train classifiers</span>

## 4.1 Format training sequences</span>

 get Gi numbers from training sequences</span>

grep '>' OH_rbcL_TrainSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_ohio.gis</span>

 use the following two Perl scripts provided in Sickel at al. (2015) and the NCBI Taxonomy module to produce required taxonomy and fasta files for later use with RDP and UTAX. See https://github.com/iimog/meta-barcoding-dual-indexing for the Perl scripts and instructions on installing the taxonomy module</span>

perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_ohio.gis --out rbcL_ohio.tax --species rbcL_ohio.species.taxids --genus rbcL_ohio.genus.taxids</span>

perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/tax2rdp_utax.pl rbcL_ohio.tax OH_rbcL_TrainSet.fasta OH_rbcL_trainRDP</span>

 remove duplicate sequences</span>

java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-dupseq --infile OH_rbcL_trainRDP.rdp.fa --outfile OH_rbcL_trainRDP.rmDS.rdp.fa --duplicates --min_seq_length 50</span>

 remove partial sequences</span>

java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-partialseq OH_rbcL_trainRDP.rdp.fa OH_rbcL_trainRDP.rmDS.rdp.fa OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --alignment-mode overlap --min_gaps 50 --knn 20</span>

 remove entries with incomplete taxonomy data</span>

sed -i '/undef/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

sed -i '/ sp/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

sed -i '/incertae/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

sed -i '/incerti/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

sed -i '/unknown/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

## 4.2 Train RDP classifier</span>

java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar train --out_dir rbcL_trained --seq OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --tax_file OH_rbcL_trainRDP.rdp.tax</span>

 move rbcL.properties, produced during training, from data directory to rdp_trained directory</span>

mv ./rbcL.properties ./rbcL_trained/rbcL.properties</span>

## 4.3 Train UTAX classifier</span>

 use training sequences curated for RDP, OH_rbcL_trainRDP.FinalTrainSet.rdp.fa from previous page, as downstream classifier training data for UTAX and as reference data for RTAX</span>

 reformat fasta headers from RDP format to UTAX format</span>

cat OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)\tRoot;/>$1;tax=/g' | perl -pe 's/(.)__/$1:/g' | perl -pe 's/_\d+;/,/g' | perl -pe 's/,$/;/g' > rbcL.utax.fa</span>

 make UTAX training database</span>

~/PATH/TO/Usearch/usearch_v8.1 -makeudb_utax rbcL.utax.fa -output rbcL.udb</span>

## 4.4 Re-format references for RTAX classifier</span>

 reformat the taxonomic lineage headers of the RDP training fasta file to make a taxonomy file which meets the requirements of the RTAX classifier</span>

grep '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)(\t)Root;/$1$2/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__/ /g' | perl -pe 's/_\d+//g' > rbcL_training.rtax.tax</span>

 reformat RDP training fasta file for use with RTAX</span>

grep -A 1 '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+).*/>$1/g' > rbcL_training.rtax.fasta</span>

# 5\. Format testing sequences</span>

## 5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages</span>

 get file of Gi numbers for each test sequence entry</span>

grep '>' OH_rbcL_TestSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_test_ohio.gis</span>

 use Gi numbers along with the Perl script from Sickel et al. (2015) and the NCBI taxonomy module to get the taxonomic lineage of each entry</span>

perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_test_ohio.gis --out rbcL_test_ohio.tax --species rbcL_test_ohio.species.taxids --genus rbcL_test_ohio.genus.taxids</span>

</span>

## 5.2 Replace Gi numbers in fasta with lineages</span>

 use the following bash script to replace fasta Gi headers with the taxonomic lineage for each entry</span>

#!/bin/bash</span>

cat rbcL_test_ohio.gis | while read Gi</span>

do</span>

grep -w "$Gi" rbcL_test_ohio.tax | perl -pe 's/(\d+)\tRoot;(.*)/>$1;$2/g' >> rbcL_test_RFMT.fasta</span>

grep -A 1 -w "$Gi" OH_rbcL_TestSet.fasta | grep -v '>' >> rbcL_test_RFMT.fasta  </span>

        done</span>

## 5.3 Format sequences to facilitate later analyses</span>

 remove taxonomically incomplete entries</span>

sed -i '/undef/,+1 d' rbcL_test_RFMT.fasta  </span>

sed -i '/ sp/,+1 d' rbcL_test_RFMT.fasta  </span>

sed -i '/incertae/,+1 d' rbcL_test_RFMT.fasta  </span>

sed -i '/incerti/,+1 d' rbcL_test_RFMT.fasta  </span>

sed -i '/unknown/,+1 d' rbcL_test_RFMT.fasta  </span>

# 6\. Classify testing sequences</span>

## 6.1 Classify sequences using RDP</span>

</span> replace tabs in header with “-“ so classifier doesn't delete the characters after the space.</span>

grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/         /-/g' > rbcL_test_rdpFMT.fasta</span>

java -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar classify -t ~/PATH/TO/rbcL.properties --outputFile rbcL.rdp.output rbcL_test_rdpFMT.fasta</span>

## 6.2 Classify sequences using RTAX</span>

 replace “;” with “–“ in header so that whole lineage is retained in RTAX output</span>

grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/;/-/g' > rbcL_test_rtaxFMT.fasta</span>

 classify sequences</span>

~/PATH/TO/rtax -r rbcL_training.rtax.fasta -t rbcL_training.rtax.tax -a rbcL_rtaxFMT.fasta -o rbcL.rtax.output</span>

## 6.3 Classify sequences using UTAX</span>

~/local/src/Usearch/usearch_v8.1 -utax rbcL_test_RFMT.fasta -db ~/PATH/TO/rbcL.udb -rdpout rbcL.utax.output -strand plus</span>

# 7\. Format classifier output and perform accuracy and sensitivity analysis</span>

## 7.1 Format RDP output</span>

 replace '\t-\t' with '\t' and replace tabs with commas</span>

cat rbcL.rdp.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_clean.rdp.output</span>

</span>

 replace RDP assignments with less than 90% bootstrap support with “NULL” using the following bash script (this facilitates later analysis in R)</span>

</span>

#!/bin/bash</span>

 this for loop can be applied to all five loci if all testing sequences have been classified</span>

for i in rbcL</span>

do</span>

cat ${i}_clean.rdp.output | while read LINE</span>

do</span>

ACTUAL_LINEAGE_FINAL=$(echo $LINE | cut -d',' -f1)</span>

KINGDOM=$(echo $LINE | cut -d',' -f5,6,7)</span>

PHYLUM=$(echo $LINE | cut -d',' -f8,9,10)</span>

CLASS=$(echo $LINE | cut -d',' -f11,12,13)</span>

ORDER=$(echo $LINE | cut -d',' -f14,15,16)</span>

FAMILY=$(echo $LINE | cut -d',' -f17,18,19)</span>

GENUS=$(echo $LINE | cut -d',' -f20,21,22)</span>

SPECIES=$(echo $LINE | cut -d',' -f23,24,25)</span>

</span>

K_VAL=$(echo $KINGDOM | cut -d',' -f3)</span>

P_VAL=$(echo $PHYLUM | cut -d',' -f3)</span>

C_VAL=$(echo $CLASS | cut -d',' -f3)</span>

O_VAL=$(echo $ORDER | cut -d',' -f3)</span>

F_VAL=$(echo $FAMILY | cut -d',' -f3)</span>

G_VAL=$(echo $GENUS | cut -d',' -f3)</span>

S_VAL=$(echo $SPECIES | cut -d',' -f3)</span>

echo {$K,$P,$C,$O,$G,$S}_VAL >> QC_CHECK_${i}</span>

if [ $K_VAL \< 0.90 ]</span>

then</span>

K_FINAL=NULL</span>

else</span>

K_FINAL=$(echo $LINE | cut -d',' -f5,6,7)</span>

fi</span>

if [ $P_VAL \< 0.90 ]</span>

then</span>

P_FINAL=NULL</span>

else</span>

P_FINAL=$(echo $LINE | cut -d',' -f8,9,10)</span>

firbcL_0.5DB_OUT</span>

if [ $C_VAL \< 0.90 ]</span>

then</span>

C_FINAL=NULL</span>

else</span>

C_FINAL=$(echo $LINE | cut -d',' -f11,12,13)</span>

fi</span>

if [ $O_VAL \< 0.90 ]</span>

thenmatK_Test_OUTPUT</span>

O_FINAL=NULL</span>

else</span>

O_FINAL=$(echo $LINE | cut -d',' -f14,15,16)</span>

fi</span>

if [ $F_VAL \< 0.90 ]</span>

then</span>

F_FINAL=NULL</span>

else</span>

F_FINAL=$(echo $LINE | cut -d',' -f17,18,19)</span>

fi</span>

if [ $G_VAL \< 0.90 ]</span>

then</span>

G_FINAL=NULL</span>

else</span>

G_FINAL=$(echo $LINE | cut -d',' -f20,21,22)</span>

fi</span>

if [ $S_VAL \< 0.90 ]</span>

then</span>

S_FINAL=NULL</span>

else</span>

S_FINAL=$(echo $LINE | cut -d',' -f23,24,25)</span>

fi</span>

echo {$ACTUAL_LINEAGE,$K,$P,$C,$O,$F,$G,$S}_FINAL >> ${i}_RDP_OUTPUT</span>

done</span>

done</span>

</span>

 miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RDP assignments (the predicted identity)</span>

cat rbcL_RDP_OUTPUT | perl –pe ‘s/NULL/NULL,NULL,NULL/g’ | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi'</span> > rbcL_RDP_FINAL</span>

</span>

## 7.2 Format RTAX output</span>

 miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RTAX assignments (the predicted identity)</span>

cat rbcL.rtax.output | perl -pe 's/(\d+)-(.__)/$1,$2/g' | perl -pe 's/\t/,/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' | perl -pe 's/, /,/g' | perl -pe 's/-//g' > rbcL_RTAX_FINAL</span>

</span>

 import data into excel and replace empty cells with “NULL”, save as .csv and import to R</span>

## 7.3 Format UTAX output</span>

 miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the UTAX assignments (the predicted identity)</span>

cat rbcL.utax.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_UTAX_output</span>

</span>

 use bash shell from RDP output reformatting (section 7.1) to remove assignments below 0.90 confidence and replace empty cells with “NULL”</span>

</span>

 further reformatting for later analysis</span>

cat rbcL_UTAX_OUTPUT</span> | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g'</span> > rbcL_UTAX_FINAL</span>

## 7.4 Perform accuracy and sensitivity analysis</span>

 isolate classified sequences to determine sensitivity in R (the example given is for RTAX at the genus level but this would have to be changed for other ranks and classifiers)</span>

genus_assigned <- as.data.frame(rbcL_RTAX_FINAL[which(rbcL_RTAX_FINAL $V23 != 'NULL'),c(7,23)])</span>

 compare testing sequence taxonomy with classifier predicted taxonomy to isolate mis-classifications in R</span>

genus_misIDs <- as.data.frame(rbcL_RTAX_FINAL[which(as.character(rbcL_RTAX_FINAL$V7) != as.character(rbcL_RTAX_FINAL$V15) & as.character(rbcL_RTAX_FINAL$V15) != 'NULL'),c(7,15)])</span>

 visually analyze the resulting data frames to ensure that the correct columns were compared and that there are no obviously confounding artifacts in the data</span>

# 8\. Calculate DB coverage</span>

</span> Obtain lists of genera from Ohio and both the training sets and the testing sets</span>

cat OhioSpecies.txt | cut -d' ' -f1 | sort | uniq > rbcL_Ohio_genera</span>

cat rbcL_test_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi’ | sort | uniq  > rbcL_test_genera</span>

cat rbcL_train_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi' | sort | uniq > rbcL_train_genera</span>

 import these to excel and save as .csv</span>

 use %in% function in R to find overlap of genera lists</span>

rbcL_test_train <- as.data.frame(rbcL_test[which(rbcL_test[,1] %in% rbcL_train[,1]),1])</span>

</span>

# 9\. Isolate sequences with no genus level reference representation</span>

 find genera present in the testing set but not in the training set</span>

grep -v -f rbcL_train_genera rbcL_test_genera > rbcL_in_test_not_train</span>
