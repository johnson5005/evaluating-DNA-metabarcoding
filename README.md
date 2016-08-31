###<span class="c4 c0 c22">Appendix S1</span>

<span class="c12 c24 c0">Evaluating the Relative Performance of DNA Metabarcoding Sequence Classifiers</span>

<span class="c12 c0 c24">Rodney T. Richardson, Johan Bengtsson-Palme and Reed M. Johnson</span>

<span class="c15 c4 c0">Table of Contents</span>

<span class="c3 c0">[1\. Introduction](#h.gjdgxs)</span><span class="c5">[        ](#h.gjdgxs)</span>

<span class="c3 c0">[2\. Retrieve data and filter by Ohio plant species list](#h.30j0zll)</span><span class="c5">[        ](#h.30j0zll)</span>

<span class="c3 c0">[2.1 Obtain data from NCBI](#h.1fob9te)</span><span class="c5">[        ](#h.1fob9te)</span>

<span class="c3 c0">[2.2 Preliminary formatting and quality curation](#h.3znysh7)</span><span class="c5">[        ](#h.3znysh7)</span>

<span class="c3 c0">[2.3 Filter by geography](#h.2et92p0)</span><span class="c5">[        ](#h.2et92p0)</span>

<span class="c3 c0">[3\. Randomly sample testing sequences from training sequences](#h.tyjcwt)</span><span class="c5">[        ](#h.tyjcwt)</span>

<span class="c3 c0">[4\. Format reference sequences and train classifiers](#h.3dy6vkm)</span><span class="c5">[        ](#h.3dy6vkm)</span>

<span class="c3 c0">[4.1 Format training sequences](#h.1t3h5sf)</span><span class="c5">[        ](#h.1t3h5sf)</span>

<span class="c3 c0">[4.2 Train RDP classifier](#h.4d34og8)</span><span class="c5">[        ](#h.4d34og8)</span>

<span class="c3 c0">[4.3 Train UTAX classifier](#h.2s8eyo1)</span><span class="c5">[        ](#h.2s8eyo1)</span>

<span class="c3 c0">[4.4 Re-format references for RTAX classifier](#h.17dp8vu)</span><span class="c5">[        ](#h.17dp8vu)</span>

<span class="c3 c0">[5\. Format testing sequences](#h.3rdcrjn)</span><span class="c5">[        ](#h.3rdcrjn)</span>

<span class="c0 c3">[5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages](#h.26in1rg)</span><span class="c5">[        ](#h.26in1rg)</span>

<span class="c3 c0">[5.2 Replace Gi numbers in fasta with lineages](#h.lnxbz9)</span><span class="c5">[        ](#h.lnxbz9)</span>

<span class="c3 c0">[5.3 Format sequences to facilitate later analyses](#h.35nkun2)</span><span class="c5">[        ](#h.35nkun2)</span>

<span class="c3 c0">[6\. Classify testing sequences](#h.1ksv4uv)</span><span class="c5">[        ](#h.1ksv4uv)</span>

<span class="c3 c0">[6.1 Classify sequences using RDP](#h.44sinio)</span><span class="c5">[        ](#h.44sinio)</span>

<span class="c3 c0">[6.2 Classify sequences using RTAX](#h.z337ya)</span><span class="c5">[        ](#h.z337ya)</span>

<span class="c3 c0">[6.3 Classify sequences using UTAX](#h.3j2qqm3)</span><span class="c5">[        ](#h.3j2qqm3)</span>

<span class="c3 c0">[7\. Format classifier output and perform accuracy and sensitivity analysis](#h.1y810tw)</span><span class="c5">[        ](#h.1y810tw)</span>

<span class="c3 c0">[7.1 Format RDP output](#h.4i7ojhp)</span><span class="c5">[        ](#h.4i7ojhp)</span>

<span class="c3 c0">[7.2 Format RTAX output](#h.1ci93xb)</span><span class="c5">[        ](#h.1ci93xb)</span>

<span class="c3 c0">[7.3 Format UTAX output](#h.2bn6wsx)</span><span class="c5">[        ](#h.2bn6wsx)</span>

<span class="c3 c0">[7.4 Perform accuracy and sensitivity analysis](#h.3as4poj)</span><span class="c5">[        ](#h.3as4poj)</span>

<span class="c3 c0">[8\. Calculate DB coverage](#h.1pxezwc)</span><span class="c5">[        ](#h.1pxezwc)</span>

<span class="c3 c0">[9\. Isolate sequences with no genus level reference representation](#h.49x2ik5)</span><span class="c5">[        ](#h.49x2ik5)</span>

<span class="c4">[](#_Toc453838845)</span>

<span>[](#_Toc453838845)</span>

# <span class="c0">1\. Introduction</span>

<span class="c12 c0">This document contains the commands used for evaluation of classifier performance. These commands are provided such that readers can work through our analyses independently and apply the approach to their own research endeavors. It should be noted however that the syntax, commands and software used here may not be entirely transferrable for future applications given differences in computational architecture, software updates, etc. Further, these commands are given without guidance in terms of directory organization, which we leave at the discretion of the reader. Additionally, due to the ease of transferring commands from one analysis to another, we do not provide commands for all the analyses performed in the paper. Lastly, while these commands were performed on files representing five loci, for ease of organization we generally give only the commands used for the</span> <span class="c12 c0 c32">rbcL</span><span class="c12 c0"> locus as an example.</span>

# <span class="c0">2\. Retrieve data and filter by Ohio plant species list</span>

## <span class="c12 c0">2.1 Obtain data from NCBI</span>

<span>### search the following in NCBI</span>

*   <span class="c0">ITS2 AND vascular plants [ORGN]</span>
*   <span class="c0">rbcL AND vascular plants [ORGN]</span>
*   <span class="c0">matK AND vascular plants [ORGN]</span>
*   <span class="c0">trnH AND psbA AND vascular plants [ORGN]</span>
*   <span class="c0">trnL AND vascular plants [ORGN]</span>

## <span class="c12 c0">2.2 Preliminary formatting and quality curation</span>

<span class="c0">### remove newline instances from NCBI formatted fasta files</span>

<span class="c0">awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNextLine.fasta</span>

<span class="c0">### filter by minimum length threshold</span>

<span class="c0">awk '!/^>/ { next } { getline seq } length(seq) >= 400 { print $0 "\n" seq }' GenBank_rbcL_030416_rmNextLine.fasta > GenBank_rbcL_030416_400.fasta</span>

<span class="c0">### filter by maximum length threshold</span>

<span class="c0">awk '!/^>/ { next } { getline seq } length(seq) <= 1500 { print $0 "\n" seq }' GenBank_rbcL_030416_500.fasta > GenBank_rbcL_030416_400_1500bp.fasta</span>

<span class="c0">### Length trimming thresholds used</span>

*   <span class="c0">200 – 800 bp for ITS2</span>
*   <span class="c0">500 – 1000 bp for matK</span>
*   <span class="c0">400 – 1500 bp for rbcL</span>
*   <span class="c0">150 – 1100 bp for trnL</span>
*   <span class="c0">200 – 800 bp for trnH</span>

<span class="c0">### remove sequences with more than two consecutive uncalled base pairs</span>

<span class="c0">sed -i '/NN/I,+1 d' GenBank_rbcL_030416_400_1500bp.fasta</span>

## <span class="c12 c0">2.3 Filter by geography</span>

<span class="c0">### get Gi numbers from NCBI entries for species found in Ohio (OhioSpecies.txt file is provided in Appendix S2). It should be noted that this approach to filtering the sequences is conservative as only one sequence representative per species is retained. Furthermore, one should be aware of that there could be naming inconsistencies between the list of species used for filtering and, e.g., the NCBI naming convention. Thus, this step should only be carried out if relevant for the questions investigated.</span>

<span class="c0">grep -f OhioSpecies.txt GenBank_rbcL_030416_400_1500bp.fasta | awk '!seen[$2$3]++' | perl -pe 's/^>gi\|(\d+)\|.*/$1/' > OH_GB_030416_rbcL.txt</span>

<span class="c0">### make blast database with reference sequences</span>

<span class="c0">~/PATH/TO/ncbi-blast-2.2.29+/bin/makeblastdb -in GenBank_rbcL_030416_400_1500bp.fasta -dbtype nucl -parse_seqids</span>

<span class="c0">### use blastdbcmd to retrieve sequences matching the Gi numbers of plants known to occur in Ohio</span>

<span class="c0">~/PATH/TO/ncbi-blast-2.2.29+/bin/blastdbcmd -db GenBank_rbcL_030416_400_1500bp.fasta -entry_batch OH_GB_030416_rbcL.txt -out OH_Spcs_rbcL.fa</span>

<span class="c0">### remove everything from sequence headers except “>Gi” to meet formatting requirements of downstream applications</span>

<span class="c0">cat OH_Spcs_rbcL.fa | perl -pe 's/^>gi\|(\d+)\|.*/>$1/g' > OH_Spcs_rbcL_400_1500bp_CLEAN.fasta</span>

<span class="c0">### remove next line instances which were added during blast database usage</span>

<span class="c0">awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNL.fasta</span>

# <span class="c0">3\. Randomly sample testing sequences from training sequences</span>

<span class="c0">### sample the headers from 15% of the sequences to serve as the testing set</span>

<span class="c0">awk 'BEGIN {srand()} !/^$/ { if (rand() <= .15) print $0}' < OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TestSet.txt</span>

<span class="c0">### use reverse grep to get headers of sequences which were not sampled for the testing set</span>

<span class="c0">grep -v -f OH_GB_030416_rbcL_TestSet.txt OH_GB_030416_rbcL.txt > OH_GB_030416_rbcL_TrainSet.txt</span>

<span class="c0"></span>

<span class="c0"></span>

<span class="c0">### make a fasta file with the header and sequence for the testing set entries</span>

<span class="c0">grep -A 1 -f OH_GB_030416_rbcL_TestSet.txt OH_Spcs_rbcL_400_1500bp_CLEAN_2.fasta > OH_rbcL_TestSet.fasta</span>

<span class="c0">### make a fasta file with the header and sequence for the training set entries</span>

<span class="c0">grep -A 1 -f OH_GB_030416_rbcL_TrainSet.txt OH_Spcs_rbcL_200_800bp_CLEAN_2.fasta > OH_rbcL_TrainSet.fasta</span>

<span class="c0">### delete grep “--” artifact</span>

<span class="c0">sed -i '/--/d' OH_rbcL_*Set.fasta</span>

# <span class="c0">4\. Format reference sequences and train classifiers</span>

## <span class="c12 c0">4.1 Format training sequences</span>

<span class="c0">### get Gi numbers from training sequences</span>

<span class="c0">grep '>' OH_rbcL_TrainSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_ohio.gis</span>

<span class="c0">### use the following two Perl scripts provided in Sickel at al. (2015) and the NCBI Taxonomy module to produce required taxonomy and fasta files for later use with RDP and UTAX. See https://github.com/iimog/meta-barcoding-dual-indexing for the Perl scripts and instructions on installing the taxonomy module</span>

<span class="c0">perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_ohio.gis --out rbcL_ohio.tax --species rbcL_ohio.species.taxids --genus rbcL_ohio.genus.taxids</span>

<span class="c0">perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/tax2rdp_utax.pl rbcL_ohio.tax OH_rbcL_TrainSet.fasta OH_rbcL_trainRDP</span>

<span class="c0">### remove duplicate sequences</span>

<span class="c0">java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-dupseq --infile OH_rbcL_trainRDP.rdp.fa --outfile OH_rbcL_trainRDP.rmDS.rdp.fa --duplicates --min_seq_length 50</span>

<span class="c0">### remove partial sequences</span>

<span class="c0">java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-partialseq OH_rbcL_trainRDP.rdp.fa OH_rbcL_trainRDP.rmDS.rdp.fa OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --alignment-mode overlap --min_gaps 50 --knn 20</span>

<span class="c0">### remove entries with incomplete taxonomy data</span>

<span class="c0">sed -i '/undef/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

<span class="c0">sed -i '/ sp/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

<span class="c0">sed -i '/incertae/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

<span class="c0">sed -i '/incerti/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

<span class="c0">sed -i '/unknown/,+1 d' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa</span>

## <span class="c12 c0">4.2 Train RDP classifier</span>

<span class="c0">java -Xmx4g -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar train --out_dir rbcL_trained --seq OH_rbcL_trainRDP.FinalTrainSet.rdp.fa --tax_file OH_rbcL_trainRDP.rdp.tax</span>

<span class="c0">### move rbcL.properties, produced during training, from data directory to rdp_trained directory</span>

<span class="c0">mv ./rbcL.properties ./rbcL_trained/rbcL.properties</span>

## <span class="c12 c0">4.3 Train UTAX classifier</span>

<span class="c0">### use training sequences curated for RDP, OH_rbcL_trainRDP.FinalTrainSet.rdp.fa from previous page, as downstream classifier training data for UTAX and as reference data for RTAX</span>

<span class="c0">### reformat fasta headers from RDP format to UTAX format</span>

<span class="c0">cat OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)\tRoot;/>$1;tax=/g' | perl -pe 's/(.)__/$1:/g' | perl -pe 's/_\d+;/,/g' | perl -pe 's/,$/;/g' > rbcL.utax.fa</span>

<span class="c0">### make UTAX training database</span>

<span class="c0">~/PATH/TO/Usearch/usearch_v8.1 -makeudb_utax rbcL.utax.fa -output rbcL.udb</span>

## <span class="c12 c0">4.4 Re-format references for RTAX classifier</span>

<span class="c0">### reformat the taxonomic lineage headers of the RDP training fasta file to make a taxonomy file which meets the requirements of the RTAX classifier</span>

<span class="c0">grep '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+)(\t)Root;/$1$2/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__/ /g' | perl -pe 's/_\d+//g' > rbcL_training.rtax.tax</span>

<span class="c0">### reformat RDP training fasta file for use with RTAX</span>

<span class="c0">grep -A 1 '>' OH_rbcL_trainRDP.FinalTrainSet.rdp.fa | perl -pe 's/^>(\d+).*/>$1/g' > rbcL_training.rtax.fasta</span>

# <span class="c0">5\. Format testing sequences</span>

## <span class="c12 c0">5.1 Get sequence Gi numbers and use them to retrieve NCBI Taxonomy database lineages</span>

<span class="c0">### get file of Gi numbers for each test sequence entry</span>

<span class="c0">grep '>' OH_rbcL_TestSet.fasta | perl -pe 's/>(\d+).*/$1/g' > rbcL_test_ohio.gis</span>

<span class="c0">### use Gi numbers along with the Perl script from Sickel et al. (2015) and the NCBI taxonomy module to get the taxonomic lineage of each entry</span>

<span class="c0">perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis rbcL_test_ohio.gis --out rbcL_test_ohio.tax --species rbcL_test_ohio.species.taxids --genus rbcL_test_ohio.genus.taxids</span>

<span class="c0"></span>

## <span class="c12 c0">5.2 Replace Gi numbers in fasta with lineages</span>

<span class="c0">### use the following bash script to replace fasta Gi headers with the taxonomic lineage for each entry</span>

<span class="c0">#!/bin/bash</span>

<span class="c0">cat rbcL_test_ohio.gis | while read Gi</span>

<span class="c0">do</span>

<span class="c0">grep -w "$Gi" rbcL_test_ohio.tax | perl -pe 's/(\d+)\tRoot;(.*)/>$1;$2/g' >> rbcL_test_RFMT.fasta</span>

<span class="c0">grep -A 1 -w "$Gi" OH_rbcL_TestSet.fasta | grep -v '>' >> rbcL_test_RFMT.fasta  </span>

<span class="c0">        done</span>

## <span class="c12 c0">5.3 Format sequences to facilitate later analyses</span>

<span class="c0">### remove taxonomically incomplete entries</span>

<span class="c0">sed -i '/undef/,+1 d' rbcL_test_RFMT.fasta  </span>

<span class="c0">sed -i '/ sp/,+1 d' rbcL_test_RFMT.fasta  </span>

<span class="c0">sed -i '/incertae/,+1 d' rbcL_test_RFMT.fasta  </span>

<span class="c0">sed -i '/incerti/,+1 d' rbcL_test_RFMT.fasta  </span>

<span class="c0">sed -i '/unknown/,+1 d' rbcL_test_RFMT.fasta  </span>

# <span class="c0">6\. Classify testing sequences</span>

## <span class="c12 c0">6.1 Classify sequences using RDP</span>

<span class="c0">###</span> <span class="c0 c12">replace tabs in header with “-“ so classifier doesn't delete the characters after the space.</span>

<span class="c0">grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/         /-/g' > rbcL_test_rdpFMT.fasta</span>

<span class="c0">java -jar ~/PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar classify -t ~/PATH/TO/rbcL.properties --outputFile rbcL.rdp.output rbcL_test_rdpFMT.fasta</span>

## <span class="c12 c0">6.2 Classify sequences using RTAX</span>

<span class="c0">### replace “;” with “–“ in header so that whole lineage is retained in RTAX output</span>

<span class="c0">grep -A 1 '>' rbcL_test_RFMT.fasta | perl -pe 's/;/-/g' > rbcL_test_rtaxFMT.fasta</span>

<span class="c0">### classify sequences</span>

<span class="c0">~/PATH/TO/rtax -r rbcL_training.rtax.fasta -t rbcL_training.rtax.tax -a rbcL_rtaxFMT.fasta -o rbcL.rtax.output</span>

## <span class="c12 c0">6.3 Classify sequences using UTAX</span>

<span class="c0 c23">~/local/src/Usearch/usearch_v8.1 -utax rbcL_test_RFMT.fasta -db ~/PATH/TO/rbcL.udb -rdpout rbcL.utax.output -strand plus</span>

# <span class="c0">7\. Format classifier output and perform accuracy and sensitivity analysis</span>

## <span class="c12 c0">7.1 Format RDP output</span>

<span class="c0">### replace '\t-\t' with '\t' and replace tabs with commas</span>

<span class="c0">cat rbcL.rdp.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_clean.rdp.output</span>

<span class="c0"></span>

<span class="c0">### replace RDP assignments with less than 90% bootstrap support with “NULL” using the following bash script (this facilitates later analysis in R)</span>

<span class="c0"></span>

<span class="c0">#!/bin/bash</span>

<span class="c0">### this for loop can be applied to all five loci if all testing sequences have been classified</span>

<span class="c0">for i in rbcL</span>

<span class="c0">do</span>

<span class="c0">cat ${i}_clean.rdp.output | while read LINE</span>

<span class="c0">do</span>

<span class="c0">ACTUAL_LINEAGE_FINAL=$(echo $LINE | cut -d',' -f1)</span>

<span class="c0">KINGDOM=$(echo $LINE | cut -d',' -f5,6,7)</span>

<span class="c0">PHYLUM=$(echo $LINE | cut -d',' -f8,9,10)</span>

<span class="c0">CLASS=$(echo $LINE | cut -d',' -f11,12,13)</span>

<span class="c0">ORDER=$(echo $LINE | cut -d',' -f14,15,16)</span>

<span class="c0">FAMILY=$(echo $LINE | cut -d',' -f17,18,19)</span>

<span class="c0">GENUS=$(echo $LINE | cut -d',' -f20,21,22)</span>

<span class="c0">SPECIES=$(echo $LINE | cut -d',' -f23,24,25)</span>

<span class="c0"></span>

<span class="c0">K_VAL=$(echo $KINGDOM | cut -d',' -f3)</span>

<span class="c0">P_VAL=$(echo $PHYLUM | cut -d',' -f3)</span>

<span class="c0">C_VAL=$(echo $CLASS | cut -d',' -f3)</span>

<span class="c0">O_VAL=$(echo $ORDER | cut -d',' -f3)</span>

<span class="c0">F_VAL=$(echo $FAMILY | cut -d',' -f3)</span>

<span class="c0">G_VAL=$(echo $GENUS | cut -d',' -f3)</span>

<span class="c0">S_VAL=$(echo $SPECIES | cut -d',' -f3)</span>

<span class="c0">echo {$K,$P,$C,$O,$G,$S}_VAL >> QC_CHECK_${i}</span>

<span class="c0">if [ $K_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">K_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">K_FINAL=$(echo $LINE | cut -d',' -f5,6,7)</span>

<span class="c0">fi</span>

<span class="c0">if [ $P_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">P_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">P_FINAL=$(echo $LINE | cut -d',' -f8,9,10)</span>

<span class="c0">firbcL_0.5DB_OUT</span>

<span class="c0">if [ $C_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">C_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">C_FINAL=$(echo $LINE | cut -d',' -f11,12,13)</span>

<span class="c0">fi</span>

<span class="c0">if [ $O_VAL \< 0.90 ]</span>

<span class="c0">thenmatK_Test_OUTPUT</span>

<span class="c0">O_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">O_FINAL=$(echo $LINE | cut -d',' -f14,15,16)</span>

<span class="c0">fi</span>

<span class="c0">if [ $F_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">F_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">F_FINAL=$(echo $LINE | cut -d',' -f17,18,19)</span>

<span class="c0">fi</span>

<span class="c0">if [ $G_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">G_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">G_FINAL=$(echo $LINE | cut -d',' -f20,21,22)</span>

<span class="c0">fi</span>

<span class="c0">if [ $S_VAL \< 0.90 ]</span>

<span class="c0">then</span>

<span class="c0">S_FINAL=NULL</span>

<span class="c0">else</span>

<span class="c0">S_FINAL=$(echo $LINE | cut -d',' -f23,24,25)</span>

<span class="c0">fi</span>

<span class="c0">echo {$ACTUAL_LINEAGE,$K,$P,$C,$O,$F,$G,$S}_FINAL >> ${i}_RDP_OUTPUT</span>

<span class="c0">done</span>

<span class="c0">done</span>

<span class="c0"></span>

<span class="c0">### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RDP assignments (the predicted identity)</span>

<span class="c0">cat rbcL_RDP_OUTPUT | perl –pe ‘s/NULL/NULL,NULL,NULL/g’ | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi'</span> <span class="c0">> rbcL_RDP_FINAL</span>

<span class="c0"></span>

## <span class="c12 c0">7.2 Format RTAX output</span>

<span class="c0">### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the RTAX assignments (the predicted identity)</span>

<span class="c0">cat rbcL.rtax.output | perl -pe 's/(\d+)-(.__)/$1,$2/g' | perl -pe 's/\t/,/g' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' | perl -pe 's/, /,/g' | perl -pe 's/-//g' > rbcL_RTAX_FINAL</span>

<span class="c0"></span>

<span class="c0">### import data into excel and replace empty cells with “NULL”, save as .csv and import to R</span>

## <span class="c12 c0">7.3 Format UTAX output</span>

<span class="c0">### miscellaneous reformatting to ensure that taxonomic values in the header (the actual taxonomic identity) match the values in the UTAX assignments (the predicted identity)</span>

<span class="c0">cat rbcL.utax.output | perl -pe 's/(\t)-(\t)/$1/g' | perl -pe 's/\t/,/g' | perl -pe 's/,,/,/g' | perl -pe 's/([a-z])-([a-z])/$1 $2/g' > rbcL_UTAX_output</span>

<span class="c0"></span>

<span class="c0">### use bash shell from RDP output reformatting (section 7.1) to remove assignments below 0.90 confidence and replace empty cells with “NULL”</span>

<span class="c0"></span>

<span class="c0">### further reformatting for later analysis</span>

<span class="c0">cat rbcL_UTAX_OUTPUT</span><span class="c0"> | perl -pe 's/(\d) ([a-z])/$1,$2/gi' | perl -pe 's/NULL ([a-z]+)/NULL,$1/gi' | perl -pe 's/([a-z]+) NULL/$1,NULL/gi' | perl -pe 's/sub__/sub_/g' | perl -pe 's/.__//g' | perl -pe 's/_\d+//g'</span><span class="c0"> > rbcL_UTAX_FINAL</span>

## <span class="c12 c0">7.4 Perform accuracy and sensitivity analysis</span>

<span class="c0">### isolate classified sequences to determine sensitivity in R (the example given is for RTAX at the genus level but this would have to be changed for other ranks and classifiers)</span>

<span class="c0">genus_assigned <- as.data.frame(rbcL_RTAX_FINAL[which(rbcL_RTAX_FINAL $V23 != 'NULL'),c(7,23)])</span>

<span class="c0">### compare testing sequence taxonomy with classifier predicted taxonomy to isolate mis-classifications in R</span>

<span class="c0">genus_misIDs <- as.data.frame(rbcL_RTAX_FINAL[which(as.character(rbcL_RTAX_FINAL$V7) != as.character(rbcL_RTAX_FINAL$V15) & as.character(rbcL_RTAX_FINAL$V15) != 'NULL'),c(7,15)])</span>

<span class="c0">### visually analyze the resulting data frames to ensure that the correct columns were compared and that there are no obviously confounding artifacts in the data</span>

# <span class="c0">8\. Calculate DB coverage</span>

<span class="c0">###</span> <span class="c12 c0">Obtain lists of genera from Ohio and both the training sets and the testing sets</span>

<span class="c0">cat OhioSpecies.txt | cut -d' ' -f1 | sort | uniq > rbcL_Ohio_genera</span>

<span class="c0">cat rbcL_test_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi’ | sort | uniq  > rbcL_test_genera</span>

<span class="c0">cat rbcL_train_ohio.tax | perl -pe 's/.*;g__([a-z]+)_\d+;.*/$1/gi' | sort | uniq > rbcL_train_genera</span>

<span class="c0">### import these to excel and save as .csv</span>

<span class="c0">### use %in% function in R to find overlap of genera lists</span>

<span class="c0">rbcL_test_train <- as.data.frame(rbcL_test[which(rbcL_test[,1] %in% rbcL_train[,1]),1])</span>

<span class="c0"></span>

# <span class="c0">9\. Isolate sequences with no genus level reference representation</span>

<span class="c0">### find genera present in the testing set but not in the training set</span>

<span class="c0">grep -v -f rbcL_train_genera rbcL_test_genera > rbcL_in_test_not_train</span>
