
#Qiime2
source activate qiime2-2018.4

#DDB fasta
#Transforme les caracteres minuscule des sequences
fasta_formatter -t -i DAIRYdb_v1.1_10290_20180806_Final_STX.fasta -o tmp_DDBv1.1_ok.tsv
cut -f1 -d\; DDBv1.1_ok.tsv > tmp_seqIDs.txt
paste tmp_seqIDs.txt tmp_DDBv1.1_ok.tsv > DDBv1.1_ok.tsv

cat DDBv1.1_ok.tsv | awk -F "\t" '{print ">"$1 "\n" toupper($3)}' > DAIRYdb_v1.1_ok.fasta

#Taxonomy
grep "^>" DAIRYdb_v1.1_10290_20180806_Final_STX.fasta | cut -f2 -d\; |sed "s/tax=//g" #| sed "s/d:/k__/g;s/,p:/; p__/g;s/,c:/; c__/g;s/,o:/; o__/g;s/,f:/; f__/g;s/,g:/; g__/g;s/,s:/; s__/g" > tmp_tax2
grep "^>" DAIRYdb_v1.1_10290_20180806_Final_STX.fasta | cut -f1 -d\; | cut -c2- > tmp_tax1
paste tmp_tax1 tmp_tax2  > DDB_taxonomy.txt
rm tmp_*


#Qiime2 classifier train
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path DAIRYdb_v1.1_ok.fasta \
  --output-path DAIRYdb_v1.1_ok.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path DDB_taxonomy.txt \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads DAIRYdb_v1.1_ok.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier DDBv1.1_classifier.qza

