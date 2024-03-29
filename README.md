<img src="https://raw.githubusercontent.com/marcomeola/DAIRYdb/master/images/logo_ddb2.png" width="200">


# DAIRYdb: a manually curated reference database for improved taxonomy annotation of 16S rRNA gene sequences from dairy products
Marco Meola, Etienne Rifa, Noam Shani, Céline Delbes, Hélène Berthoud, Christophe Chassard. (2019) BMC Genomics 20(1):560. https://doi.org/10.1186/s12864-019-5914-8


## Update (april 2023): DAIRYdb v3.0.0

The new version v3.0.0 of DAIRYdb is [available here](https://nextcloud.inrae.fr/s/5ne3ooAsM7zQtcW). All changes are informed in the file `changes_v3.0.0.html`, main changes are listed in [#7](https://github.com/marcomeola/DAIRYdb/issues/7).

## Update (april 2022): DAIRYdb v2.1

The new version v2.1 of DAIRYdb is [available here](https://nextcloud.inrae.fr/s/5ne3ooAsM7zQtcW). See issue [#6](https://github.com/marcomeola/DAIRYdb/issues/6) for details.

## Update (march 2021): DAIRYdb v2.0
Latest versions of the DAIRYdb v1.2.5 (with the old nomenclature) and the **DAIRYdb v2.0** (containing the same sequences as v1.2.5 but with the new nomenclature for Lactic acid bacteria are [deposited here](https://nextcloud.inrae.fr/s/5ne3ooAsM7zQtcW). For differences between v1.2.5 and v2.0 see `changed_taxa_DBB2.0.txt`. See issue [#3](https://github.com/marcomeola/DAIRYdb/issues/3) for detailed updates.

## Description
DAIRYdb provides 10'439 sequences of full-length 16S ribosomal RNA (V1-V9) from microbial species (10'332 bacteria (50 more than v1.2.0), 107 archaea) of dairy products. The taxonomy has been automatically and manually curated on the 7 ranks.
DAIRYdb is able to assign sequences to the species rank whereas classical Databases are less accurate.


## Installation
### Download DAIRYdb
DAIRYdb_v1.2.4 is available here as newick tree file and adapted to different classification tools: Metax2.2, Blast+ and SINTAX, Qiime2, FROGS. For request to adapt DAIRYdb to other classifier please do not hesitate drop me a line at mmeola@imm.uzh.ch.

#### Sintax (Usearch32bit)
DAIRYdb_v1.2.4_20200604_STX.udb was generated using usearch v10.0. If the available .udb file is not working properly on your system it is recommended to recreate the .udb datbase with your usearch version and OS using following code:
```
usearch -makeudb_sintax DAIRYdb_v1.2.4_20200604_STX.fasta -output DAIRYdb_v1.2.4_20200604_STX.udb
```
Command to call the taxonomy predictor Sintax
```
usearch -sintax otus.fasta -db DAIRYdb_v1.2.4_20200604.udb -tabbedout out.sintax -strand both -sintax_cutoff 0.6
```
#### Metaxa2
SSU_DAIRYdb_v1.2.4_20200604_MTX was generated using Metaxa2 v2.2. If the available Metaxa2 adapted DAIRYdb SSU_DAIRYdb_v1.2.4_20200604_MTX is not working properly on your system it is recommended to recreate the Metaxa2 database with your Metaxa2 version and OS using following code:
```
metaxa2_dbb -o SSU_DAIRYdb_v1.2.4_20200604_MTX -g SSU_DAIRYdb_v1.2.4_20200604_MTX -t DAIRYdb_v1.2.4_20200604_TAX.txt --auto_rep T --cpu 4 --cutoffs 0,75,78.5,82,86.5,94.5,98.65 --save_raw T -a DAIRYdb_v1.2.4_20200604_Archaea.fasta -b DAIRYdb_v1.2.4_20200604_Bacteria.fasta --filter_uncultured F --correct_taxonomy F --evaluate F --plus T --divergent T
```
Unpack the tarball with
```
tar -xvfz SSU_DAIRYdb_v1.2.4_20200604_MTX.tar.gz
```
and copy the unpacked folder into the folder metaxa2_db (usually located at /usr/local/bin/metaxa2_db)
or follow the instructions on http://microbiology.se/2018/09/13/dairydb-added-to-metaxa2/

Command to call the taxonomy predictor Metaxa2.2 using the DAIRYdb
```
metaxa2 -i otus.fasta -g SSU_DAIRYdb_v1.2.4_20200604_MTX -o test --cpu 4 --taxonomy T --plus T -T 0,75,78.5,82,86.5,94.5,98.65 -taxlevel 7 -d blast -t b,a
```
#### Blast+
Database generated using Blast+
```
makeblastdb -in DAIRYdb_v1.2.4_20200604_blast.fasta -dbtype nucl
```
Command to call the taxonomy predictor Blast+
```
blastn -query otus.fasta -db DAIRYdb_v1.2.4_20200604_blast.fasta -num_threads 5 -out OUT_tax.txt -evalue 1 -outfmt 6 -perc_identity 97 -max_target_seqs 50
```
#### Qiime2
Database generated using Qiime2 classifier train
For more explanation check qiime2 tutorial (https://docs.qiime2.org/2018.6/tutorials/feature-classifier/)

##### Importing reference data sets
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path DAIRYdb_v1.2.4_ok.fasta \
  --output-path DAIRYdb_v1.2.4_ok.qza


qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path DDB_taxonomy.txt \
  --output-path ref-taxonomy.qza
```

##### Train the classifier
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads DAIRYdb_v1.2.4_.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier DAIRYdb_v1.2.4_20200604_qiime2_classifier.qza
```

##### Test the classifier
```
qiime feature-classifier classify-sklearn \
  --i-classifier DAIRYdb_v1.2.4_20200604_qiime2_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

#### IDTAXA
Open the R file DAIRYdb_v1.2.4_20200604_IDTAXA.R and run the commands or open the workspace with the trained classifier with DAIRYdb DAIRYdb_v1.2.4_20200604_IDTAXA.RData and import your fasta with the otus as described in the R script.

#### mothur
```
mothur "#classify.seqs(fasta=OTUS.fasta, template=DAIRYdb_v1.2.4_20200604_mothur.fasta , taxonomy=DAIRYdb_v1.2.4_20200604_mothur.tax)"
```

#### Kraken2
This version was added upon request although Kraken2 was developed for shotgun sequencing. We have no test or validation run with Kraken2 and the usage of DAIRYdb with Kraken2 goes without warranty.

```
kraken2 --db path_to/DAIRYdb_v1.2.4_20200604_kraken2 OTUS.fasta --use-names --report res.report > res.csv
```

## Usage recommendations for real samples

We recommend to use the taxonomy classification predicted coherently by both, Metaxa2 and SINTAX using the Excel file Taxonomy.template.xlsx. Classification errors should be reduced over selecting only coherent classification at any rank between both tools.

1) Classify your OTUs with Metaxa2 (see Metaxa2 manual for options)

Metaxa2 performance is highly influenced by the values given for classification in -T
```
metaxa2 -i otus.fasta -g DAIRYdb_v1.2.4_20200604_MTX -o out_metaxa2 --cpu 4 --taxonomy T --plus T -T 0,75,78.5,82,86.5,94.5,98.65 -taxlevel 7
```

2) Classify your OTUs with SINTAX
```
usearch -sintax otus.fasta -db DAIRYdb_v1.2.4_20200604.udb -tabbedout out.sintax -strand both -sintax_cutoff 0.6
```
Although lowering the sintax_cutoff might lead to an increased number of false positives at lower ranks, the final risk of over-classification is lower due to high quality of the DAIRYdb and the comparison with Metaxa2.
We suggest to use the Template.taxonomy.xlsx file for final taxonomic classification using the results from both tools. With the DAIRYdb and this approach, about 90% of all OTUs from dairy samples should obtain a confident species annotation.

Alternatively use the python script crossvalid_tax (https://github.com/erifa1/crossvalid_tax.git) for an automated process, which does the same as the Excel file. A cross-check with the Excel file might give you more control over the final annotation.

## Note
DAIRYdb is under active development and validation. Please independently confirm the DAIRYdb predictions by manually inspecting the tree and bringing any discrepancies to our attention.
Moreover, please let us know if you want DAIRYdb to be adapted to a specific classifier not yet available here. Also, if you adapted DAIRYdb to any other classifier, do not hesitate to send us the files so that we can push them on github.


## Licence
 <img src="https://raw.githubusercontent.com/marcomeola/DAIRYdb/master/images/etalab.png" width="50"> [ETALAB](https://www.etalab.gouv.fr/wp-content/uploads/2017/04/ETALAB-Licence-Ouverte-v2.0.pdf)
GPL 3.0
## Copyright
2019 Agroscope, INRA

## Disclaimer
DAIRYdb is released under the ETALAB and GPL 3.0 licenses. The software is therefore open-source and free to use, as long as any modification to the source code will be exclusively for your sole purpose, or released within the terms of the license. Any commercial sale (standalone or as part of a package) is forbidden. DAIRYdb is made available to the community is delivered without any warranty, as expressed by the terms of this disclaimer. It is implied that you agree with the terms of the license and the disclaimer, if you decide to use the DAIRYdb.


## Citation
If you use the DAIRYdb, please cite:

Marco Meola, Etienne Rifa, Noam Shani, Céline Delbes, Hélène Berthoud, Christophe Chassard. (2019) BMC Genomics 20(1):560. https://doi.org/10.1186/s12864-019-5914-8

## References
If you use the DAIRYdb implemented with one of the mentioned classification tools, please cite accordingly:

#### SINTAX
Edgar, R.: SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences. bioRxiv, 074161(2016). doi:10.1101/074161

#### Metaxa2
Bengtsson-Palme, J., Hartmann, M., Eriksson, K.M., Pal, C., Thorell, K., Larsson, D.G.J., Nilsson, R.H.: Metaxa2: improved identification and taxonomic classification of small and large subunit rrna in metagenomic data. Mol Ecol Resour, 15(6), 1403–14 (2015). doi:10.1111/1755-0998.12399

Bengtsson-Palme, J., Richardson, R.T., Meola, M., Wurzbacher, C., Tremblay, E.D., Thorell, K., Kanger, K., Eriksson, K.M., Bilodeau, G.J., Johnson, R.M., Hartmann, M., Henrik Nilsson, R.: Metaxa2 database builder: Enabling taxonomic identification from metagenomic or metabarcoding data using any genetic marker. Bioinformatics, 482 (2018). doi:10.1093/bioinformatics/bty482

#### Blast+
Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L.: BLAST+: architecture and applications. BMC Bioinformatics, 10, 421 (2009). doi:10.1186/1471-2105-10-421

#### Qiime2
Bokulich, N.A., Kaehler, B.D., Rideout, J.R., Dillon, M., Bolyen, E., Knight, R., Huttley, G.A. and Caporaso, J.G.: Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome, 6(1), 90 (2018). doi:10.1186/s40168-018-0470-z

#### FROGS
Escudié, F., Auer, L., Bernard, M., Mariadassou, M., Cauquil, L., Vidal, K., Maman, S., Hernandez-Raquet, G., Combes, S., Pascal, G.: FROGS: Find, Rapidly, OTUs with Galaxy Solution, Bioinformatics, 34(8), 1287–1294 (2018). doi: 10.1093/bioinformatics/btx791

#### IDTAXA
Murali, A., Bhargava, A., Wright, E. S.: IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences, Microbiome, 6:140, (2018). doi: 10.1186/s40168-018-0521-5

#### Kraken2
Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). doi: 10.1186/s13059-019-1891-0


## Archive

Previous versions of the DAIRYdb are available at [this link.](https://frama.link/dairydb_archive)

## Contact
marco.meola@agroscope.admin.ch
