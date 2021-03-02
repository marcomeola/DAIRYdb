# Create DB

```bash
# search taxid base on taxon names1
sed 1d ../DAIRYdb_v1.2.4_20200221_OK.csv|cut -f16|cut -f1 -d_ > taxon_list.txt

# short ids for blastdb (50characters max)
sed 1d ../DAIRYdb_v1.2.4_20200221_OK.csv|cut -f1|awk '{print "seq" $1}' > newids.txt
sed 1d ../DAIRYdb_v1.2.4_20200221_OK.csv|awk -F"\t" '{print ">seq" $1 "\n" $21}' > DDB_blast_nanoclust.fasta

#taxonkit to determine taxids
taxonkit name2taxid -j 6 -i 1 --data-dir /home-local/rifa/bank/blast/taxdump/ taxon_list.txt > taxid.txt

# sed -i 1d taxid.txt

# Taxid_map file
cut -f2 taxid.txt > taxid2.txt
paste newids.txt taxid2.txt > taxid_map.txt

# Some missing taxids were fixed to "2" (Bacteria)

# Makeblast db
~/ssd/softs/ncbi-blast-2.10.1+/bin/makeblastdb -in DDB_blast_nanoclust.fasta -parse_seqids -blastdb_version 5 -taxid_map ../taxid_map.txt -title "ddb_test" -dbtype nucl
```


# test NanoCLUST
```bash
cd ~/ssd/projets/full4best/dumas/nanopore/test_nanoclust_ddb

nextflow run /home-local2/rifa/repository/NanoCLUST/main.nf \
            -profile conda \
            --min_read_length 1000 \
            --reads '../reads/reads_ok_5k/*5k.fastq' \
            --db "db/ddb/DDB_blast_nanoclust.fasta" \
            --tax "db/taxdb/" \
            -resume

```
