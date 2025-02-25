## Genome analysis with Blast
Here we show how we used blast to search for orthologous genes (Table 2) in the genome of Prototheca ciferrii (Accession number: GCA_003613005.1)
```
#conda activate EDTA

#prepare db
makeblastdb -in GCA_003613005.1_ASM361300v1_genomic.fna -dbtype nucl -out prototheca_db

#database cointainig the Acession numbers: 
XM_009314279.1
AJ_866271.1
XM_002181496.1
XM_002180472.1
XM_002286473.1
X_77132.1
XM_001385419.1
NM_001020457.1

#perform blast
run-blast (){

    fasta=$1
        blastn -query $fasta \
               -db /RAID/Data/databases/nt/nt \
               -outfmt "6 qseqid staxids bitscore std sscinames scomnames" \
               -max_hsps 1 \
               -evalue 1e-25 \
               -num_threads 20 \
               -out ${fasta}.blast.out
}

run-blast GCA_003613005.1_ASM361300v1_genomic.fna
blastn -query genes_of_interest2.fasta -db prototheca_db -out results.out -outfmt 6 -num_threads 32

#Results
XM_009314279.1	PEIA01006889.1	71.637	684	170	24	253	924	5071	5742	1.30e-40	167
XM_009314279.1	PEIA01006686.1	71.736	651	158	24	288	924	3682	4320	2.17e-38	159

```
