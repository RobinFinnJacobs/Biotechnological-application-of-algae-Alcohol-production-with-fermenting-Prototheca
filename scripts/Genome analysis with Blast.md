##
Here we show how we used blast to search for orthologous genes (Table 2) in the genome of Prototheca ciferrii 
```
#conda activate EDTA

#prepare db
makeblastdb -in GCA_003613005.1_ASM361300v1_genomic.fna -dbtype nucl -out prototheca_db

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

#Look into GFF
#grep -i "alcohol dehydrogenase" genome.gff


#Use Orthology - Orthofinder
#you need protein sequences (Fasta) from genome assembly
#orthofinder -f input_folder -o output_folder
```
