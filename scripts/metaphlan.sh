#Run Metaphlan against the paired end reads
for R1 in *_R1.fastq.gz
do
R2=${R1//_R1.fastq.gz/_R2.fastq.gz}
bn=$R1
filename="${bn%_*.*.*}"
metaphlan $R1,$R2 --bowtie2out METAPHLAN/${filename}.bowtie2.bz2 --bowtie2db /store/nthukus/DATABASES_METAPHLAN/metaphlan_vJune23_DB --nproc 32 --input_type fastq -o METAPHLAN/${filename}.tsv -s METAPHLAN/${filename}.sam.bz2
done


mkdir output

merge_metaphlan_tables.py *.tsv > output/merged_abundance_table.tsv
