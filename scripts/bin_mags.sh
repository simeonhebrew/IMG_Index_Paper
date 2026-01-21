#Paths
READ_DIR="/scratch/nthukus/SEQUENCING_READS"
OUTPUT_DIR="/scratch/nthukus/bins"


cd $READ_DIR

for R1 in *_1.fastq
do
R2=${R1//_1.fastq/_2.fastq}
bn=$R1
filename="${bn%_*.*}"
megahit -1 $R1 -2 $R2 -t 32 -o ${filename}_megahit --out-prefix ${filename}
done

mkdir CONTIG_DIR

mv *_megahit/*contigs.fa CONTIG_DIR/


# Create output directory
mkdir -p "$OUTPUT_DIR"

#Metagenomic binning through coverage depth estimation
for contigs in "CONTIG_DIR"/*.contigs.fa; do
    sample=$(basename "$contigs" .contigs.fa)

    r1="$READ_DIR/${sample}_1.fastq.gz"
    r2="$READ_DIR/${sample}_2.fastq.gz"

    echo "Processing $sample"

    # Build Bowtie2 index
    bowtie2-build "$contigs" "${sample}_bt2_index"

    # Align reads to contigs
    bowtie2 -x "${sample}_bt2_index" -1 "$r1" -2 "$r2" -p 32 | samtools sort -@ 4 -o "${sample}.bam"

    samtools index "${sample}.bam"

    # Compute coverage depth
    jgi_summarize_bam_contig_depths --outputDepth "${sample}.depth.txt" "${sample}.bam"

    # Run MetaBAT2
    mkdir -p "$OUTPUT_DIR/$sample"
    metabat2 -i "$contigs" -a "${sample}.depth.txt" -o "$OUTPUT_DIR/$sample/bin" --minContig 1500 -t 32

    echo "Finished $sample"
done
