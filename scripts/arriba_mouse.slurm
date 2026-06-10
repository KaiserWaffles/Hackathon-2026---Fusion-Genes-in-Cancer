#!/bin/bash
#SBATCH --job-name=arriba_PRJNA1062304
#SBATCH --partition=compute-low-priority
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=$HOME/shared_data/bbermude/hackathon_2026/logs/%x.%j.out
#SBATCH --error=$HOME/shared_data/bbermude/hackathon_2026/logs/%x.%j.err

# IMPORTANT: replace $HOME for logs with your whole home path.

hostname

date

module load STAR
module load samtools
module load R-4.5.3

A_DIR="$HOME/shared_data/bbermude/arriba_v2.5.1"
WD=$HOME/shared_data/bbermude/hackathon_2026/PRJNA1062304
IN_DIR="$WD/fastq"
OUT_DIR="$WD/arriba"
mkdir -p $OUT_DIR

# Change to working directory
cd $OUT_DIR
for R1 in "$IN_DIR"/*_R1.fastq.gz
do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    echo "Sample: $SAMPLE"
    echo "  R1: $R1"
    echo "  R2: $R2"

    # Example:
    # bwa mem ref.fa "$R1" "$R2" > "${SAMPLE}.sam"
    $A_DIR/run_arriba.sh $A_DIR/STAR_index_GRCm39_GENCODEM38/ $A_DIR/GENCODEM38.gtf $A_DIR/GRCm39.fa $A_DIR/database/blacklist_mm39_GRCm39_v2.5.1.tsv.gz $A_DIR/database/known_fusions_mm39_GRCm39_v2.5.1.tsv.gz $A_DIR/database/protein_domains_mm39_GRCm39_v2.5.1.gff3 16 $R1 $R2

    # now rename output files
    mv Aligned.sortedByCoord.out.bam $SAMPLE.Aligned.sortedByCoord.out.bam
    mv Aligned.sortedByCoord.out.bam.bai $SAMPLE.Aligned.sortedByCoord.out.bam.bai
    mv fusions.tsv $SAMPLE.fusions.tsv
    mv fusions.discarded.tsv $SAMPLE.fusions.discarded.tsv
    mv SJ.out.tab $SAMPLE.SJ.out.tab
    mv Log.out $SAMPLE.Log.out

    # now make arriba visualization fusion.pdf
    $A_DIR/draw_fusions.R --fusions=$SAMPLE.fusions.tsv --alignments=$SAMPLE.Aligned.sortedByCoord.out.bam --output=$SAMPLE.fusions.pdf --annotation=$A_DIR/GENCODEM38.gtf --cytobands=$A_DIR/database/cytobands_mm39_GRCm39_v2.5.1.tsv --proteinDomains=$A_DIR/database/protein_domains_mm39_GRCm39_v2.5.1.gff3

done

date
