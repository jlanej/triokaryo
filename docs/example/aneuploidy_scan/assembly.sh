#!/bin/bash
# How a diploid long-read assembly represents an aneuploid genome's sex chromosomes, held against the trio reading.
# HGSVC3 ships verkko contigs aligned to GRCh38 as BAMs, which samtools reads remotely by chromosome (no download of the
# assembly); HPRC ships the haplotype FASTAs, aligned here to the GRCh38 X and Y with minimap2. Usage: assembly.sh WORK_DIR
set -euo pipefail
W=${1:?work dir}; HERE="$(cd "$(dirname "$0")" && pwd)"; EX="$HERE/.."; mkdir -p "$W"; cd "$W"
H="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230927_verkko_batch2/alignments/contig_to_ref/hg38"
for s in HG03456 HG00731; do                      # the 47,XYY genome and a 46,XY control from the same batch
  mkdir -p "$s"
  for h in hap1 hap2; do
    : > "$s/$h.tsv"
    for c in chrX chrY chr21; do samtools view "$H/$s/$s.vrk-ps-sseq.asm-$h.hg38.sort.bam" $c | python3 "$EX/bam_alignments.py" $h >> "$s/$h.tsv"; done
  done
  python3 "$EX/assembly_sexchrom.py" --chroms chrX,chrY,chr21 --out "$s/summary.tsv" hap1="$s/hap1.tsv" hap2="$s/hap2.tsv"
done
# HG02300: XX by depth in a pedigree male; HPRC release 2, hifiasm trio-binned haplotypes
S="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG02300/assemblies/release2"
for h in pat mat; do [ -s HG02300_$h.fa.gz ] || curl -sS --fail -o HG02300_$h.fa.gz "$S/HG02300_${h}_hprc_r2_v1.0.1.fa.gz"; done
for c in chrX chrY; do [ -s $c.fa.gz ] || curl -sS --fail -o $c.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz"; done
[ -s hg38_XY.fa ] || gzip -dc chrX.fa.gz chrY.fa.gz > hg38_XY.fa
for h in pat mat; do [ -s HG02300_$h.XY.paf ] || minimap2 -t 8 -x asm5 --secondary=no hg38_XY.fa HG02300_$h.fa.gz > HG02300_$h.XY.paf 2> /dev/null; done
python3 "$EX/assembly_sexchrom.py" --chroms chrX,chrY --out HG02300/summary.tsv pat=HG02300_pat.XY.paf mat=HG02300_mat.XY.paf
