#!/bin/bash
# Whole-chromosome aneuploidies across the 1000 Genomes DRAGEN cohort, and which of them sit in a trio and have an assembly.
# Downloads the public inputs (about 6 MB: one small ploidy VCF per sample, the pedigree, the assembly indexes) into WORK
# and runs scan.py. Usage: fetch.sh WORK_DIR [OUT_DIR]
set -euo pipefail
W=${1:?work dir}; OUT=${2:-$(dirname "$0")}
HERE="$(cd "$(dirname "$0")" && pwd)"
B="https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.7.6/hg38-graph-based"
mkdir -p "$W/ploidy"; cd "$W"
[ -s 3202.ped ] || curl -sS --fail -o 3202.ped "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
[ -s hprc_y1.index ] || curl -sS --fail -o hprc_y1.index "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index"
[ -s hprc_r2.csv ] || curl -sS --fail -L -o hprc_r2.csv "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_release2_v1.0.index.csv"
[ -s hgsvc3_accessions.tsv ] || curl -sS --fail -L -o hgsvc3_accessions.tsv "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Assembly_Info/v2.0/supp-table_accessions-assembly-data.2024v2.0.tsv"
# one ploidy VCF per sample (about 2 kB each); a few samples are absent from the bucket and are reported
awk 'NR>1{print $2}' 3202.ped | xargs -P 16 -I{} sh -c '[ -s "ploidy/{}.ploidy.vcf.gz" ] || curl -sS --fail --retry 3 -o "ploidy/{}.ploidy.vcf.gz" "'"$B"'/{}/{}.ploidy.vcf.gz" || echo "absent: {}"'
python3 "$HERE/scan.py" --ploidy-dir ploidy --ped 3202.ped --hprc-y1 hprc_y1.index --hprc-r2 hprc_r2.csv --hgsvc3 hgsvc3_accessions.tsv --out "$OUT"
