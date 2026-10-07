#!/bin/bash
# fetch_1kg_trio.sh KID DAD MOM KID_SEX DAD_SEX MOM_SEX OUT_DIR [triokaryo run args...]
#
# A trio of the 1000 Genomes Project, from Illumina's public DRAGEN 3.7.6 re-analysis of the 3,202 high-coverage genomes
# (AWS Open Data, bucket 1000genomes-dragen, data/dragen-3.7.6/hg38-graph-based/<sample>/): the three per-sample
# hard-filtered small-variant VCFs (about 430 MB each, with AD, DP and GQ), the per-sample CNV and ploidy calls (for the
# comparison), reduced to PASS biallelic SNVs, merged into one trio VCF with `bcftools merge -0` (a sample without a
# record at a site is homozygous reference: triokaryo takes a parent's such call as confident), and run.
#
#   docs/example/fetch_1kg_trio.sh NA12739 NA12748 NA12749 1 1 2 out/NA12739 --panel 1kg-dragen --events docs/example/events.external.tsv
#
# Needs bcftools and a python with triokaryo (TK_PY, default python3). The downloads land in OUT_DIR/download (1.3 GB per
# trio): delete it when done. Pass --panel 1kg-dragen (the shipped panel of twelve of these genomes) or a panel built with
# `triokaryo panel`: without one, the reference's own depth structure reads as events.
set -euo pipefail
kid=$1; dad=$2; mom=$3; ks=$4; ds=$5; ms=$6; out=$7; shift 7
B="https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.7.6/hg38-graph-based"
TK=${TK_PY:-python3}
d="$out/download"; mkdir -p "$d"
for s in "$kid" "$dad" "$mom"; do
    for f in hard-filtered.vcf.gz hard-filtered.vcf.gz.tbi cnv.vcf.gz ploidy.vcf.gz; do
        [ -s "$d/$s.$f" ] || curl -sS --fail --retry 3 -o "$d/$s.$f" "$B/$s/$s.$f"
    done
    [ -s "$d/$s.snv.vcf.gz" ] || { bcftools view -f PASS -v snps -m2 -M2 -Oz -o "$d/$s.snv.vcf.gz" "$d/$s.hard-filtered.vcf.gz" && bcftools index -t "$d/$s.snv.vcf.gz"; }
done
[ -s "$d/trio.vcf.gz" ] || { bcftools merge -0 -m none -Oz -o "$d/trio.vcf.gz" "$d/$kid.snv.vcf.gz" "$d/$dad.snv.vcf.gz" "$d/$mom.snv.vcf.gz" && bcftools index -t "$d/trio.vcf.gz"; }
printf '#kid\tdad\tmom\tkid_sex\tdad_sex\tmom_sex\n%s\t%s\t%s\t%s\t%s\t%s\n' "$kid" "$dad" "$mom" "$ks" "$ds" "$ms" > "$d/trios.tsv"
"$TK" -m triokaryo.cli run --vcf "$d/trio.vcf.gz" --pedigree "$d/trios.tsv" --child "$kid" --out "$out" "$@"
echo "done: $out/index.html  (downloads in $d: delete when finished)"
