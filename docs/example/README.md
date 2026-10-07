# Example: two 1000 Genomes trios

All data here derive from public genomes: Illumina's DRAGEN 3.7.6 re-analysis of the 1000 Genomes Project's 3,202
high-coverage genomes (AWS Open Data, bucket `1000genomes-dragen`, `data/dragen-3.7.6/hg38-graph-based/<sample>/`).
The per-sample `hard-filtered.vcf.gz` files (with `AD`, `DP`, `GQ`) were reduced to PASS biallelic SNVs and merged
into one VCF per trio. The comparison calls are NGS-DOSE's alignment-free karyotype of the same cohort and DRAGEN's
own CNV calls. Sample IDs are the public ones. The landing page [`../index.html`](../index.html) shows the figures;
this file lists the results.

Both runs used the shipped reference panel (`--panel 1kg-dragen`) and the shipped GC track (the default at 1 Mb on
GRCh38; `gc_corrected = 1` in the summaries), with the merged trio VCFs as the only input. Run time was 88 s (NA12739,
5.4 million sites) and 100 s (HG01103, 5.9 million sites) after the merge.

## Results

| trio (child) | NGS-DOSE | triokaryo |
| --- | --- | --- |
| CEU 1444, **NA12739** (father NA12748, mother NA12749) | `47,XY,+12` (2.99 copies, every cell) | gain of the whole of chromosome 12, f = 0.98 by depth, 0.94 by the folded bands, 0.95 by the phased track. Extra copy **paternal** (17,630 phased sites; the likelihood ratio over 17,678 informative sites agrees). The two paternal copies are **one homologue throughout** (hetero_share 0.004), isodisomic at the centromere with no crossover: stage "mitotic, or meiosis II without a crossover", not a meiosis I nondisjunction. De novo. DRAGEN: CN 3 in thirteen segments |
| the father NA12748 | `46,XY` | gain of 13q from 88 Mb (87.62 Mb at site resolution) to the telomere, f = 0.68 by depth, 0.67 by the folded bands, 0.63 by the phased track. The son's paternal chromosome 13 carries a **crossover at 112.7 Mb**: up to it the duplicated homologue is the one the son received, beyond it the other, and the transmitted-allele fraction flips sign there (`crossovers.tsv`); read with the sign, the phased estimate had been 0.56. The son did not inherit the gain. DRAGEN: CN 3 in four segments. NGS-DOSE reported 1.999 copies of chromosome 13 |
| the father NA12748 | not reported | **loss of the terminal 21q, 42.8 Mb (42.64 Mb at site resolution) to the telomere at 46.7 Mb, f = 0.60, from the phased scan alone**: four bins, under the depth segmentation's five-bin floor, although the depth deviates in the same direction (−93 SE). On the transmitted homologue, which the son received intact. The segment's Mendelian-error rate (0.077 against 0.010 genome-wide) is consistent with the father's heterozygous sites being called homozygous where the lost allele's read fraction falls to 0.29. DRAGEN: CN 1 over 44.2–46.7 Mb (segment mean 0.68) |
| the mother NA12749 | `46,XX` | a loss of Xq28 (147.9–155.7 Mb) at f = 0.09 from the phased scan, the depth leaning the same way (−8 SE). Her X reads 1.97 copies and the folded phased shift runs at about 0.012 along the whole X, flipping sign at the son's maternal crossovers: a loss of one X in about 5% of her cells, below the depth's 10% floor and at the scan's 6% one, which reports the one block that exceeds it. That block lies over the Xq28 segmental duplications, where the other trio's mother shows a shift of 0.019 too (rejected there, the depth flat), so the reported block is in part the region's own structure |
| PR 19, **HG01103** (father HG01101, mother HG01102) | `46,XX,-2(140-179Mb)[0.68]` | loss of 2q from 141 to 173 Mb (140.68–173.16 Mb at site resolution), f = 0.67 by depth, 0.61 by the folded bands, 0.57 by the phased track, followed by 173–180 Mb at f = 0.32 (a subclone); the **paternal** copy is lost (4,060 phased sites). De novo. DRAGEN: CN 1 in five segments over 140.7–179.5 Mb, and a sixth of 2.1 Mb at 200.4–202.4 Mb, below triokaryo's size floor |
| HG01103 | not reported | loss of 14q 35–41 Mb, f = 0.62 by all three estimates, paternal copy. DRAGEN: CN 1 in two segments over 35.3–41.2 Mb |
| HG01103, HG01101, HG01102 | | five homozygous segments of 5 to 15 Mb (the child's 10p; the father's 11p, pericentromeric 11 and 17q; the mother's 7q) classified as **runs of homozygosity**: heterozygosity rate ≤ 0.09 of the member's own, no excess of Mendelian errors, few phased sites. Not uniparental disomies |
| all six members | | nothing else: no depth call below 10% anywhere (an earlier version, without the GC track, the tighter panel mask and the spike smoothing, had flagged gains of 10–25% at 1p36, 9q34, 19p13, 22q11 and the pericentromeres of 9, 10 and X, which the phased track put at 1–6%) |
| the two children | | karyotype strings `47,XY,+12pat(MII/mit)` (NA12739) and `mos 46,XX,del(2)(q22.1q31.1)pat[0.67]/46,XX,del(2)(q31.1q31.2)pat[0.32]/46,XX,del(14)(q13.2q21.1)pat[0.62]/46,XX,roh(10)(p12.2p11.22)` (HG01103); the parents' `mos 46,XY,dup(13)(q31.2q34)[0.68]/46,XY,del(21)(q22.3)[0.60]/46,XY` and `mos 46,XX,del(X)(q28)[0.09]/46,XX`; the bands come from the site-resolution boundaries |
| all six members | | sex-chromosome complements XY, XY, XX (NA12739 trio) and XX, XY, XX (HG01103 trio), all agreeing with the pedigree; the males' Y reads 1.00–1.01 copies against the panel's male Y level (father/son Y depth ratio −0.018 log2 over 4,190 sites in NA12739) and the females' 0; no whole-X or whole-Y event (X copy numbers 1.01, 1.01, 1.97 and 2.01, 1.01, 2.00) |

The pages: [`NA12739/index.html`](NA12739/index.html), [`HG01103/index.html`](HG01103/index.html), and the cohort
report with the concordance against the supplied events, [`cohort/index.html`](cohort/index.html): 27 of the 28
supplied events are matched; the exception is the 2.1-Mb DRAGEN segment. Each page embeds its figures; the PNGs,
their sidecars (title, caption, key) and legend images are under `figures/`.

In the figures: on chromosome 12 of NA12739 the maternal-allele fraction sits at 1/3 along the whole chromosome and
the copy-number track reads paternal 2, maternal 1. On chromosome 21 the father's transmitted-allele fraction steps
from 1/2 to 0.29 at 42.8 Mb while the transmitted homologue's copy number falls to 0.4 and the other stays at 1. On
chromosome 2 of HG01103 the maternal fraction steps to 0.70 over 141–173 Mb and the paternal copy number falls to
0.35. In the genome-wide figure the son's X reads maternal 1, paternal 0.

Segments the phased scan rejected (`phased_rejected.tsv`): the MHC (chr6:31.6–32.7 Mb), chr19:42.8–43.0 Mb (the PSG
cluster), chr22:16.1–16.4 Mb (the 22q11 low-copy repeats), chr15:31.5–34.4 Mb, chr3:192–198 Mb and chr8:8–12 Mb,
where a dense cluster of sites or paralogous sequence shifts the pooled fraction over a short span or in a few windows
only; and bands-only shifts on the mothers' X under the 0.025 floor (NA12749 over 9–25 Mb, 0.017, a block of the
whole-X signal above; HG01102 over 146–156 Mb, 0.019, the Xq28 structure). An earlier version of the method, which
tested the band deviation before the heterozygosity rate, had called a 42% CN-LOH of 7q in HG01102 and mosaic LOH
segments in HG01101; with the heterozygosity rate tested first they are constitutional runs of homozygosity, as the
phased fraction (at 1/2, from few sites) confirms.

## The reference panel

Each genome carries the reference's depth structure (centromere flanks, segmental duplications, acrocentric short
arms), and paralogous sequence shifts the BAF bands in every genome. Normalised against its own median alone, each
member shows dozens of spurious events there. The panel subtracts the median LRR, band deviation and heterozygosity
rate per 1-Mb bin over other genomes processed the same way. The shipped panel (`--panel 1kg-dragen`,
`panel.1kg_dragen_3.7.6.1mb.tsv`) comprises twelve such genomes: six unrelated (HG00096, HG00097, HG00099, HG00100,
HG00101, HG00102) and the six members of the two trios shown here. The example is therefore not a fully independent
evaluation: with 12 genomes the per-bin median is insensitive to an event present in one or two of them, but the
trios' own genomes contributed to their reference. Bins the panel cannot characterise (fewer than five genomes, or a
robust SD above 0.10, where the panel's genomes disagree by more than the smallest reportable event: about 3% of the
bins, at pericentromeres, 22q11 and Xq28) are masked; bins whose bands are shifted across the panel are excluded from
the CN-LOH search and the phased tracks. For another pipeline (GATK's `DP`, another depth filter) build a panel from that cohort with
`triokaryo panel --vcfs ...` or, with many trios, from a first pass's `bins.tsv` files.

## Reproduce

```bash
# one trio, end to end (downloads 1.3 GB per trio into OUT/download; delete it afterwards)
docs/example/fetch_1kg_trio.sh NA12739 NA12748 NA12749 1 1 2 out/NA12739 --panel 1kg-dragen --events docs/example/events.external.tsv
docs/example/fetch_1kg_trio.sh HG01103 HG01101 HG01102 2 1 2 out/HG01103 --panel 1kg-dragen --events docs/example/events.external.tsv
triokaryo cohort --runs out/NA12739 out/HG01103 --events docs/example/events.external.tsv --out out/cohort
```

The script's bcftools steps (PASS biallelic SNVs per sample, `merge -0`, index) are what `triokaryo merge --child --father
--mother --out` does in one command.

`events.external.tsv` holds NGS-DOSE's karyotype events for the two trios and DRAGEN's CNV calls of 2 Mb and more
(PASS `DRAGEN:GAIN`/`LOSS` records from each sample's `cnv.vcf.gz`, with the integer copy number and the segment mean)
in the generic `sample chrom start end label type` form. The panel was built with
`triokaryo panel --vcfs <the twelve hard-filtered VCFs> --thin 2`, and the shipped GC track with
`triokaryo gc-track --fasta hg38.fa` from the UCSC hg38 sequence.
