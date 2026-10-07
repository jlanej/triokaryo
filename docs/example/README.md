# Example: two 1000 Genomes trios

All data here derive from public genomes: Illumina's DRAGEN 3.7.6 re-analysis of the 1000 Genomes Project's 3,202
high-coverage genomes (AWS Open Data, bucket `1000genomes-dragen`, `data/dragen-3.7.6/hg38-graph-based/<sample>/`).
The per-sample `hard-filtered.vcf.gz` files (with `AD`, `DP`, `GQ`) were reduced to PASS biallelic SNVs and merged
into one VCF per trio. The comparison calls are NGS-DOSE's alignment-free karyotype of the same cohort and DRAGEN's
own CNV calls. Sample IDs are the public ones. The landing page [`../index.html`](../index.html) shows the figures;
this file lists the results.

Both runs used the shipped reference panel (`--panel 1kg-dragen`) and no GC track (`gc_corrected = 0` in the
summaries). Run time was 60 s (NA12739, 5.4 million sites) and 72 s (HG01103, 5.9 million sites). The meiotic-stage
and sex-chromosome columns added later were computed from the committed phased windows and bins (`triokaryo report`
rebuilds the pages from the tables); the father/son Y depth ratio needs site-level depths and is not available for
these runs.

## Results

| trio (child) | NGS-DOSE | triokaryo |
| --- | --- | --- |
| CEU 1444, **NA12739** (father NA12748, mother NA12749) | `47,XY,+12` (2.99 copies, every cell) | gain of the whole of chromosome 12, f = 0.96 by depth, 0.94 by the folded bands, 0.95 by the phased track. Extra copy **paternal** (17,630 phased sites; the likelihood ratio over 17,678 informative sites agrees). The two paternal copies are **one homologue throughout** (hetero_share 0.004), isodisomic at the centromere with no crossover: stage "mitotic, or meiosis II without a crossover", not a meiosis I nondisjunction. De novo. DRAGEN: CN 3 in thirteen segments |
| the father NA12748 | `46,XY` | gain of 13q from 88 Mb (87.62 Mb at site resolution) to the telomere, f = 0.69 by depth, 0.67 by the folded bands, 0.56 by the phased track, on the homologue transmitted to the child, who did not inherit the gain. DRAGEN: CN 3 in four segments. NGS-DOSE reported 1.999 copies of chromosome 13 |
| the father NA12748 | not reported | **loss of the terminal 21q, 42.8 Mb (42.64 Mb at site resolution) to the telomere at 46.7 Mb, f = 0.59, from the phased scan alone**: four bins, under the depth segmentation's five-bin floor, although the depth deviates in the same direction (−55 SE). On the transmitted homologue, which the child received intact. The segment's Mendelian-error rate (0.077 against 0.010 genome-wide) is consistent with the father's heterozygous sites being called homozygous where the lost allele's read fraction falls to 0.29. DRAGEN: CN 1 over 44.2–46.7 Mb (segment mean 0.68) |
| NA12739 | not reported | a 12% gain of chr1:0–5 Mb by depth that the phased track puts at 4%: flagged as a possible depth artefact |
| PR 19, **HG01103** (father HG01101, mother HG01102) | `46,XX,-2(140-179Mb)[0.68]` | loss of 2q from 141 to 173 Mb (140.68–173.16 Mb at site resolution), f = 0.67 by depth, 0.62 by the folded bands, 0.57 by the phased track, followed by 173–180 Mb at f = 0.31 (a subclone); the **paternal** copy is lost (4,060 phased sites). De novo. DRAGEN: CN 1 in five segments over 140.7–179.5 Mb, and a sixth of 2.1 Mb at 200.4–202.4 Mb, below triokaryo's size floor |
| HG01103 | not reported | loss of 14q 35–41 Mb, f = 0.64 (0.62 by both allelic estimates), paternal copy. DRAGEN: CN 1 in two segments over 35.3–41.2 Mb |
| HG01103, HG01101, HG01102 | | five homozygous segments of 5 to 15 Mb (the child's 10p; the father's 11p, pericentromeric 11 and 17q; the mother's 7q) classified as **runs of homozygosity**: heterozygosity rate ≤ 0.09 of the member's own, no excess of Mendelian errors, few phased sites. Not uniparental disomies |
| HG01101 | | two gains by depth (9q34 f = 0.10, 19p13 f = 0.14) that the phased track puts at 0.01–0.03: flagged |
| all six members | | sex-chromosome complements XY, XY, XX (NA12739 trio) and XX, XY, XX (HG01103 trio), all agreeing with the pedigree; the males' Y reads 0.98–1.02 copies against the panel's male Y level and the females' 0; no whole-X or whole-Y event (X copy numbers 1.00, 1.00, 1.94 and 1.98, 1.00, 1.97) |

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
cluster), the 22q11 low-copy repeats, chr15:31.5–34.4 Mb, chr3:192–198 Mb, chr8:8–12 Mb and chr10:131–134 Mb, where a
dense cluster of sites or paralogous sequence shifts the pooled fraction over a short span or in a few windows only;
and a bands-only shift of 0.019 on the mother's X (chrX:145–156 Mb), under the 0.025 floor. An earlier version of the
method, which tested the band deviation before the heterozygosity rate, had called a 42% CN-LOH of 7q in HG01102 and
mosaic LOH segments in HG01101; with the heterozygosity rate tested first they are constitutional runs of
homozygosity, as the phased fraction (at 1/2, from few sites) confirms.

## The reference panel

Each genome carries the reference's depth structure (centromere flanks, segmental duplications, acrocentric short
arms), and paralogous sequence shifts the BAF bands in every genome. Normalised against its own median alone, each
member shows dozens of spurious events there. The panel subtracts the median LRR, band deviation and heterozygosity
rate per 1-Mb bin over other genomes processed the same way. The shipped panel (`--panel 1kg-dragen`,
`panel.1kg_dragen_3.7.6.1mb.tsv`) comprises twelve such genomes: six unrelated (HG00096, HG00097, HG00099, HG00100,
HG00101, HG00102) and the six members of the two trios shown here. The example is therefore not a fully independent
evaluation: with 12 genomes the per-bin median is insensitive to an event present in one or two of them, but the
trios' own genomes contributed to their reference. Bins the panel cannot characterise (fewer than five genomes, robust
SD above 0.25) are masked; bins whose bands are shifted across the panel are excluded from the CN-LOH search and the
phased tracks. For another pipeline (GATK's `DP`, another depth filter) build a panel from that cohort with
`triokaryo panel --vcfs ...` or, with many trios, from a first pass's `bins.tsv` files.

## Reproduce

```bash
# one trio, end to end (downloads 1.3 GB per trio into OUT/download; delete it afterwards)
docs/example/fetch_1kg_trio.sh NA12739 NA12748 NA12749 1 1 2 out/NA12739 --panel 1kg-dragen --events docs/example/events.external.tsv
docs/example/fetch_1kg_trio.sh HG01103 HG01101 HG01102 2 1 2 out/HG01103 --panel 1kg-dragen --events docs/example/events.external.tsv
triokaryo cohort --runs out/NA12739 out/HG01103 --events docs/example/events.external.tsv --out out/cohort
```

`events.external.tsv` holds NGS-DOSE's karyotype events for the two trios and DRAGEN's CNV calls of 2 Mb and more
(PASS `DRAGEN:GAIN`/`LOSS` records from each sample's `cnv.vcf.gz`, with the integer copy number and the segment mean)
in the generic `sample chrom start end label type` form. The panel was built with
`triokaryo panel --vcfs <the twelve hard-filtered VCFs> --thin 2`.
