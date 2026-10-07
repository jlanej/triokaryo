# Example: two 1000 Genomes trios, from public data

Everything here was made from public data: Illumina's DRAGEN 3.7.6 re-analysis of the 1000 Genomes Project's 3,202
high-coverage genomes (AWS Open Data, bucket `1000genomes-dragen`, `data/dragen-3.7.6/hg38-graph-based/<sample>/`),
read as the per-sample `hard-filtered.vcf.gz` (with `AD`, `DP`, `GQ`) and merged into one trio VCF each, and NGS-DOSE's
alignment-free karyotype of the same cohort (NGS-DOSE-1000G, `docs/data/cohort.tsv`). The sample IDs are the public
ones. The landing page [`../index.html`](../index.html) shows the figures; this file lists the readings.

| trio (child) | what NGS-DOSE read | what triokaryo reads from the VCF |
| --- | --- | --- |
| CEU 1444, **NA12739** (father NA12748, mother NA12749) | `47,XY,+12` (2.99 copies, every cell) | a gain of the whole of chromosome 12 in 96% of cells by depth, 95% by the phased bands; the extra copy is **paternal** (17,630 phased sites; the opposite-homozygote likelihood ratio agrees), and the two paternal copies are **one homologue throughout**: a mitotic duplication or a meiosis II error, not a meiosis I nondisjunction. New: neither parent carries it. DRAGEN: CN 3 in thirteen pieces along the chromosome |
| the father NA12748 | `46,XY`, nothing | a gain of 13q from 88 Mb (87.62 Mb to site resolution) to the end in 69% of cells by depth, 67% by the folded bands, 56% by the phased bands; it lies on the homologue he passed to the child, who did not inherit the gain. DRAGEN: CN 3 in four pieces. NGS-DOSE read 1.999 copies of chromosome 13 |
| the father NA12748 | nothing | **a loss of the end of 21q, 42.8 to 46.1 Mb, in 59% of cells, from the phased bands alone**: four bins long, under the depth scan's five-bin floor, though the depth leans the same way. On the transmitted homologue, which the child received whole. DRAGEN: CN 1 over 44.2 to 46.7 Mb (segment mean 0.68) |
| NA12739 | nothing | a 12% gain of the first 5 Mb of chromosome 1 by depth, which the phased bands put at 4%: flagged in the note as a likely GC artefact |
| PR 19, **HG01103** (father HG01101, mother HG01102) | `46,XX,-2(140-179Mb)[0.68]` | a loss of 2q from 141 to 173 Mb (173.16 Mb to site resolution) in 67% of cells by depth, 62% by the folded bands, 57% by the phased bands, then 173 to 180 Mb in 31%; the **paternal** copy is the one lost (4,060 phased sites). New. DRAGEN: CN 1 in five pieces over 140.7 to 179.5 Mb, and a sixth of 2.1 Mb at 200 Mb under triokaryo's floor |
| HG01103 | nothing on 14 | a loss of 14q from 35 to 41 Mb in 64% of cells, the paternal copy. DRAGEN: CN 1 over 35.3 to 41.2 Mb |
| HG01103, HG01101, HG01102 | | five homozygous stretches of 5 to 15 Mb (the child's 10p; the father's 11p, 11q and 17q; the mother's 7q) read as **runs of homozygosity**: no heterozygous calls, no Mendelian errors, the phased fraction at one half. Not uniparental disomies |
| HG01101 | | two 10–14% gains by depth (9q34, 19p13) that the phased bands put at 1–3%: flagged |

The pages: [`NA12739/index.html`](NA12739/index.html), [`HG01103/index.html`](HG01103/index.html), and the two gathered
with the concordance against the given events in [`cohort/index.html`](cohort/index.html): 28 of the 29 given events
are matched, the one left the 2.1-Mb DRAGEN piece. Each page embeds its figures; the PNGs, their sidecars (title,
caption, key) and legend images are beside them under `figures/`.

What the figures show that the numbers summarise: on chromosome 12 of NA12739 the maternal-allele fraction sits at 1/3
along the whole chromosome and the copies read paternal 2, maternal 1; on chromosome 21 the father's transmitted-allele
fraction steps from one half to 0.29 at 42.8 Mb and the copies show the homologue he passed on dropping to 0.4 while the
other stays at one. On chromosome 2 of HG01103 the maternal fraction steps to 0.70 from 141 to 173 Mb and the paternal
copies fall to 0.35. The genome-wide copies row shows the son's X as maternal 1, paternal 0, as it should.

What the phased bands set aside, with the reason in `phased_rejected.tsv`: the MHC (6p21.3), the PSG cluster (19q13.2),
the 22q11 low-copy repeats, 15q13 and 3q29 - regions where a cluster of sites or paralogous sequence parts the bands in
a member or in everyone. The two trios' first reading had called a 42% copy-neutral LOH of 7q in HG01102 and mosaic
LOH stretches in HG01101: with the heterozygosity rate read first they are runs of homozygosity in every cell, which the
phased fraction (one half, from few sites) confirms.

## The panel

Real genomes carry a depth structure that is the reference's, not theirs - centromere flanks, segmental duplications,
the acrocentric short arms - and paralogous sequence parts the B-allele bands in everyone. Read against its own median
alone, each member of a trio shows dozens of "events" there. The panel takes it out: the median LRR, band deviation and
heterozygosity rate per 1-Mb bin over other genomes called the same way. The one shipped (`--panel 1kg-dragen`,
`panel.1kg_dragen_3.7.6.1mb.tsv`) is twelve of these genomes: six unrelated (HG00096, HG00097, HG00099, HG00100,
HG00101, HG00102) and the two trios. Bins the panel cannot pin (fewer than five genomes, a spread beyond 0.25) are left
out of the calls; bins whose bands are parted in everyone are left out of the LOH search and of the phased tracks. For
another pipeline (GATK's `DP`, another depth filter) build the panel from that cohort: `triokaryo panel --vcfs ...` or,
with many trios, from a first pass's `bins.tsv` files.

## Reproduce

```bash
# one trio, end to end (downloads 1.3 GB per trio into OUT/download; delete it afterwards)
docs/example/fetch_1kg_trio.sh NA12739 NA12748 NA12749 1 1 2 out/NA12739 --panel 1kg-dragen --events docs/example/events.external.tsv
docs/example/fetch_1kg_trio.sh HG01103 HG01101 HG01102 2 1 2 out/HG01103 --panel 1kg-dragen --events docs/example/events.external.tsv
triokaryo cohort --runs out/NA12739 out/HG01103 --events docs/example/events.external.tsv --out out/cohort
```

`events.external.tsv` holds NGS-DOSE's karyotype events for the two trios and DRAGEN's CNV calls of 2 Mb and more
(from each sample's `cnv.vcf.gz`: PASS `DRAGEN:GAIN`/`LOSS` records, with the integer copy number and the segment
mean), in the generic `sample chrom start end label type` form. The panel was built with
`triokaryo panel --vcfs <the twelve hard-filtered VCFs> --thin 2`.
