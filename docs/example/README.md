# Example: two 1000 Genomes trios, from public data

Everything here was made from public data: Illumina's DRAGEN 3.7.6 re-analysis of the 1000 Genomes Project's 3,202
high-coverage genomes (AWS Open Data, bucket `1000genomes-dragen`, `data/dragen-3.7.6/hg38-graph-based/<sample>/`),
read as the per-sample `hard-filtered.vcf.gz` (with `AD`, `DP`, `GQ`), and NGS-DOSE's alignment-free karyotype of the
same cohort (NGS-DOSE-1000G, `docs/data/cohort.tsv`). The sample IDs are the public ones.

| trio (child) | what NGS-DOSE read | what triokaryo reads from the VCF |
| --- | --- | --- |
| CEU 1444, **NA12739** (father NA12748, mother NA12749) | `47,XY,+12` (2.99 copies, every cell) | a gain of the whole of chromosome 12 in 96% of cells by depth and 94% by the bands; the extra copy is **paternal** (log-likelihood ratio 1.8 × 10⁵ over 17,678 informative sites): one homologue of the father's doubled, as a line's trisomy is; new, neither parent carries it |
| the father NA12748 | `46,XY`, nothing | a gain of 13q from 88 Mb to the end in 69% of cells by depth and 67% by the bands; DRAGEN's own CNV caller reads CN 3 there too; not passed to the child |
| PR 19, **HG01103** (father HG01101, mother HG01102) | `46,XX,-2(140-179Mb)[0.68]` | a loss of 2q from 141 to 173 Mb in 67% of cells by depth and 62% by the bands, then 173 to 180 Mb in 31%; the **paternal** copy is lost (ratio 3.9 × 10⁴ over 4,060 sites); DRAGEN reads CN 1 over 140.7 to 172 Mb; new |
| HG01103 | nothing on 14 | a loss of 14q from 35 to 41 Mb in 64% of cells, the paternal copy; DRAGEN reads CN 1 over 37.3 to 41.2 Mb |

The pages: [`NA12739/index.html`](NA12739/index.html), [`HG01103/index.html`](HG01103/index.html), and the two gathered
with the concordance against the given events in [`cohort/index.html`](cohort/index.html). Each page embeds its figures;
the PNGs, their sidecars (title, caption, key) and legend images are beside them under `figures/`.

What the figures show that the numbers summarise: on chromosome 12 of NA12739 the heterozygous sites sit on two bands at
1/3 and 2/3, and the sites where the parents are opposite homozygotes split by the parent of the alt allele - the
father's alleles on the 2/3 band, the mother's on 1/3. On chromosome 2 of HG01103 the bands sit at 0.24 and 0.76 from
141 to 173 Mb with the mother's alleles on the upper band (the paternal copy is the one lost), and the depth and the
bands give the same share of cells. The within-family track (the child over the parents' mean) shows the child's
chromosome 12 and, mirrored, the father's 13q.

Beyond the events in the table, the pages list a few low-level calls (a 5-Mb gain of 1p in 12% of NA12739's cells; in
HG01101 a 9% gain of 9q and 14% of 19p and three mosaic copy-neutral LOH stretches; in HG01102 a copy-neutral LOH of
8 Mb of 7q in 42%) and two runs of homozygosity read as "in every cell, no Mendelian errors" (HG01103 on 10p, HG01101
on 11p): identical-by-descent segments, not uniparental disomies. None of these has a second method behind it here.

## The panel

Real genomes carry a depth structure that is the reference's, not theirs - centromere flanks, segmental duplications,
the acrocentric short arms - and paralogous sequence parts the B-allele bands in everyone. Read against its own median
alone, each member of a trio shows dozens of "events" there. The panel takes it out: the median LRR, band deviation and
heterozygosity rate per 1-Mb bin over other genomes called the same way. The one shipped (`--panel 1kg-dragen`,
`panel.1kg_dragen_3.7.6.1mb.tsv`) is twelve of these genomes: six unrelated (HG00096, HG00097, HG00099, HG00100,
HG00101, HG00102) and the two trios. Bins the panel cannot pin (fewer than five genomes, a spread beyond 0.25) are left
out of the calls; bins whose bands are parted in everyone are left out of the LOH search. For another pipeline
(GATK's `DP`, another depth filter) build the panel from that cohort: `triokaryo panel --vcfs ...` or, with many
trios, from a first pass's `bins.tsv` files.

## Reproduce

```bash
# one trio, end to end (downloads 1.3 GB per trio into OUT/download; delete it afterwards)
docs/example/fetch_1kg_trio.sh NA12739 NA12748 NA12749 1 1 2 out/NA12739 --panel 1kg-dragen --events docs/example/events.external.tsv
docs/example/fetch_1kg_trio.sh HG01103 HG01101 HG01102 2 1 2 out/HG01103 --panel 1kg-dragen --events docs/example/events.external.tsv
triokaryo cohort --runs out/NA12739 out/HG01103 --events docs/example/events.external.tsv --out out/cohort
```

`events.external.tsv` holds NGS-DOSE's karyotype events for the two trios and DRAGEN's CNV calls of 3 Mb and more
(from each sample's `cnv.vcf.gz`, `DRAGEN:GAIN`/`LOSS` records with PASS), in the generic `sample chrom start end
label type` form. The panel was built with `triokaryo panel --vcfs <the twelve hard-filtered VCFs> --thin 2`.
