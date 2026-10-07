# Roadmap

Notes carried forward from the review of 2026-10-07, for the next rounds of work. Items are grouped by theme and
ordered by expected value; "done" items are kept for the record with the commit that closed them.

## Done

- Documentation in concise scientific prose, checked against the code; `triokaryo report` rebuilds captions from the
  tables.
- Inheritance requires reciprocal overlap of at least 50%; supplied events keep the lenient shorter-segment rule.
- A run of homozygosity receives no parent of origin or homologue count from the phased track.

## Method

- **Joint segmentation of LRR and the phased fraction.** Shared breakpoints and one likelihood over (type, cell
  fraction, parent) per segment, replacing the post-hoc typing of phased-scan segments by the depth's deviation. This
  also yields confidence intervals for the cell fraction instead of three point estimates (`f_lrr`, `f_baf`,
  `f_phase`), and a combined estimate with weights.
- **Reference-bias estimator.** The median of alt/depth − 1/2 over heterozygous sites is quantised to ±1/(2·depth) and
  returned exactly 0 for all six example members, so the correction is a no-op at 30×. Replace with a depth-weighted
  mean or a binomial maximum-likelihood estimate.
- **Conditional likelihood ratio.** `origin_llr` assumes the event exists; a flagged depth artefact still receives a
  large ratio. Report it only when the event is supported by the phased track, or add a model-comparison term against
  "no event".
- **Panel-independent control of reference structure.** The within-trio depth tracks cancel shared structure; a
  within-trio segmentation (child over parental mean) could call de novo events without a panel.
- **Join note.** "joined across n bin(s) the panel left out" is written even when no panel masked the gap.

## Calibration and evaluation

- **Sensitivity curves from the simulator.** Sweep cell fraction (2–20%), depth (15–60×), event size (2–20 Mb) and
  site density; report detection per source (depth, bands, phased) and the error of each cell-fraction estimate. The
  mock machinery (`triokaryo mock`) already supports planted events, dense contigs and low cell fractions.
- **False-discovery control for the phased scan.** A tag-permutation null (shuffle the parental tags within windows)
  gives an empirical distribution of segment shifts under no event; report an FDR per threshold.
- **Leave-trio-out evaluation.** Run more 1000 Genomes trios with a panel that excludes the evaluated trio; the
  shipped panel contains the two example trios. Known mosaic aneuploidies in 1kGP lymphoblastoid lines are catalogued
  in the literature and can serve as truth.

## Sex chromosomes

- **Y segmental events.** The Y is analysed only as a whole chromosome (copy number, mosaic loss); segmental Y events
  would need a male-only panel with Y bins aligned to one copy and a mappability-aware bin mask.
- **Panel Y rows.** `triokaryo panel` aligns each genome's Y to the two-copy level like the X; the shipped panel
  predates this and its Y rows, if present, are re-aligned at load time by a heuristic. Rebuild the shipped panel from
  the twelve VCFs when they are next available.
- **X inactivation and PAR.** The pseudoautosomal regions are excluded throughout; PAR1 could be analysed as autosomal
  (both parents contribute).

## Outputs and usability

- **Karyotype string per trio** in ISCN-like form (e.g. `47,XY,+12 pat`, `mos 46,XY,del(2)(q22q31)[0.67]`), and a
  standard CNV export (BED or VCF) for downstream tools.
- **Crossover map table** per meiotic event: position, parent, and the auxiliary track that resolved it (the
  centromere-anchored stage reading computes these; a dedicated table would make them queryable).
- **`triokaryo merge`**: a wrapper around `bcftools merge -0` for per-sample VCFs, with the PASS/biallelic SNV
  reduction, so that the example recipe is a single command.
- **Throughput for cohorts.** Parallelise the scan and the binning by chromosome; stream sites instead of holding every
  chromosome's arrays; a `--regions` option for a quick look at one chromosome.
- **Caller coverage.** Document which callers' VCFs carry `AD`, `DP` and `GQ` as required (GATK, DRAGEN, DeepVariant
  with `--vcf_stats`); accept `AD`-only VCFs by deriving `DP`.
