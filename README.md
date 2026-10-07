# triokaryo

Detection and characterisation of large chromosomal events in a parent–offspring trio from the trio's small-variant
VCF: whole-chromosome and segmental gains and losses, copy-neutral loss of heterozygosity (CN-LOH) and uniparental
disomy (UPD), constitutional or mosaic.

For each member (child, father, mother) triokaryo computes three tracks from the PASS biallelic SNVs: a binned log R
ratio (LRR) from read depth, the B-allele frequency (BAF) at heterozygous sites, and the heterozygosity rate. It
segments each track and, using the pedigree, phases the child's alleles by transmission to obtain a signed
allelic-imbalance track. Each event is reported with its type; its mosaic cell fraction *f* estimated independently
from depth, from the folded BAF bands and from the phased track; its parent of origin; whether the two copies from
the named parent are one homologue or two (mitotic versus meiotic origin); its boundaries at site resolution; its
Mendelian-error rate; and whether it is inherited or de novo. Figures are written as PNG, SVG and PDF with separate
legends, and calls from other methods (a depth-based karyotype such as NGS-DOSE, a CNV caller, a clinical
karyotype) can be supplied for comparison.

All normalisation and all tests are within the trio; no cohort is required. The input is any VCF holding the three
members with per-sample `GT`, `AD`, `DP` and `GQ`: a joint-called family VCF, or per-sample VCFs merged with
`bcftools merge -0`.

> **No real data in this repository.** The tests and the demo use a simulated trio written by `triokaryo mock`
> (GRCh38 chromosome lengths, simulated sites, planted events). Do not commit VCFs, pedigrees, CRAMs or sample IDs.

## Example: two 1000 Genomes trios

[`docs/example/`](docs/example/README.md) applies the method to public genomes (Illumina's DRAGEN 3.7.6 re-analysis
of the 1000 Genomes high-coverage cohort) and compares the calls with NGS-DOSE's alignment-free karyotype and with
DRAGEN's CNV calls. In trio NA12739: a trisomy 12 in the child (paternal extra copy, one homologue, 17,630 phased
sites) and a mosaic gain of 13q in the father (69% of cells) that NGS-DOSE did not report. In trio HG01103: a mosaic
loss of 2q (141–173 Mb, 67% of cells, paternal copy) with a subclonal extension to 180 Mb, and a mosaic loss of 14q.
Pages: [NA12739](docs/example/NA12739/index.html), [HG01103](docs/example/HG01103/index.html), and the
[cohort report](docs/example/cohort/index.html). They are self-contained HTML files; open them from a clone (GitHub
displays their source). The fetch script in `docs/example/` reproduces them.

## Transmission phasing

The BAF at a heterozygous site is 1/2 ± d with the sign of d unknown, so BAF is conventionally read folded as
|BAF − 1/2|, an estimator biased upward by binomial noise (expectation 0.07 at 30× when d = 0). The trio supplies the
sign. Where the parents are opposite homozygotes (father 0/0 and mother 1/1, or the reverse), the parental origin of
each of the child's alleles is known regardless of the child's genotype call, and the fraction of the child's reads
carrying the maternal allele (the *maternal fraction*) can be computed directly. It is 1/2 where both homologues are
present in equal copy, 1/2 + d or 1/2 − d along an event according to the parent of origin, and 1 where no paternal
copy is present (paternal deletion, maternal isodisomy or heterodisomy). In a parent, at heterozygous sites where the
child is homozygous, the transmitted allele is known, and the fraction of reads carrying it shows whether an event in
that parent lies on the transmitted homologue.

These main tracks are pooled in windows of a fixed number of sites and fitted by exact one-dimensional total-variation
denoising. The LRR-derived copy number split by the fitted fraction gives per-homologue copy number along the child
(maternal, paternal) and along each parent (transmitted, untransmitted). A phasing error flips one site's
contribution and so attenuates a signal rather than creating one.

Sites where only one parent is homozygous also phase the child, but only under the assumption that each parent
contributed one homologue; where the child carries two different homologues of one parent (meiotic trisomy,
heterodisomy), sites tagged by the other parent's homozygous genotype are misassigned. These are therefore drawn as
auxiliary tracks: they depart from the main track where the two homologues differ and rejoin it after each crossover,
which distinguishes meiotic from mitotic (or meiosis II) origin and maps the crossovers. A segmentation scan of the
main track detects events below the depth detection limit (gains or losses in a few per cent of cells) and
uniparental heterodisomy, which depth, folded BAF and heterozygosity rate all miss. The raw BAF is retained in every
figure.

## Install

```bash
pip install git+https://github.com/jlanej/triokaryo         # or: make venv && make test
# Apptainer / Singularity (HPC): the image CI publishes to GHCR
apptainer build triokaryo.sif docker://ghcr.io/jlanej/triokaryo:latest
```

Python ≥ 3.9 with pysam, numpy and matplotlib.

## Quickstart

```bash
# a GC track, once per reference and bin size (1 Mb)
triokaryo gc-track --fasta GRCh38.fa --out gc.grch38.1mb.tsv

# one trio: the VCF holds the three members; the trios file names them (#kid dad mom kid_sex dad_sex mom_sex)
triokaryo run --vcf family.vcf.gz --pedigree trios.tsv --child KID --gc-track gc.grch38.1mb.tsv \
              --panel 1kg-dragen --events ngsdose/karyotype/events.tsv --out out/KID
# or name the members directly
triokaryo run --vcf family.vcf.gz --child KID --father DAD --mother MOM --sex M,M,F --out out/KID

# many trios: all events, per-trio quality metrics, concordance with the supplied events
triokaryo cohort --runs 'out/*' --events ngsdose/karyotype/events.tsv --out cohort

# the demo on the simulated trio
make demo && open mock_out/KID/index.html
```

A trio VCF of about 5.5 million sites runs in 60–75 s in under 1 GB of memory. `--thin n` reads every n-th site for
a quick pass (1-Mb bins still hold hundreds of sites).

**Use a reference panel on real data.** Normalised against its own median alone, a genome shows the reference's
depth structure (centromere flanks, segmental duplications, acrocentric short arms) and BAF bands split by
paralogous sequence, which segment into spurious events. `--panel` subtracts the per-bin median LRR, band deviation
and heterozygosity rate of other genomes processed the same way and masks bins the panel cannot characterise.
Build a panel from any genomes called the same way (`triokaryo panel --vcfs ...`; five or more), or, with many trios,
from a first pass's `bins.tsv` files (`--runs 'out/*'`). `--panel 1kg-dragen` is a shipped panel of twelve public
1000 Genomes genomes for data called by DRAGEN 3.7.6. Per-sample VCFs are merged into a trio VCF with
`bcftools merge -0`; a sample with no record at a site is then written as homozygous reference, which the trio
analysis accepts as a confident parental genotype.

## Signals

| signal | definition | expected value under an event in a cell fraction *f* |
| --- | --- | --- |
| **LRR** | log2 of the bin's median depth over the member's autosomal median, GC-corrected against a GC track by a running median | one-copy gain: log2(1 + f/2); one-copy loss: log2(1 − f/2) |
| **BAF** | alt-allele read fraction at the member's heterozygous sites | gain: bands at 1/(2+f) and (1+f)/(2+f) (1/3, 2/3 at f = 1); loss: (1−f)/(2−f) and 1/(2−f); CN-LOH: (1−f)/2 and (1+f)/2 |
| **heterozygosity rate** | heterozygous calls over confident calls per bin | zero under a constitutional loss of heterozygosity (isodisomy, deletion, run of homozygosity); unchanged under trisomy |
| **phased fraction** | fraction of reads carrying the maternal allele (child) or the transmitted allele (parent), at phased sites | 1/2 ± d with the sign giving the parent of origin; 1 where only one parent's copies are present |
| **child vs parental mean** | log2 of the child's depth over the parents' mean depth, per site, median per bin | zero where the child's copy number equals the parents' mean; drawn, not called |

## Calling

Each member's LRR is segmented per chromosome by binary segmentation on a robust noise scale (median absolute
deviation of first differences), with boundary refinement and merging of segments whose means differ by less than
three standard errors. A segment is a gain or loss when |mean LRR| ≥ max(0.07, 3 SE), that is a cell fraction of
about 10% or more, and spans at least five bins (5 Mb at the default bin size). The band deviation d of a segment is
the maximum-likelihood estimate under alt ~ Binomial(depth, 1/2 ± d) over its heterozygous sites, which, unlike the
mean of |BAF − 1/2|, is unbiased at low depth; it gives a second estimate of *f* (`f_baf` beside `f_lrr`). CN-LOH is
sought where depth called nothing: a band deviation above the member's own baseline, or a heterozygosity rate below
0.35 of the member's own (constitutional). The phased scan then segments the pooled phased track outside the called
events and reports shifts of ≥ 0.015 (gain or loss in about 6% of cells) over ≥ 2 Mb, typed by the direction of the
depth over the same bins, or as CN-LOH or heterodisomy when the depth is flat. The X is analysed relative to the
member's own X copy number (1 or 2, compared with the pedigree sex: 47,XXY is reported as "X copies 2 in a reported
male"), so X events are mosaic changes of that state; the Y is not analysed. A segment covering ≥ 90% of a
chromosome (arm) is `whole` (`p`, `q`), otherwise a `stretch`.

**Parent of origin** is estimated twice. (i) At informative sites (parents opposite homozygotes), the child's alt
read count is modelled as binomial with a success probability determined by the event type, *f*, and which parent
contributed the alt allele; the log-likelihood ratio of the two parental assignments over the event's sites gives
`origin`, `origin_llr` and `origin_n`. (ii) The sign of the phased shift gives `origin_phase`. A disagreement is
flagged. **Inheritance:** a child's event is inherited when a parent carries an event of the same type on the same
chromosome with reciprocal overlap of at least 50% (the intersection covers at least half of each segment), "in a
share of the parent's cells" when the parent's *f* < 0.8; otherwise it is de novo. A parent's event is marked
transmitted or not by the same rule.
**Mendelian errors** are counted per event and compared with the genome-wide rate: a constitutional deletion or an
isodisomy produces errors at informative sites, a trisomy or a mosaic event does not, and a run of homozygosity has
none.

## Output (per trio, under `--out`)

- `events.tsv`: one row per event with sample, role, coordinates, span, type (gain, loss, LOH, UPD), source (depth,
  bands, phased), `f`, `f_lrr`, `f_baf`, `f_phase`, LRR, band deviation and its likelihood ratio, phased shift and
  site counts, homologue classification, site-resolution boundaries, heterozygosity rate, Mendelian-error rate,
  parent of origin (both estimates), inheritance, overlapping supplied events, and notes.
- `phased.tsv`: the pooled windows of each member's phased track (fraction, error, step fit, per-homologue copies,
  auxiliary tracks); `phased_rejected.tsv`: segments the phased scan rejected, with the reason.
- `bins.tsv`: per bin, GC, the panel's values and masks, and per member depth, LRR, corrected LRR, call counts,
  heterozygosity rate and band deviation; the within-trio depth tracks.
- `summary.tsv`, `summary.json`: the trio, the sites used, the X copy number per member with the sex check, the
  genome-wide Mendelian-error rate, event counts, parameters.
- `external.tsv`: the supplied events (`--events`) and whether each was matched.
- `guide.html`: how to read every figure row, colour, call and column, with pattern cards of each event type (also
  `triokaryo guide --out`).
- `figures/`: `genome` (per member LRR with step fit and BAF with the pooled phased fraction; child vs parental mean;
  the child's maternal and paternal copy number; calls; supplied events) and `chrom_<chrom>` for every chromosome
  with an event (per member: LRR with calls, BAF with informative sites coloured by parent, phased fraction with
  windows, step fit and auxiliary tracks, per-homologue copy number, heterozygosity rate). Each has a `.txt` sidecar
  (title, caption, key) and a legend image under `legends/`.
- `index.html`: a self-contained page with the events, the supplied events, and every figure with its caption.

`triokaryo cohort --runs 'out/*' --events ...` gathers many runs into a cohort report: counts, a landscape figure
(events per chromosome; one row per trio), a sortable table of every event linked to its trio's page and figure,
per-trio quality metrics, concordance with the supplied events, and the rejected phased segments by region; with
`events.all.tsv`, `summary.all.tsv`, `concordance.tsv`, `flags.tsv` and `rejected.all.tsv`. `triokaryo report
--runs 'out/*'` rebuilds a run's page, sidecars, legends and guide from its tables without the VCF.

## Supplied events (`--events`)

A TSV with columns `sample chrom start end label` (optionally `type`), or NGS-DOSE's `karyotype/events.tsv`
(`sample chrom span start_mb end_mb label kind`), recognised by its columns. Supplied events are drawn as brackets
above the LRR and matched to triokaryo's events of the same sample when the intersection covers at least half of
the shorter segment. Copy-neutral events are never counted as unmatched against a depth-based method.

## Simulated trio and tests

`triokaryo mock --out dir` writes a trio VCF (GRCh38 lengths, 60 sites per Mb, Poisson depth with a per-sample GC
bias and a matching synthetic GC track, genotypes called from the simulated read counts) with eleven planted events
(nine lesions, two of them inherited): a maternal meiotic trisomy 21; a 30% mosaic trisomy 12 (paternal homologue
duplicated); a constitutional gain of 10q in the father, inherited by the child; a de novo deletion of 18q (paternal
copy); a maternal isodisomy 7; a 40% mosaic CN-LOH of 6p; a 15% mosaic trisomy 8 in the mother, not transmitted; a
6-Mb deletion in the father, inherited; and a maternal heterodisomy 15. `truth.json` lists them and
`events.external.tsv` lists the copy-number ones as an external caller would. `pytest` verifies that every planted
event is detected with its type, cell fraction (all three estimates), parent of origin, homologue classification and
inheritance; that nothing else is called; that a null trio yields no calls; that the GC correction removes the GC
bias; that a 47,XXY child reads two X copies; that the step fit equals the exact total-variation solution; and that
8% events planted in a dense simulation are recovered by the phased scan. `--xxy`, `--no-events` and `--low-share`
produce those variants.

## Limitations

- Only PASS biallelic SNVs are used; regions without calls (centromeres, large gaps) have no bins. Depth-based calls
  need at least five 1-Mb bins (`--min-len`, `--bin`) and a cell fraction of about 10% (`--min-f`); the phased scan
  extends this to 2 Mb and about 6% where the phased site density allows.
- Without a reference panel the reference's own depth structure is called as events. A panel of fewer than five
  genomes masks every bin; a panel built from the trio alone follows the trio's own events.
- The Y chromosome is not analysed. A 46,XY/47,XXY mosaic is reported as a mosaic X gain. Parent of origin is estimated on
  autosomes only.
- A constitutional deletion has no heterozygous sites, so `f_baf` is undefined there (`f` comes from depth); its
  Mendelian-error rate and phased fraction of 1 identify it.
- The parent-of-origin likelihood assumes one event per region and is conditional on an event being present; it does
  not test for the event. A trisomy from a meiosis II or mitotic error (one homologue duplicated) still yields its
  parent, since informative sites resolve the parent, not the homologue.
- Supplied events are matched by intersection over the shorter segment, deliberately lenient so that a caller's
  fragmented segments match one event; inheritance requires reciprocal overlap of at least 50%.

## Development

```bash
make venv && make test        # pytest on the simulated trio (a few minutes; figures are rendered)
make demo                     # the simulated trio end to end into mock_out/
make docker-test              # the image, run on the simulated trio
```

CI runs the tests on two Python versions, builds the image, runs the simulated trio inside it, and pushes to GHCR on
`main` and tags. See `docs/METHODS.md`, `docs/OUTPUT.md` and, for planned work, `docs/ROADMAP.md`.

## License

MIT.
