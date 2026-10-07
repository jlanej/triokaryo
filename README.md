# triokaryo

Large chromosomal events in a **trio**, read from its own **VCF**: for the child, the father and the mother, the
depth along the genome (LRR), the B-allele fraction at heterozygous sites (BAF) and the heterozygosity rate;
segmentation of each; every gain, loss and copy-neutral loss of heterozygosity with the **share of cells carrying
it** estimated twice (from the depth and, independently, from the B-allele bands); the **parent of origin** of
the child's events from the sites where the parents are opposite homozygotes; what is **inherited** and what is
**new**; figures fit for a manuscript; and the match against events called elsewhere (a depth tool such as
NGS-DOSE's karyotype, a clinical karyotype).

No cohort is needed and nothing is compared across families: every normalisation and every test is within the
trio. The VCF is the one a joint caller writes for the family (GATK genotype refinement, VQSR-filtered), with
`GT`, `AD`, `DP` and `GQ` per member.

> **No real data in this repository.** The tests and the demo run on a mock trio that `triokaryo mock` writes
> (GRCh38 lengths, simulated sites, planted events). Never commit a VCF, a pedigree, a CRAM or a sample ID.

## On real data: two 1000 Genomes trios

[`docs/example/`](docs/example/README.md) runs the method on public genomes (Illumina's DRAGEN re-analysis of the
1kGP high-coverage cohort) beside NGS-DOSE's alignment-free karyotype and DRAGEN's own CNV calls: NA12739's trisomy 12
(the extra copy paternal, over 17,678 informative sites), a 69% mosaic gain of 13q in the father that only the VCF and
DRAGEN see, and HG01103's loss of 2q (141 to 173 Mb in 67% of cells, the paternal copy) with a 14q loss beside it. The
pages are [NA12739](docs/example/NA12739/index.html), [HG01103](docs/example/HG01103/index.html) and the
[cohort](docs/example/cohort/index.html) (self-contained pages: open them from a clone, GitHub shows their
source); the recipe reproduces them.

## The phased bands

The B-allele fraction's noise is binomial and cannot be smoothed away site by site, but the trio gives it a sign.
Wherever the parents are opposite homozygotes the child's two alleles have known parents, and the fraction of the
mother's allele - the maternal fraction - sits at one half along a normal chromosome, at 1/2 + d or 1/2 − d along
an event by the parent of origin, at 1 where no paternal copy is left. Along a parent, the fraction of the allele
passed to the child (read where the child is homozygous) says whether an event lies on the homologue the child got.
Pooled over windows and fitted by total-variation denoising, this signed track is the clean reading the folded bands
cannot give; the LRR's copies split by it are the maternal and paternal copies along the child's genome (a trisomy's
extra copy, a deletion's missing one, a disomy's two from one parent, each with its parent named). The sites where
only one parent is homozygous phase the child too, but read the wrong way wherever the child carries two different
homologues of one parent - so they are drawn as auxiliary tracks that part from the main one along a meiotic
trisomy or heterodisomy and return at each crossover, which tells a meiotic error from a mitotic one. A scan of the
main track finds what the depth cannot: a gain or loss in a few per cent of cells, and a uniparental heterodisomy.
The raw B-allele fraction stays in every figure beside it.

## Install

```bash
pip install git+https://github.com/jlanej/triokaryo         # or: make venv && make test
# Apptainer / Singularity (HPC): the image CI publishes to GHCR
apptainer build triokaryo.sif docker://ghcr.io/jlanej/triokaryo:latest
```

## Quickstart

```bash
# a GC track, once per reference and bin size (1 Mb)
triokaryo gc-track --fasta GRCh38.fa --out gc.grch38.1mb.tsv

# one trio: the VCF holds the three members; the trios file names them (#kid dad mom kid_sex dad_sex mom_sex)
triokaryo run --vcf family.vcf.gz --pedigree trios.tsv --child KID --gc-track gc.grch38.1mb.tsv \
              --events ngsdose/karyotype/events.tsv --out out/KID
# or name the members outright
triokaryo run --vcf family.vcf.gz --child KID --father DAD --mother MOM --sex M,M,F --out out/KID

# many trios, gathered: every event, every X reading, and the concordance with the events given from elsewhere
triokaryo cohort --runs 'out/*' --events ngsdose/karyotype/events.tsv --out cohort

# the demo on the mock
make demo && open mock_out/KID/index.html
```

A trio of a 5-million-site WGS VCF takes about a minute and under 1 GB; `--thin 3` reads every third site for a
quick look (the bins keep hundreds of sites each).

**On real data use a panel.** A genome read against its own median shows the reference's depth structure - centromere
flanks, segmental duplications, the acrocentric short arms - and the bands parted by paralogous sequence as dozens of
events. `--panel` takes out what other genomes called the same way share (the median LRR, band deviation and
heterozygosity rate per bin), and leaves out the bins it cannot pin. Build one from any genomes counted the same way
(`triokaryo panel --vcfs ...`, five or more), or from a first pass's `bins.tsv` files when there are many trios
(`--runs 'out/*'`: the cohort is its own panel); `--panel 1kg-dragen` is a shipped one of twelve public 1000 Genomes
genomes for data called by DRAGEN. Per-sample VCFs (one genome each) are merged into a trio VCF with
`bcftools merge -0`: a sample without a record at a site is homozygous reference, which the trio reading takes as a
confident parental call.

## What it finds, and how

| signal | what it reads | what it shows |
| --- | --- | --- |
| **LRR** | log2 of the bin's median depth over the member's autosomal median; GC-corrected (a running median against the GC track) | a gain of one copy in a share *f* of the cells reads log2(1 + f/2); a loss log2(1 − f/2) |
| **BAF** | the alt-allele fraction at the member's heterozygous sites | a gain parts the bands to 1/(2+f) and (1+f)/(2+f) (1/3 and 2/3 at f = 1); a loss to (1−f)/(2−f) and 1/(2−f); a copy-neutral LOH to (1−f)/2 and (1+f)/2 |
| **het rate** | heterozygous calls over confident calls, per bin | zero under a loss of heterozygosity in every cell (a uniparental disomy, a deletion); unchanged under a trisomy |
| **child over parents' mean** | log2 of the child's depth over the parents' mean, site by site, per bin | the within-family difference: zero where the child inherited what the parents carry |

Each member's LRR is segmented chromosome by chromosome (binary segmentation on a robust noise scale, segments
merged when alike); a segment is a gain or loss when its mean is beyond 0.07 (a share of cells of about 10%) and
three standard errors. The band deviation *d* at the segment's heterozygous sites is a maximum-likelihood fit of
alt ∼ Binomial(depth, ½ ± d), which does not inflate at low depth as |BAF − ½| does; it gives *f* a second time
(`f_baf` beside `f_lrr`). Copy-neutral LOH is searched where the depth called nothing: the band deviation above the
member's own, or no heterozygous calls at all. The X is read against the member's own X copy state (one or two
copies, checked against the pedigree's sex: a 47,XXY reads "X copies 2 in a reported male"), so its events are
mosaic changes of what the member has; the Y has too few sites in a VCF and is left to a depth tool. A segment
covering 90% of a chromosome (arm) is `whole` (`p`, `q`), else a `stretch`.

**Parent of origin.** At a site where the father is 0/0 and the mother 1/1 (or the reverse) the child's two
alleles have known parents, and the child's alt fraction says whose copy is extra (a gain: the duplicated
parent's allele at (1+f)/(2+f)), lost (a loss: the retained parent's at 1/(2−f)) or doubled (an LOH: the retained
parent's at (1+f)/2). The log-likelihood ratio of the two assignments over the event's informative sites names the
parent and says how sure (`origin`, `origin_llr`, `origin_n`). **Inheritance:** the same event in a parent
(reciprocal overlap ≥ 50%, the same type), with "in a share of the parent's cells" where the parent is mosaic;
a parent's event is "passed to the child" or not. **Mendelian errors** within each event are reported beside the
genome's baseline: a deletion and a uniparental disomy break Mendel at the informative sites, a trisomy and a
mosaic do not.

## Output (per trio, under `--out`)

- `events.tsv` — one row per event: sample, role, chrom, start, end, span, type (gain, loss, LOH, UPD), `source`
  (depth, bands, phased), `f` (the depth's where there is one, else the bands'), `f_lrr`, `f_baf`, `f_phase`,
  `lrr`, `d_hat`, `llr_baf`, `phase_shift`, `n_phased`, `homologues` (one or two, for a child's gain or disomy),
  `start_fine`, `end_fine` (the edges at site resolution), `het_rate`, `het_rate_rel`, `mie_rate`, `origin`,
  `origin_phase`, `inheritance`, `external` (the given events it overlaps), `note` (an LOH in every cell with no
  Mendelian errors is a run of homozygosity, with them a uniparental isodisomy; pieces joined across a masked gap;
  the phased share against the depth's).
- `phased.tsv` — the pooled windows of each member's phased track (fraction, error, step fit, copies, the auxiliary
  tracks); `phased_rejected.tsv` the scan's segments set aside, with the reason.
- `bins.tsv` — every bin: GC, the panel's median, spread and masks, and per member the depth, LRR, GC- and
  panel-corrected LRR, calls, heterozygous calls, het rate (and relative to the panel), band deviation; the two
  within-trio tracks.
- `summary.tsv`, `summary.json` — the trio, the sites used, the X copies per member with the sex check, the
  genome's Mendelian-error rate, the event counts, the parameters.
- `external.tsv` — the given events (`--events`) and whether a triokaryo event matched each.
- `guide.html` — how to read every row, colour, call and column, with pattern cards of what each kind of event
  looks like; beside every page, and `triokaryo guide --out`.
- `figures/` — `genome` (LRR with its step fit and BAF with the pooled phased fraction per member, the child over
  the parents' mean, the child's maternal and paternal copies, the calls, the given events) and `chrom_<chrom>` for
  every chromosome with an event (per member: LRR with its step fit and calls, the raw BAF with the informative sites
  coloured by parent, the phased fraction with its windows, step fit and auxiliary tracks, the two homologues' copies,
  the het rate), each with a `.txt` sidecar (title, caption, key) and `legends/<name>_legend.*`.
- `cohort/` (`triokaryo cohort --runs 'out/*' --events ...`) — the cohort report: the counts, the landscape figure
  (events per chromosome; one row per trio), every event sortable and filterable with links to its trio's page and
  figure, the trios with their quality readings, the concordance with the given events, the segments the phased scan
  set aside by region, and the guide; `events.all.tsv`, `summary.all.tsv`, `concordance.tsv`, `flags.tsv`,
  `rejected.all.tsv`. `triokaryo report --runs 'out/*'` rebuilds a run's page from its tables after a change.
- `index.html` — self-contained: the events, the given events, every figure with its caption and key.

`triokaryo cohort` writes `events.all.tsv`, `summary.all.tsv`, `concordance.tsv` (each given event: matched or
not) and an `index.html` linking the trios.

## Events given from elsewhere (`--events`)

A TSV with `sample chrom start end label` (and optionally `type`); or NGS-DOSE's `karyotype/events.tsv`
(`sample chrom span start_mb end_mb label kind`), recognised by its columns. They are drawn as brackets above the
LRR and matched to triokaryo's events by sample and reciprocal overlap ≥ 50%. A depth tool cannot see a
copy-neutral LOH, so those are never "unmatched" against it.

## The mock and the tests

`triokaryo mock --out dir` writes a trio VCF (GRCh38 lengths, 60 sites per Mb, Poisson depth with a GC bias per
sample and a synthetic GC track, genotypes called as a caller would) with ten planted events: a maternal
trisomy 21; a mosaic +12 (30%, paternal homologue); +10q in the father, inherited; a de novo 18q deletion with the
paternal copy lost; a maternal isodisomy 7; a mosaic copy-neutral LOH of 6p (40%); a mosaic +8 in the mother (15%),
not passed; a 6-Mb deletion in the father, inherited. `truth.json` lists them; `events.external.tsv` lists the
copy-number ones as another caller would. `pytest` checks every event is found with its type, share of cells
(both estimates), parent of origin and inheritance, that nothing else is called, that a null trio calls nothing,
that the GC correction removes the GC bias, and that a 47,XXY child reads two X copies. `--xxy` and `--no-events`
make those mocks.

## Limits

- Sites only where the caller wrote a PASS biallelic SNV: a region without calls (centromeres, large gaps) has no
  bins. Bins are 1 Mb; events under about 5 Mb are not sought (`--min-len`, `--bin`). Without a panel the reference's
  own structure reads as events; a panel of fewer than five genomes, or built from the trio alone, follows the trio's
  own events (an inherited event sits in two of three genomes and vanishes).
- A gain in a share of cells under about 10% is below the reporting floor (`--min-f`); the band deviation needs
  depth to see it.
- The Y is not read. A 46,XY/47,XXY mosaic in a male reads as an X gain in a share of the cells; the parent of
  origin is read on autosomes only.
- A deletion in all cells has no heterozygous sites, so `f_baf` is undefined there (the depth gives `f`); its
  Mendelian-error rate is the tell.
- The parent-of-origin models assume one event per region; a trisomy from a meiosis II or mitotic error (one
  homologue doubled) still reads its parent, since the opposite-homozygote sites see the parent, not the homologue.

## Development

```bash
make venv && make test        # pytest on the mock (about 30 s)
make demo                     # the mock end to end into mock_out/
make docker-test              # the image, run on the mock
```

CI runs the tests on two Pythons, builds the image, runs the mock inside it, and pushes to GHCR on `main` and
tags. See `docs/METHODS.md` and `docs/OUTPUT.md`.

## License

MIT.
