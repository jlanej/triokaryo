# Output

## events.tsv (per trio) and events.all.tsv (cohort)

| column | meaning |
| --- | --- |
| sample, role | the member (child, father, mother); the cohort table adds `trio` |
| chrom, start, end | the event in bp: the first bin's start (0-based) and the last bin's end; for a phased-scan event, its first and last window sites |
| start_fine, end_fine, edge_sites | boundaries at site resolution from the phased sites (NA at a chromosome end or with too few sites), and the smaller of the two site counts used |
| span | whole, p, q or stretch |
| type | gain, loss, LOH (copy-neutral loss of heterozygosity), UPD (uniparental heterodisomy) |
| source | depth (LRR segmentation), bands (band deviation or heterozygosity rate), phased (the phased scan) |
| f | cell fraction: f_lrr for a depth-called gain or loss, otherwise f_baf |
| f_lrr, f_baf, f_phase | cell fraction from depth, from the folded band deviation and from the phased shift |
| lrr, lrr_se, n_bins | the segment's mean LRR, its standard error from the chromosome's noise scale, and its bin count |
| d_hat, llr_baf | band deviation at the segment's heterozygous sites (maximum likelihood) and its log-likelihood ratio against d = 0; for a phased-scan event d_hat is the absolute phased shift and llr_baf is NA |
| phase_shift, phase_se, n_phased | the phased fraction's shift from 1/2 over the event (maternal allele in the child, transmitted allele in a parent), its binomial error, and the site count |
| origin, origin_llr, origin_n | parent of origin from the informative sites (child, autosomes), the log-likelihood ratio and the site count |
| origin_phase | from the sign of the phased shift: the parent of the extra, lost or retained copy (child), or whether the event lies on the transmitted homologue (parent) |
| homologues, hetero_share | for a child's gain, LOH or UPD with a parent named: whether that parent's two copies are one homologue or two, and the share of the event's windows in which they differ |
| stage, centromere, n_crossovers, crossovers | for a child's whole-chromosome gain or heterodisomy with a parent named: the meiotic stage (meiosis I, meiosis II, or mitotic / meiosis II without a crossover) from the state of the two copies nearest the centromere, that centromeric state, and the crossovers as changes of state along the chromosome (positions in Mb) |
| het_rate, het_rate_rel, n_het, n_called | heterozygosity rate in the segment, its ratio to the member's own (or to the panel's), and the counts |
| mie_rate | Mendelian-error rate among the segment's confident sites (the genome-wide rate is in the summary) |
| inheritance | child: inherited from the father / the mother, or new; parent: passed to the child or not |
| external | labels of the supplied events (`--events`) overlapping the segment |
| note | annotations: run of homozygosity; isodisomy; no heterozygous calls; pieces joined; phased-versus-depth disagreement; phased-scan typing; heterodisomy; parent-of-origin disagreement; loss on a single X; whole-chromosome X or Y against the pedigree sex (47,XXY, 45,X, mosaics, loss of Y) with the source of the Y estimate |

## phased.tsv

One row per pooled window of each member's main phased track: role, chrom, start, end, mid (the first, last and
median site positions), n_sites, depth (summed), frac (the pooled fraction), se, shared (deviating in two or more
members), step (the step fit), lrr (the window's bin), copies_tagged and copies_other (per-homologue copy number),
and the auxiliary tracks' pooled fraction and site count per window (`aux_mother_hom`, `aux_father_hom` for the
child; `aux_child_het` for a parent).

## phased_rejected.tsv

Segments of the phased scan that were not reported, with the reason: span under 2 Mb; windows inconsistent with the
segment mean; or, with flat depth, a shift under the bands-only floor of 0.025.

## bins.tsv

Per bin: chrom, start, end, gc; the panel's median LRR and robust SD, the mask, the panel's band deviation and the
band mask; per member `n_sites`, `depth`, `lrr`, `lrr_gc` (GC- and panel-corrected), `n_called`, `n_het`, `het_rate`,
`bdev`, `het_rel` (relative to the panel); `child_vs_mid` and `father_vs_mother`.

## summary.tsv and summary.json

The trio, its members and their sexes; records read, sites used and records skipped by reason; the genome-wide
Mendelian-error rate; per member the event counts by type, the X and Y copy numbers (`x_copies`, `y_copies`; raw
values `x_copies_raw`, `y_copies_raw`), the sex-chromosome complement (`sex_karyotype`) with its check against the
pedigree sex (`sex_check`; `x_check` the X-only wording), the median autosomal depth, the phased-site count, the
window size, the reference bias, the number of phased-scan events and of shared windows; the father/son Y depth
ratio (`y_father_son_log2`, `y_father_son_sites`) and whether the panel carried Y rows (`y_panel`); the parameters.
`summary.json` also holds every event, every supplied event and the per-member sex-chromosome state.

## figures/

`genome.{png,svg,pdf}`: eight rows; per member the LRR with its step fit and calls and the BAF with the pooled
phased fraction, then the child's depth over the parents' mean, and the child's maternal and paternal copy number.
`chrom_<chrom>.{png,svg,pdf}`: five rows per member; LRR with step fit and calls, BAF with the child's informative
sites coloured by parent, the phased fraction (sites, windows, step fit, auxiliary tracks), per-homologue copy
number, and heterozygosity rate. Each figure has a sidecar `<name>.txt` (title, caption, key with hex colours) and a
legend image `legends/<name>_legend.{png,svg,pdf}`; the image itself carries no legend or title.

Colours (Okabe–Ito): LRR points dark grey; BAF light grey; pooled phased fraction bluish green `#009E73` (hollow where
shared); step fits black; gain vermillion `#D55E00`; loss blue `#0072B2`; LOH and UPD reddish purple `#CC79A7`. In the
child's panels, maternal alleles, copies and the auxiliary track that departs under two maternal homologues are
`#CC79A7` and the paternal equivalents `#E69F00`; a parent's transmitted-homologue copies are `#D55E00` and the
untransmitted `#56B4E9`; supplied events are black brackets; reference lines dotted grey.

## guide.html

Written beside every trio page and the cohort report: the meaning of every figure row, colour, event type (pattern
cards drawn from idealised tracks: gain, loss, CN-LOH, isodisomy, heterodisomy, run of homozygosity, meiotic and
mitotic trisomy), calling rule and column. Also `triokaryo guide --out guide.html [--figures dir]`.

## cohort/ (`triokaryo cohort`)

`index.html`: counts (trios, events by type, member and source, de novo and inherited, phased-scan events,
depth calls doubted by the phased track, runs of homozygosity, X copy numbers disagreeing with the pedigree,
concordance), `figures/landscape.*` (events per chromosome by type; one row per trio with each event a bar coloured
by type, the child's thick, a parent's thin above or below, a black line over a bar for a matched event), a sortable
and filterable table of every event linked to its trio's page and chromosome figure, per-trio quality metrics,
the concordance with the supplied events, and the rejected phased segments by region. Tables: `events.all.tsv`,
`summary.all.tsv` (per trio: members, sexes, sites, Mendelian-error rate, X copy numbers and checks, depths, phased
sites, shared windows, rejected segments, events, flagged events, run time, run directory), `concordance.tsv`,
`flags.tsv` (every event with a note), `rejected.all.tsv`; plus `guide.html` and `figures/patterns_*.*`.

## triokaryo report --runs

Rebuilds a run's page, figure sidecars, legend images and guide from its tables and existing figures, without the
VCF: after a change to the page or the key, or to add the guide to an older run.
