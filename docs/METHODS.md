# Methods

Notation: *f* is the fraction of cells carrying an event (cell fraction); *d* is the deviation of an allelic fraction
from 1/2. Bins are 1 Mb unless `--bin` is set. Defaults are given in parentheses; those with a flag are command-line
options.

## Input and site selection

The VCF is read once with pysam, restricted to the three members. A record is used if it is a biallelic SNV whose
FILTER is PASS or empty. Per member the following are stored: depth (`DP`, else the sum of `AD`), alt-allele depth
(`AD[1]`), genotype class (0/0, 0/1, 1/1 or missing; a haploid call is the homozygote of its allele) and `GQ`.
Records on contigs absent from the assembly table (alternate loci, decoys, chrM) are skipped. `--thin n` keeps every
n-th usable record. A homozygous-reference genotype with no `AD`, `DP` or `GQ`, as written by `bcftools merge -0` for
a sample without a record, is flagged and accepted as a confident genotype wherever a genotype rather than an
allelic fraction is required (parental genotypes at informative sites, Mendelian-error counting, phasing tags).

A call is *confident* when DP ≥ `--min-dp` (8) and GQ ≥ `--min-gq` (20).

## Bins and per-bin statistics

Fixed-width bins from each chromosome start. Per bin and member: the number of sites with depth; the **depth**, the
median DP over them; the numbers of confident calls and of confident heterozygous calls; the **heterozygosity rate**,
heterozygous over confident calls; and the **band deviation**, the median |BAF − 1/2| over the heterozygous sites
when there are at least five. **LRR** = log2(depth / m), where m is the member's median depth over autosomal bins
with ≥ 20 sites.

**Within-trio tracks.** For sites with depth in all three members, log2 of the child's depth over the parents' mean
and log2 of the father's over the mother's, summarised per bin as the median over its sites (≥ 5 sites). They cancel
site-level depth effects shared by the three libraries and are plotted, not segmented.

**GC correction** (`--gc-track`; `triokaryo gc-track` writes each bin's GC fraction from the reference). Autosomal
bins with ≥ 20 sites are sorted by GC and a running median over 41 bins gives the expected LRR at each GC value; it
is subtracted from every bin by interpolation and the result is re-centred so the autosomal median is zero. The
median is insensitive to the events themselves. The correction is skipped when fewer than 123 bins are available;
without a track the LRR is used as it is.

**Reference panel** (`--panel`). `triokaryo panel` processes other genomes identically (per-sample or multi-sample
VCFs, or earlier runs' `bins.tsv`) and records, per bin, the number of genomes, the median self-normalised LRR (each
genome's X shifted to two copies by its own X median) and its robust SD (1.4826 × MAD), the median band deviation and
the median heterozygosity rate. With a panel, each member's LRR is its self-normalised LRR minus the panel median,
re-centred on the autosomes; bins with fewer than five genomes, no median, or a robust SD above 0.25 are masked and
excluded from calling. The band deviation is read net of the panel's regional excess over its genome-wide median;
bins whose panel band deviation exceeds that median by more than 0.03 (paralogous sequence) are excluded from the
CN-LOH search and from the phased tracks; the heterozygosity rate is expressed relative to the panel's. GC
correction, when requested, is applied to the panel-corrected LRR.

## Segmentation of the depth track

Per member and chromosome, over unmasked bins with ≥ 20 sites (`min_sites`). The noise scale σ is 1.4826 × MAD of the
first differences, divided by √2. Binary segmentation splits recursively at the position maximising
|mean(left) − mean(right)| / (σ √(1/n_L + 1/n_R)) while the statistic is ≥ `--z` (5) and both sides hold ≥ `--min-len`
(5) bins, to at most 20 segments per chromosome. Each boundary is then moved within ±4 bins to the position minimising
the two adjacent segments' within-segment sums of squares, and adjacent segments whose means differ by less than
three standard errors are merged.

**Copy-number calls.** A segment of ≥ `min_len` bins with |mean LRR| ≥ max(`--min-abs` (0.07), 3 SE), SE = σ/√n, is a
**gain** (mean > 0) or a **loss**, with f_lrr = 2(2^LRR − 1) or 2(1 − 2^LRR); 0.07 corresponds to f ≈ 0.10. Events
with f < `--min-f` (0.10) are not reported. Adjacent events of one member and type that are separated by at most
three bins and differ in LRR by at most 0.08 (for LOH, in f_baf by at most 0.15) are joined, with the LRR re-averaged
over the pieces.

**Span.** A segment covering ≥ 90% of the chromosome's usable bins is `whole`; ≥ 90% of one arm's usable bins and
< 50% of the other arm's, `p` or `q`; otherwise `stretch`. Arm boundaries are the GRCh38 centromere starts.

## B-allele band deviation

For a segment, d is the maximum-likelihood value on a grid (0 to 0.475 in steps of 0.005) of
Σ log[½ Binom(k; n, ½ + d) + ½ Binom(k; n, ½ − d)] over the member's confident heterozygous sites in the segment, with
k the alt depth and n the depth (≥ 5 sites); `llr_baf` is the log-likelihood gain over d = 0. Each site's contribution
is weighted by its depth through the binomial, so the estimate does not inflate at low depth as the mean of
|BAF − 1/2| does. f_baf = 4d/(1 − 2d) for a gain, 4d/(1 + 2d) for a loss and 2d for CN-LOH.

## Copy-neutral loss of heterozygosity

Sought in bins the depth did not call. Two tracks are segmented as above: the band deviation (bins with ≥ 5
heterozygous calls and not excluded by the panel; other bins set to the member's autosomal median) and the
heterozygosity rate (bins with ≥ `min_called` confident calls, where `min_called` is 50 or, in sparser data, half the
member's typical count per bin, at least 10). With a panel, the panel-relative versions of both tracks are used. A
candidate is a segment of ≥ `min_len` bins, less than half covered by depth calls and less than half excluded by the
panel, whose mean band deviation is ≥ 0.03 above the member's median or whose mean heterozygosity rate is ≤ 0.35 of
the member's median. It is called **LOH** when its site-level d ≥ 0.04 with `llr_baf` ≥ 10 (f = 2d, about 8% of cells
or more), or when its heterozygosity rate is ≤ 0.35 of the member's own over ≥ max(40, 2 × `min_called`) confident
calls (constitutional; f = 1; noted "no heterozygous calls"). With a panel, d is first reduced by the panel's regional
band-deviation excess.

## The X chromosome

Bins overlapping the pseudoautosomal regions are excluded. The member's X copy number is round(2 × 2^median LRR) over
the X (≥ 5 bins), reported beside the pedigree sex (`x_check`), and the X LRR is re-centred on that median, so X events
are mosaic changes of the member's constitutional X state. A single X has no heterozygous sites, so no CN-LOH is sought
there, and a loss on a single X is annotated. The Y is not analysed. The phased scan includes the X only for a member
recorded as female.

## Transmission phasing

**Site classes.** For the child, the *main* class is the sites where the parents are opposite homozygotes (both
confident) and the child has a confident call of any genotype; the tagged allele is the mother's. Two *auxiliary*
classes are the sites where exactly one parent is homozygous and the child and the other parent are heterozygous:
`mother_hom` (mother homozygous, father heterozygous; tag the mother's allele) and `father_hom` (father homozygous,
mother heterozygous; tag the allele the father did not transmit). For a parent, the main class is the parent's
confident heterozygous sites where the child is confidently homozygous, excluding sites where the other parent is
homozygous for the opposite allele; the tagged allele is the child's, that is the transmitted one. The auxiliary class
`child_het` is the sites where the parent and the child are heterozygous and the other parent homozygous; the tag is
the allele the other parent did not transmit. Sites in pseudoautosomal regions and in bins excluded by the panel for
paralogy are dropped, and the member whose fraction is read must have DP ≥ `min_dp` at the site.

**Reference bias.** The member's genome-wide median of (alt fraction − 1/2) over every seventh confident autosomal
heterozygous site (≥ 200 sites) is subtracted from every site's alt fraction with the sign of the tag.

**Phased fraction.** Per site, the fraction of reads carrying the tagged allele: the maternal fraction along the
child, the transmitted-allele fraction along a parent. Expected values: 1/2 where the two homologues are present in
equal copy; (1 + f)/(2 + f) or 1/(2 + f) under a gain, 1/(2 − f) or (1 − f)/(2 − f) under a loss, and (1 ± f)/2 under
CN-LOH, according to which parent's copy is affected; 1 or 0 where only one parent's copies are present. A phasing
error (a genotype error elsewhere in the trio) flips one site's contribution and so attenuates a signal.

**Auxiliary tracks and homologue count.** The auxiliary classes are correctly phased only while each parent
contributed one homologue. Where the child carries two different homologues of one parent (meiotic trisomy,
heterodisomy), the class tagged by the other parent's homozygous genotype is misassigned: along a maternal meiotic
trisomy the `father_hom` track reads 1/3 where the two maternal homologues differ and returns to the main track's 2/3
where a crossover has made them identical. For a child's gain, LOH or UPD with a parent named and |shift| ≥ 0.05, the
share of the event's windows (≥ 4 windows with ≥ 5 auxiliary sites each) in which the relevant auxiliary track
disagrees in sign with the main track (for UPD: lies within half the main track's deviation) is `hetero_share`:
≥ 0.9 "two different homologues throughout (meiotic)", ≤ 0.1 "one homologue throughout (mitotic, or a meiosis II error
without a crossover)", otherwise the share is reported (meiotic with crossovers).

**Windows and step fit.** Main-class sites are pooled in consecutive windows of w sites, w chosen so that a window
spans about 400 kb at the member's genome-wide phased-site density (20 ≤ w ≤ 200); a window never spans a gap of more
than 3 Mb, and a trailing window less than half full is merged into the previous one. Each window carries its
depth-weighted pooled fraction and the binomial standard error at 1/2, √(1/(4 Σ depth)). A window is marked *shared*
when it deviates from 1/2 by ≥ 0.05 in two or more members at once (another member's window counts when its deviation
is ≥ 0.05 and at least half of this member's); shared windows indicate paralogous sequence or an imbalance common to
the family and are excluded from the fits and the scan. The step fit is the exact one-dimensional total-variation
denoising of the window fractions, min_x ½ Σ (y_i − x_i)² + λ Σ |x_{i+1} − x_i|, by Condat's direct algorithm, with
λ = 2.5 × the track's robust noise SD from its first differences (one λ per member). The same fit is applied to the
LRR. **Per-homologue copy number** is the LRR step fit's copy number, 2 × 2^LRR, multiplied by the fraction's step
fit, both fits taken over a running median of five points: maternal and paternal along the child, transmitted and
untransmitted along a parent. No copy number is drawn in windows the scan rejected.

**Per-event phased statistics.** The depth-weighted pooled fraction over the event's main-class sites gives the
shift from 1/2 (`phase_shift`), its binomial error (`phase_se`) and the site count (`n_phased`); f_phase follows from
|shift| by the band-deviation formulas above. With ≥ 40 sites and |shift| ≥ 3 SE, the sign gives `origin_phase`: for
the child, the parent of the extra, lost or retained copy (for UPD, the parent of both copies); for a parent, whether
the duplicated, lost or retained homologue is the transmitted one. A disagreement with the likelihood-based `origin`
is noted. For a depth-called event with ≥ 500 phased sites, f_phase < f_lrr/2 (and lower by more than 0.05) is noted
as a possible depth artefact; f_phase > 1.5 f_lrr is also noted.

**Boundaries at site resolution.** For each event boundary not at a chromosome end, the main-class sites within 1.5
bins of the bin edge are taken (the window is widened by factors of 1.5, to at most 6 bins, until at least 10 sites
lie on each side of the bin edge), signed by the event's shift, and split at the position (≥ 3 sites on each side)
that most reduces the squared error of a two-mean fit with the inside mean above the outside mean; the boundary is
the midpoint between the two sites flanking the split (`start_fine`, `end_fine`). `edge_sites` is the smaller of the
two site counts used.

**Phased scan.** On the main-class windows that lie outside the member's events and are not shared, the same binary
segmentation as for the depth is run (z = 5, ≥ 8 windows per segment; a chromosome with fewer than 16 windows is one
segment). A segment is reported when its mean shift is ≥ 0.015 (about 6% of cells for a gain or loss) and ≥ 5
empirical standard errors, its median shift is ≥ 0.015 with the same sign, ≥ 75% of its windows share the sign of the
mean, and it spans ≥ 2 Mb (a dense cluster of sites can fill many windows over a short distance). It is typed by the
depth over its bins: a gain when the mean LRR is ≥ 3 SE above zero, a loss when ≥ 3 SE below; otherwise CN-LOH with
f = 2d, requiring d ≥ 0.025 when the depth is flat (a gain or loss at a cell fraction the depth cannot resolve gives
the same track, and the note says so); or **UPD** (uniparental heterodisomy) when the depth is flat, d ≥ 0.4 and the
heterozygosity rate is ≥ 0.5 of the member's own, since the child is then homozygous for one parent's allele at every
informative site yet heterozygous wherever that parent is, a state that depth, folded BAF and heterozygosity rate all
miss. Segments failing the span or consistency criteria are written to `phased_rejected.tsv` with the reason.

## Parent of origin from informative sites

Informative sites are those where the father is 0/0 and the mother 1/1, or the reverse, with all three calls
confident. With f the event's cell fraction (f_lrr for a gain or loss, f_baf for LOH; 1 if neither is available;
clipped to [0.02, 1]), the child's alt count k at depth n is modelled as Binomial(n, p):

| event | model A | model B |
| --- | --- | --- |
| gain | extra copy maternal: p = (1+f)/(2+f) if the alt allele is maternal, 1/(2+f) otherwise | extra copy paternal (mirror image) |
| loss | paternal copy lost: p = 1/(2−f) if the alt allele is maternal, (1−f)/(2−f) otherwise | maternal copy lost (mirror image) |
| LOH | maternal copy retained: p = (1+f)/2 if the alt allele is maternal, (1−f)/2 otherwise | paternal copy retained (mirror image) |

p is clipped to [0.01, 0.99]. `origin_llr` = log L(A) − log L(B); its sign names the parent and its magnitude the
support; `origin_n` is the number of sites, at least 10. The ratio is conditional on the event's presence and type and
does not test for the event.

## Inheritance

A child's event is "inherited from the father/mother" when that parent has an event of the same type on the same
chromosome whose intersection with the child's covers at least half of the shorter of the two segments, with "(in a
share of the parent's cells)" when the parent's f < 0.8; otherwise it is "new (neither parent carries it)". A parent's
event is "passed to the child" or "not passed to the child" by the same rule. The rule is not reciprocal: a small
parental event inside a large event of the child satisfies it.

## Mendelian errors

Among confident sites in a region, a Mendelian error is a child homozygous for an allele absent from a parent, or
heterozygous where both parents are the same homozygote. `mie_rate` is the fraction of such sites (NA with fewer than
20 sites); the genome-wide rate is the median of the per-autosome rates.

## Runs of homozygosity versus uniparental disomy

A child's LOH with f ≥ 0.8 is a uniparental isodisomy only when it carries Mendelian errors at the informative sites
(`mie_rate` ≥ 0.05): the child is homozygous for one parent's allele where the parents are opposite homozygotes.
Without such errors it is a run of homozygosity (both copies identical by descent, the parents sharing the haplotype)
and no parent of origin is given. A parent's LOH with f ≥ 0.8 is annotated as a run of homozygosity unless the depth
indicates otherwise.

## Supplied events

`--events` are drawn and matched to triokaryo's events of the same sample and chromosome when the intersection covers
at least half of the shorter segment, irrespective of type. NGS-DOSE's `karyotype/events.tsv` is read by its columns:
`whole` spans the chromosome, `p`/`q` the arm (GRCh38 centromeres built in), and a stretch its `start_mb`–`end_mb`.
