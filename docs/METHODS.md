# Methods

## Sites

One pass over the VCF with pysam, the three members selected. A record is used when it is a biallelic SNV whose
FILTER is PASS or empty. Per member: `DP` (else the sum of `AD`), `AD[1]`, the genotype class (0/0, 0/1, 1/1,
missing; a haploid call counts as the homozygote of its allele) and `GQ`. Records on contigs the genome does not
carry (alts, decoys, chrM) are skipped. `--thin n` keeps every n-th usable record.

## Bins and tracks

Fixed-width bins (default 1 Mb) from the chromosome start. Per bin and member: the number of sites with depth,
the **depth** as the median DP over them, confident calls (DP ≥ `--min-dp`, GQ ≥ `--min-gq`), heterozygous calls
among them, the **heterozygosity rate**, and the **band deviation** as the median |BAF − ½| at the heterozygous
sites (five or more). **LRR** = log2(depth / the member's median depth over autosomal bins with ≥ 20 sites).

**GC correction.** With a GC track (`gc-track`: the GC fraction of each bin from the reference), the LRR of the
autosomal bins is sorted by GC and a running median (41 bins) gives the LRR expected at each GC; it is subtracted
from every bin by interpolation, and the result is re-centred so the autosomal median is zero. A median is robust
to the events themselves. Without a track the LRR is used as it is.

**The panel.** `triokaryo panel` reads other genomes the same way (per-sample or multi-sample VCFs, or earlier
runs' `bins.tsv`) and writes, per bin, the number of genomes, the median self-normalised LRR (the X aligned to two
copies by each genome's own X median) and its robust SD, the median band deviation, and the median heterozygosity
rate. With `--panel`, each member's LRR is the self-normalised LRR less the panel's median, re-centred on the
autosomes; bins with fewer than five genomes, no median, or a spread beyond 0.25 are masked (no call); the band
deviation is read net of the panel's regional excess over its typical bin, and bins whose panel bands are parted by
more than 0.03 above typical are left out of the LOH search; the heterozygosity rate is read as the ratio to the
panel's. The GC correction then acts on the panel-corrected LRR.

**Within-trio tracks.** Per site with depth in all three, log2 of the child's depth over the parents' mean, and
log2 of the father's over the mother's; the bin's value is the median over its sites. They cancel what the three
libraries share at a site and are drawn, not called.

## Segmentation and calls

Per member and chromosome, over bins with ≥ `min_sites` (20) sites: the noise scale is 1.4826 × MAD of the first
differences / √2. Binary segmentation splits recursively at the position maximising |mean(left) − mean(right)| /
(sd √(1/nL + 1/nR)) while that statistic is ≥ `z` (5) and both sides hold ≥ `min_len` (5) bins; adjacent segments
whose means differ by less than three standard errors are merged. A segment with |mean LRR| ≥ max(0.07, 3 SE) is
a **gain** (mean > 0) or **loss**, with `f_lrr` = 2(2^LRR − 1) or 2(1 − 2^LRR); events with f < `min_f` (0.10)
are not reported.

The **band deviation** d of a segment is the maximum-likelihood value over a grid (0 to 0.475 by 0.005) of
Σ log[½ Binom(alt; dp, ½ + d) + ½ Binom(alt; dp, ½ − d)] over the member's heterozygous sites in the segment, with
`llr_baf` the log-likelihood gain over d = 0. `f_baf` = 4d/(1 − 2d) for a gain, 4d/(1 + 2d) for a loss, 2d for
an LOH.

**Copy-neutral LOH** is searched where the depth called nothing: the band-deviation track (bins with ≥ 5
heterozygous calls; others at the member's autosomal median) and the heterozygosity-rate track (bins with ≥
`min_called` confident calls — 50, or half the member's typical bin where the VCF is sparser) are segmented the
same way. A segment is an LOH when its site-level d ≥ 0.04 with `llr_baf` ≥ 10 (a share of cells of about 8% and
up; f = 2d), or when its heterozygosity rate is ≤ 0.35 of the member's own with enough confident calls (in every
cell; f = 1, noted "no heterozygous calls").

After the recursive split each boundary is moved within ±4 bins to where the two neighbours' squared deviations
from their means are least, and remnants shorter than `min_len` bins are not called. Adjacent events of one member
and type parted by at most three masked bins (a centromere) or by a step in level of at most 0.08 (LOH: 0.15 in
f) are joined into one, and the span re-read.

**Span.** A segment covering ≥ 90% of the chromosome's usable bins is `whole`; ≥ 90% of one arm's and under half
of the other's, `p` or `q`; else `stretch`.

**Runs of homozygosity.** A child's LOH in every cell is a uniparental isodisomy only where it breaks Mendel at
the informative sites (the child homozygous for one parent's allele: a Mendelian-error rate ≥ 0.05 in the segment);
without errors it is a run of homozygosity, both copies identical by descent, and no parent of origin is given.

**The X.** Bins in the pseudoautosomal regions are left out. The member's X copy state is round(2 × 2^median LRR)
over the X (1 or 2, or 3), reported with the pedigree's sex; the X's LRR is read against that median, so a
member's X events are changes of what the member has. A single X has no heterozygous sites, so no LOH is sought
there. The Y is not read.

## The trio

**Informative sites:** father 0/0 and mother 1/1, or the reverse, all three confident. With *f* the event's share
of cells (the depth's for a gain or loss, the bands' for an LOH; 1 where neither), the child's alt count k at
depth n follows Binomial(n, p) with p by model:

| event | model A | model B |
| --- | --- | --- |
| gain | extra copy maternal: p = (1+f)/(2+f) if alt is maternal, 1/(2+f) otherwise | extra copy paternal (the mirror) |
| loss | paternal copy lost: p = 1/(2−f) if alt is maternal, (1−f)/(2−f) otherwise | maternal copy lost (the mirror) |
| LOH | maternal copy retained: p = (1+f)/2 if alt is maternal, (1−f)/2 otherwise | paternal copy retained (the mirror) |

p is kept within [0.01, 0.99]. `origin_llr` = log L(A) − log L(B); its sign names the parent, its size says how
sure; `origin_n` the sites. At least ten informative sites are needed.

**Inheritance:** a child's event is "inherited from the father/mother" when a parent has an event of the same
type on the same chromosome with reciprocal overlap ≥ 0.5, "(in a share of the parent's cells)" when that parent's
f < 0.8; else "new (neither parent carries it)". A parent's event is "passed to the child" or "not passed".

**Mendelian errors:** among confident sites in a region, the child homozygous for an allele no parent carries, or
heterozygous where both parents are the same homozygote; the genome's rate is the median over the autosomes.

## Events given from elsewhere

`--events` is drawn and matched (same sample, reciprocal overlap ≥ 0.5). NGS-DOSE's `karyotype/events.tsv` is
read by its columns: `whole` spans the chromosome, `p`/`q` the arm (GRCh38 centromeres built in), a stretch its
`start_mb`–`end_mb`.
