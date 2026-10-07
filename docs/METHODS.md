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

**Transmission phasing: the main tracks.** The alt-allele fraction at a heterozygous site sits at 1/2 ± d with
the sign unknown, so on its own it is read folded (|BAF − 1/2|, biased upward by the noise: 0.07 at 30× where
d = 0). The trio gives the sign. At a site where the parents are opposite homozygotes the child's two alleles have
known parents whatever the child's call, and the fraction of the mother's allele is the maternal fraction: 1/2 where
the homologues are equal, 1/2 + d or 1/2 − d along an event by the parent of origin, 1 where no paternal copy is
left (a deletion of the paternal copy, a maternal isodisomy or heterodisomy). At a parent's heterozygous site where
the child is homozygous, the allele the parent passed on is the child's, and the fraction of the transmitted allele
along the parent says whether an event of the parent's lies on the homologue the child got. These are the main
tracks. A phasing error (a genotype error elsewhere in the trio) flips a site at random and so weakens a reading
rather than inventing one. The reference bias (the alt allele read a little under one half) is taken out by the
member's own genome-wide median of (alt fraction − 1/2) over its heterozygous sites, with the sign of each site's
tag. Sites in bins whose bands are parted in the panel's genomes are left out.

**The auxiliary tracks, and one homologue or two.** The sites where only one parent is homozygous phase the child
too - the homozygous parent's allele is theirs, the other allele the other parent's - but only while each parent gave
one homologue. Where the child carries two different homologues of one parent (a meiotic trisomy, a heterodisomy)
the sites tagged by the other parent's homozygosity read the wrong way. So they are pooled as auxiliary tracks beside
the main one: along a maternal meiotic trisomy the track read at the father's homozygous sites leaves the main track
(1/3 against 2/3) wherever the two maternal copies differ and returns wherever a crossover made them identical - a
map of the meiotic error and its crossovers. A parent's auxiliary track (the parent and the child heterozygous, the
other parent homozygous) does the same from the parent's side. For each gain, LOH or heterodisomy of the child's with
a parent named and a shift of at least 0.05, the share of the event's windows where the auxiliary track parts from
the main one is reported: ≥ 90% "two different homologues throughout (meiotic)", ≤ 10% "one homologue throughout
(mitotic, or a meiosis II error without a crossover)", else the share (a meiotic error with crossovers).

**Pooled windows, the step fit, the copies.** The main track's sites are pooled in windows of a fixed number of
sites (the number set so that a window spans about 400 kb at the genome's density, never fewer than 20 sites nor
more than 200; a window never crosses a gap of over 3 Mb), the fraction depth-weighted, with its binomial error. A
window shifted by ≥ 0.05 in two or more members at once, each at least half the other's, is parted in everyone -
paralogous sequence, or an imbalance the family shares - and is left out of the fits and the scan. The step fit is
the exact one-dimensional total-variation denoising of the windows (min ½ Σ (y − x)² + λ Σ |x_{i+1} − x_i|, Condat's
direct algorithm; λ = 2.5 times the track's noise sd from its first differences): a piecewise-constant reading in
which every jump has to earn its height - a lone spike under 2λ is flattened, a plateau of m points keeps its height
less 2λ/m. The same fit is drawn on the LRR. The two homologues' copies are the LRR step fit's copies (2 × 2^LRR)
split by the fraction's step fit: maternal and paternal along the child, passed and not passed along a parent.

**Each event's phased reading.** The pooled fraction over the event's main-track sites gives the shift from one
half, its binomial error, and the share of cells by the event's type (gain 4d/(1 − 2d), loss 4d/(1 + 2d), LOH 2d).
With ≥ 40 sites and a shift of ≥ 3 errors, the sign names the parent: for the child the extra copy's parent, the
lost copy's, the retained copy's; for a parent whether the duplicated, lost or retained homologue is the one passed
to the child. The opposite-homozygote reading above is a subset of the same sites read with a likelihood ratio; the
two are reported side by side and a disagreement noted. For a depth event with ≥ 500 phased sites, a phased share
under half the depth's is noted (the depth's call may be a GC or panel artefact rather than a copy change).

**Edges at site resolution.** For an event not running to the chromosome's end, the sites of the main track within
1.5 bins of each bin edge (wider where sparse), signed by the event's shift, are split where a two-mean fit (an
outside at one half, an inside away from it) reduces the squared error most, with at least 12 sites each side; the
edge is the midpoint between the two sites either side of the split.

**The phased scan.** On the main track's windows not already inside an event of the member's (and not parted in
everyone), the same recursive segmentation as the LRR's (z = 5, 8 windows a segment; a chromosome too short to split
is one segment). A segment is an event when its mean shift is ≥ 0.015 (a gain or loss in about 6% of cells, an LOH
in 3%) and ≥ 5 empirical standard errors, its median shift is ≥ 0.015 the same way and ≥ 75% of its windows lean
that way (a few parted windows do not make an event), and it spans ≥ 2 Mb (a dense cluster of sites, as a
gene-family region, makes many windows of a few hundred kb). It is typed by the depth's lean over the same bins
(≥ 3 standard errors up, a gain; down, a loss), else read as a copy-neutral LOH at a share of 2d - which a gain or
loss at a share the depth cannot resolve would mimic, and the note says so - unless the shift is ≥ 0.4 with the
heterozygosity rate kept (≥ 0.5 of the member's own): a uniparental heterodisomy, which the depth, the folded bands
and the heterozygosity rate all miss (the child homozygous for one parent's allele wherever the parents are opposite
homozygotes, heterozygous wherever that parent is). Segments set aside, with the reason, go to `phased_rejected.tsv`.

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
