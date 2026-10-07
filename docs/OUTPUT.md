# Output

## events.tsv (per trio) / events.all.tsv (cohort)

| column | meaning |
| --- | --- |
| sample, role | the member (child, father, mother) |
| chrom, start, end | the segment (bp; start 0-based bin start, end the bin end) |
| span | whole, p, q, stretch |
| type | gain, loss, LOH (copy-neutral loss of heterozygosity), UPD (a uniparental heterodisomy) |
| source | depth (the LRR), bands (the folded bands or the heterozygosity rate), phased (the phased scan: what the depth did not call) |
| f | the share of cells: `f_lrr` for a gain or loss, else `f_baf` |
| f_lrr, f_baf | from the depth; from the B-allele bands |
| lrr, lrr_se, n_bins | the segment's mean LRR, its standard error, its bins |
| d_hat, llr_baf | the band deviation at the segment's heterozygous sites and its evidence against zero |
| het_rate, het_rate_rel, n_het, n_called | the heterozygosity rate in the segment, relative to the member's own, and the counts |
| mie_rate | Mendelian errors among the segment's confident sites (summary: the genome's baseline) |
| origin, origin_llr, origin_n | the child's parent of origin from the opposite-homozygote sites (see METHODS), its log-likelihood ratio, the sites |
| phase_shift, phase_se, n_phased | the phased fraction's shift from one half over the event (child: the maternal allele; parent: the allele passed to the child), its binomial error, the sites |
| f_phase, origin_phase | the share of cells from the shift; what its sign says (the child: the extra, lost or retained copy's parent; a parent: whether the event lies on the homologue passed to the child) |
| homologues, hetero_share | a child's gain, LOH or heterodisomy: the named parent's two copies one homologue or two, and the share of the event's windows where they differ |
| start_fine, end_fine, edge_sites | the edges at site resolution from the phased sites (NA at a chromosome's end), and the sites the narrower edge used |
| inheritance | the child's: inherited from whom, or new; a parent's: passed to the child or not |
| external | the labels of the given events (`--events`) the segment overlaps |
| note | a single X, no heterozygous calls (in every cell), ... |

## phased.tsv

The pooled windows of each member's main phased track: role, chrom, start, end, mid, n_sites, depth (the summed
depth), frac (the pooled fraction), se, shared (parted in two or more members), step (the step fit), lrr (the bin's),
copies_tagged, copies_other (from the step fits), and the auxiliary tracks' pooled fractions and sites per window
(`aux_mother_hom`, `aux_father_hom` for the child; `aux_child_het` for a parent).

## phased_rejected.tsv

Segments of the phased scan set aside, with the reason (under the 2-Mb floor; the windows not agreeing with the mean).

## bins.tsv

chrom, start, end, gc; per member `n_sites`, `depth`, `lrr`, `lrr_gc`, `n_called`, `n_het`, `het_rate`, `bdev`;
`child_vs_mid`, `father_vs_mother`.

## summary.tsv / summary.json

The trio and its members and sexes; records and sites used, the skipped records by reason; the genome's
Mendelian-error rate; per member the events by type, the X copies and the check against the pedigree's sex, the
median autosomal depth; the parameters.

## figures/

`genome.{png,svg,pdf}` (eight rows: LRR with its step fit and BAF with the pooled phased fraction per member, the
child over the parents' mean, the child's maternal and paternal copies) and `chrom_<chrom>.{png,svg,pdf}` (five rows
per member: LRR with its step fit and calls, the BAF with the informative sites coloured, the phased fraction - sites,
windows, step fit, auxiliary tracks - the two homologues' copies, the heterozygosity rate), each with `<name>.txt`
(title, caption, key with hex colours) and `legends/<name>_legend.{png,svg,pdf}`. Colours (Okabe-Ito): LRR points
dark grey; BAF light grey; the pooled phased fraction bluish green `#009E73` (hollow where parted in everyone); step
fits black; gain vermillion `#D55E00`; loss blue `#0072B2`; LOH and UPD reddish purple `#CC79A7`; in the child's
panels the mother's alleles `#CC79A7`, the father's `#E69F00` (the informative sites; the copies; the auxiliary track
that parts where the child carries two different homologues of that parent); a parent's copies passed `#D55E00` and
not passed `#56B4E9`; given events black brackets; references dotted grey.

## cohort/

`events.all.tsv`, `summary.all.tsv`, `concordance.tsv` (each given event with its match), `index.html`.
