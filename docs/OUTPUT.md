# Output

## events.tsv (per trio) / events.all.tsv (cohort)

| column | meaning |
| --- | --- |
| sample, role | the member (child, father, mother) |
| chrom, start, end | the segment (bp; start 0-based bin start, end the bin end) |
| span | whole, p, q, stretch |
| type | gain, loss, LOH (copy-neutral loss of heterozygosity) |
| f | the share of cells: `f_lrr` for a gain or loss, else `f_baf` |
| f_lrr, f_baf | from the depth; from the B-allele bands |
| lrr, lrr_se, n_bins | the segment's mean LRR, its standard error, its bins |
| d_hat, llr_baf | the band deviation at the segment's heterozygous sites and its evidence against zero |
| het_rate, het_rate_rel, n_het, n_called | the heterozygosity rate in the segment, relative to the member's own, and the counts |
| mie_rate | Mendelian errors among the segment's confident sites (summary: the genome's baseline) |
| origin, origin_llr, origin_n | the child's parent of origin (see METHODS), its log-likelihood ratio, the informative sites |
| inheritance | the child's: inherited from whom, or new; a parent's: passed to the child or not |
| external | the labels of the given events (`--events`) the segment overlaps |
| note | a single X, no heterozygous calls (in every cell), ... |

## bins.tsv

chrom, start, end, gc; per member `n_sites`, `depth`, `lrr`, `lrr_gc`, `n_called`, `n_het`, `het_rate`, `bdev`;
`child_vs_mid`, `father_vs_mother`.

## summary.tsv / summary.json

The trio and its members and sexes; records and sites used, the skipped records by reason; the genome's
Mendelian-error rate; per member the events by type, the X copies and the check against the pedigree's sex, the
median autosomal depth; the parameters.

## figures/

`genome.{png,svg,pdf}` and `chrom_<chrom>.{png,svg,pdf}`, each with `<name>.txt` (title, caption, key with hex
colours) and `legends/<name>_legend.{png,svg,pdf}`. Colours (Okabe-Ito): LRR points dark grey; BAF light grey;
gain vermillion `#D55E00`; loss blue `#0072B2`; LOH reddish purple `#CC79A7`; in the child's BAF panel the
informative sites whose alt allele came from the mother `#CC79A7`, from the father `#E69F00`; given events black
brackets; references dotted grey.

## cohort/

`events.all.tsv`, `summary.all.tsv`, `concordance.tsv` (each given event with its match), `index.html`.
