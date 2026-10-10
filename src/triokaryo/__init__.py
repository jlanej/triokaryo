"""triokaryo: large chromosomal events in a parent-offspring trio from its small-variant VCF.

For child, father and mother: binned read depth (LRR, GC- and panel-corrected where tracks are given), B-allele
frequency at heterozygous sites and heterozygosity rate, each segmented; transmission phasing of the child's alleles,
giving a signed allelic-imbalance track; per event the type, the cell fraction (from depth, from the folded BAF bands
and from the phased track), the parent of origin, the homologue count, site-resolution boundaries, the Mendelian-error
rate and the inheritance; tables and figures. All comparisons are within the trio; no cohort is required.
"""
__version__ = "0.1.4"
