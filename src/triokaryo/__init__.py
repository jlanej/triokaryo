"""triokaryo - large chromosomal events in a trio, read from its VCF.

For child, father and mother: the depth of every PASS biallelic SNV, binned and normalised (LRR, GC-corrected where a
GC track is given), the B-allele fraction at heterozygous sites (BAF) and the heterozygosity rate; segmentation of
each; the copy state and the share of cells carrying each event from the depth and, independently, from the B-allele
bands; the parent of origin of the child's events from the sites where the parents are opposite homozygotes; what is
inherited and what is new; figures for review, and tables. No cohort is needed: every comparison is within the trio.
"""
__version__ = "0.1.0"
