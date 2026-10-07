"""A GC track: the GC fraction of each bin of the genome, from the reference FASTA (pysam.FastaFile), written as
chrom, start (0-based), end, gc. Built once per reference and bin size; `run --gc-track` corrects the LRR with it."""
import pysam


def write_gc_track(fasta, out, genome, bin_size=1_000_000):
    fa = pysam.FastaFile(fasta)
    names = set(fa.references)
    n = 0
    with open(out, "w") as fh:
        fh.write("chrom\tstart\tend\tgc\n")
        for chrom in genome.chroms:
            name = chrom if chrom in names else (chrom[3:] if chrom[3:] in names else None)
            if name is None:
                continue
            L = fa.get_reference_length(name)
            for s in range(0, L, bin_size):
                e = min(s + bin_size, L)
                seq = fa.fetch(name, s, e).upper()
                at = seq.count("A") + seq.count("T")
                gc = seq.count("G") + seq.count("C")
                fh.write("%s\t%d\t%d\t%s\n" % (chrom, s, e, ("%.4f" % (gc / (at + gc))) if at + gc else "NA"))
                n += 1
    fa.close()
    return n
