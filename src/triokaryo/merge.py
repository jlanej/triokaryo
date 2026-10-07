"""Per-sample VCFs into one trio VCF with bcftools: each sample's file reduced to PASS biallelic SNVs, the three merged with
`bcftools merge -0` (a sample without a record at a site is written as homozygous reference, which the trio analyses take
as a confident parental genotype), indexed. Mirrors the example recipe (docs/example/fetch_1kg_trio.sh) as one command;
needs bcftools and tabix on the PATH (the container image has them)."""
import os
import shutil
import subprocess


def _run(cmd, log=None):
    if log:
        log(" ".join(cmd))
    subprocess.run(cmd, check=True)


def merge_trio(vcfs, out, log=None, keep_intermediate=False, threads=1):
    """vcfs: the child's, the father's and the mother's VCF (bgzipped and indexed, or plain); out: the trio VCF to write (.vcf.gz).
    Returns the output path. Raises SystemExit with the reason when bcftools is missing or a step fails."""
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools is not on the PATH; install it (conda install -c bioconda bcftools) or use the container image")
    if len(vcfs) != 3:
        raise SystemExit("merge takes exactly three VCFs: the child's, the father's and the mother's")
    out = out if out.endswith(".vcf.gz") else out + ".vcf.gz"
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    tmp = []
    try:
        for i, v in enumerate(vcfs):
            if not os.path.exists(v):
                raise SystemExit("not found: %s" % v)
            snv = out[:-7] + ".%d.snv.vcf.gz" % i
            _run(["bcftools", "view", "--threads", str(threads), "-f", "PASS,.", "-v", "snps", "-m2", "-M2", "-Oz", "-o", snv, v], log)
            _run(["bcftools", "index", "-t", "-f", snv], log)
            tmp.append(snv)
        _run(["bcftools", "merge", "--threads", str(threads), "-0", "-m", "none", "-Oz", "-o", out] + tmp, log)
        _run(["bcftools", "index", "-t", "-f", out], log)
    except subprocess.CalledProcessError as e:
        raise SystemExit("bcftools failed (exit %d): %s" % (e.returncode, " ".join(e.cmd)))
    finally:
        if not keep_intermediate:
            for p in tmp:
                for q in (p, p + ".tbi"):
                    if os.path.exists(q):
                        os.remove(q)
    return out
