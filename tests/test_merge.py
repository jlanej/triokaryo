"""triokaryo merge: per-sample VCFs into a trio VCF with bcftools (skipped without bcftools)."""
import os
import shutil

import pytest


def test_merge_per_sample_vcfs(tmp_path):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools not on the PATH")
    import pysam
    from triokaryo.cli import main
    from triokaryo.mock import write_mock
    from triokaryo.pedigree import read_trios
    from triokaryo.pipeline import run_trio
    m = write_mock(str(tmp_path / "mock"), seed=2, contigs=["chr20", "chr21", "chr22"], sites_per_mb=200)
    # split the simulated trio VCF into per-sample files, dropping each sample's homozygous-reference records as a per-sample caller would
    per = []
    for sample in ("KID", "DAD", "MOM"):
        p = str(tmp_path / ("%s.vcf.gz" % sample))
        vf = pysam.VariantFile(m["vcf"])
        vf.subset_samples([sample])
        out = pysam.VariantFile(p, "wz", header=vf.header)
        for rec in vf:
            gt = rec.samples[sample]["GT"]
            if gt is not None and any(a for a in gt if a):
                out.write(rec)
        out.close()
        pysam.tabix_index(p, preset="vcf", force=True)
        per.append(p)
    out = str(tmp_path / "trio.vcf.gz")
    assert main(["merge", "--child", per[0], "--father", per[1], "--mother", per[2], "--out", out]) == 0
    assert os.path.exists(out) and os.path.exists(out + ".tbi") and not os.path.exists(str(tmp_path / "trio.0.snv.vcf.gz"))
    vf = pysam.VariantFile(out)
    assert list(vf.header.samples) == ["KID", "DAD", "MOM"]
    n = sum(1 for _ in vf)
    assert n > 1000
    # the merged trio runs, and the merge's homozygous-reference convention gives confident parental genotypes
    trio = read_trios(m["trios"])[0]
    res = run_trio(out, trio, str(tmp_path / "run"), figures=False, log=lambda s: None)
    assert res["summary"]["child_phased_sites"] > 100 and res["summary"]["mie_rate_genome"] < 0.02
