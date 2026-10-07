"""triokaryo batch: every trio of a pedigree from per-trio VCFs, in parallel, then the cohort report."""
import os

from triokaryo.cli import main
from triokaryo.mock import write_mock


def test_batch_runs_every_trio_and_the_cohort(tmp_path):
    rows = []
    for i, prefix in enumerate(("T0_", "T1_")):
        m = write_mock(str(tmp_path / ("%sKID" % prefix)), seed=60 + i, prefix=prefix, contigs=["chr20", "chr21", "chr22", "chrX", "chrY"], sites_per_mb=120,
                       no_events=(i == 1))
        rows += [l for l in open(m["trios"]).read().splitlines()[1:]]
    ped = tmp_path / "trios.tsv"
    ped.write_text("#kid\tdad\tmom\tkid_sex\tdad_sex\tmom_sex\n" + "\n".join(rows) + "\n")
    out = tmp_path / "out"
    rc = main(["batch", "--pedigree", str(ped), "--vcf-pattern", str(tmp_path / "{kid}" / "mock.vcf.gz"), "--out", str(out), "--jobs", "2", "--no-figures",
               "--cohort", str(tmp_path / "cohort")])
    assert rc == 0
    for kid in ("T0_KID", "T1_KID"):
        assert os.path.getsize(os.path.join(str(out), kid, "events.tsv")) > 0 and os.path.exists(os.path.join(str(out), kid, "summary.json"))
    assert "T0_KID" in open(os.path.join(str(tmp_path / "cohort"), "index.html")).read()
    assert len(open(os.path.join(str(tmp_path / "cohort"), "summary.all.tsv")).read().splitlines()) == 3
    # a trio whose VCF is missing is skipped with a warning, not an error
    ped2 = tmp_path / "trios2.tsv"
    ped2.write_text("#kid\tdad\tmom\nNOPE\tA\tB\n" + rows[0] + "\n")
    assert main(["batch", "--pedigree", str(ped2), "--vcf-pattern", str(tmp_path / "{kid}" / "mock.vcf.gz"), "--out", str(tmp_path / "out2"), "--no-figures"]) == 0
