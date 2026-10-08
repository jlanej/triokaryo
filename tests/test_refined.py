"""A VCF whose genotypes were refined under a pedigree prior (GATK CalculateGenotypePosteriors: PP beside PL), and one whose
homozygous-reference genotypes carry a reference block's depth (GATK GenotypeGVCFs: MIN_DP): the genotypes are re-derived
from PL, and the bin depth taken from each member's own variant genotypes."""
import json

import numpy as np

from triokaryo.mock import write_mock
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio
from triokaryo.vcfscan import GT_HET, GT_HOMALT, GT_HOMREF, _from_pl


def test_genotype_and_quality_from_pl():
    assert _from_pl([0, 45, 600]) == (GT_HOMREF, 45) and _from_pl([90, 0, 95]) == (GT_HET, 90) and _from_pl([400, 60, 0]) == (GT_HOMALT, 60)
    assert _from_pl([0, 120]) == (GT_HOMREF, 99) and _from_pl([30, 0]) == (GT_HOMALT, 30)
    assert _from_pl([0, 0, 50]) is None and _from_pl(None) is None and _from_pl([0, None, 50]) is None and _from_pl([0]) is None


def test_refined_genotypes_are_recovered_from_pl(tmp_path):
    """The family prior moves the child's Mendelian-violating genotypes to consistent ones or strips their quality: the isodisomy
    loses the Mendelian errors that identify it and the heterodisomy its informative sites. Read from PL, both are as planted."""
    contigs = ("chr7", "chr15", "chr16", "chr17", "chr18", "chrX", "chrY")
    m = write_mock(str(tmp_path / "refined"), seed=21, contigs=contigs, refined=True)
    trio = read_trios(m["trios"])[0]
    asis = run_trio(m["vcf"], trio, str(tmp_path / "asis"), gc_track=m["gc"], genotypes="vcf", figures=False, log=lambda s: None)
    auto = run_trio(m["vcf"], trio, str(tmp_path / "auto"), gc_track=m["gc"], figures=False, log=lambda s: None)
    assert auto["summary"]["genotypes_from_pl"] and not asis["summary"]["genotypes_from_pl"]
    e7 = [e for e in auto["events"] if e.sample == "KID" and e.chrom == "chr7" and e.type == "LOH"]
    assert len(e7) == 1 and "isodisomy" in e7[0].note and e7[0].mie_rate > 0.1 and e7[0].origin == "maternal copy retained (paternal replaced)", e7
    e7r = [e for e in asis["events"] if e.sample == "KID" and e.chrom == "chr7" and e.type == "LOH"]
    assert not e7r or e7r[0].mie_rate < 0.5 * e7[0].mie_rate or "isodisomy" not in e7r[0].note, [(e.mie_rate, e.note) for e in e7r]
    assert any(e.sample == "KID" and e.chrom == "chr15" and e.type == "UPD" for e in auto["events"])
    assert auto["summary"]["child_x_copies"] == 1 and auto["summary"]["child_sex_check"] == "agrees"
    found = {(e.sample, e.chrom, e.type) for e in auto["events"]}
    for t in json.load(open(m["truth"]))["events"]:
        if t["chrom"] in contigs:
            assert (t["sample"], t["chrom"], t["type"]) in found, t["label"]


def test_reference_block_depth_is_left_out_of_the_depth_track(tmp_path):
    """Homozygous-reference genotypes at a block's minimum depth read shallower than the heterozygous ones, and their share rises
    where heterozygous sites vanish (an isodisomy): taken over every genotype the depth dips there; taken over the member's own
    variant genotypes it is flat."""
    m = write_mock(str(tmp_path / "blocks"), seed=22, sites_per_mb=300, contigs=("chr7", "chr8", "chr9"), ref_blocks=True)
    trio = read_trios(m["trios"])[0]
    every = run_trio(m["vcf"], trio, str(tmp_path / "every"), gc_track=m["gc"], depth_sites="all", figures=False, log=lambda s: None)
    auto = run_trio(m["vcf"], trio, str(tmp_path / "auto"), gc_track=m["gc"], figures=False, log=lambda s: None)
    assert auto["summary"]["depth_sites"] == "variant" and every["summary"]["depth_sites"] == "all"
    assert 0.74 < auto["summary"]["child_homref_depth_ratio"] < 0.9, auto["summary"]["child_homref_depth_ratio"]
    sl = every["bins"].of("chr7")
    dip, flat = np.nanmean(every["bins"].lrr_gc[0][sl]), np.nanmean(auto["bins"].lrr_gc[0][sl])
    assert dip < -0.02 and abs(flat) < 0.015 and flat > dip + 0.015, (dip, flat)
