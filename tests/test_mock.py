import pysam


def test_mock_vcf_reads_back(mock):
    vf = pysam.VariantFile(mock["vcf"])
    assert list(vf.header.samples) == ["KID", "DAD", "MOM"]
    n = sum(1 for _ in vf)
    assert n > 100000
    t = mock["truth_data"]
    assert len(t["events"]) == 10 and t["x_copies"] == {"KID": 1, "DAD": 1, "MOM": 2}
