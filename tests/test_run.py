"""Every planted event found with its type, share of cells, parent of origin and inheritance; nothing else called."""
import numpy as np


def find(events, sample, chrom, etype):
    return [e for e in events if e.sample == sample and e.chrom == chrom and e.type == etype]


def test_every_planted_event_is_found_and_nothing_else(run, mock):
    ev = run["events"]
    truth = mock["truth_data"]["events"]
    for t in truth:
        hits = [e for e in find(ev, t["sample"], t["chrom"], t["type"]) if min(e.end, t["end"]) - max(e.start, t["start"]) >= 0.5 * (t["end"] - t["start"])]
        assert hits, "planted %s not found: %s" % (t["label"], [(e.sample, e.chrom, e.type, e.start, e.end) for e in ev])
    # nothing beyond the planted ones (two bins of slack at the ends of a 4-Mb deletion is still one event)
    assert len(ev) == len(truth), [(e.sample, e.chrom, e.type, e.span, e.f) for e in ev]


def test_share_of_cells_from_depth_and_from_bands_agree(run):
    ev = run["events"]
    (e12,) = find(ev, "KID", "chr12", "gain")
    assert abs(e12.f_lrr - 0.30) < 0.05 and abs(e12.f_baf - 0.30) < 0.06 and e12.span == "whole"
    (e21,) = find(ev, "KID", "chr21", "gain")
    assert abs(e21.f_lrr - 1.0) < 0.08 and abs(e21.f_baf - 1.0) < 0.1
    (e8,) = find(ev, "MOM", "chr8", "gain")
    assert abs(e8.f_lrr - 0.15) < 0.05 and abs(e8.f_baf - 0.15) < 0.06
    (e6,) = find(ev, "KID", "chr6", "LOH")
    assert abs(e6.f_baf - 0.40) < 0.06 and e6.span == "p" and abs(e6.lrr) < 0.03
    (e7,) = find(ev, "KID", "chr7", "LOH")
    assert e7.f_baf == 1.0 and e7.span == "whole" and e7.het_rate_rel < 0.1 and "no heterozygous calls" in e7.note and "isodisomy" in e7.note
    (e18,) = find(ev, "KID", "chr18", "loss")
    assert abs(e18.f_lrr - 1.0) < 0.08 and e18.start == 55_000_000


def test_parent_of_origin(run):
    ev = run["events"]
    o = {(e.chrom, e.type): e for e in ev if e.role == "child"}
    assert o[("chr21", "gain")].origin == "extra copy maternal" and o[("chr21", "gain")].origin_llr > 50
    assert o[("chr12", "gain")].origin == "extra copy paternal" and o[("chr12", "gain")].origin_llr < -50
    assert o[("chr10", "gain")].origin == "extra copy paternal"
    assert o[("chr18", "loss")].origin == "paternal copy lost" and o[("chr18", "loss")].origin_llr > 50
    assert o[("chr2", "loss")].origin == "paternal copy lost"
    assert o[("chr7", "LOH")].origin.startswith("maternal copy retained") and o[("chr7", "LOH")].origin_n > 100
    assert o[("chr6", "LOH")].origin.startswith("maternal copy retained")


def test_inheritance_and_mendelian_errors(run):
    ev = run["events"]
    o = {(e.sample, e.chrom, e.type): e for e in ev}
    assert o[("KID", "chr10", "gain")].inheritance == "inherited from the father"
    assert o[("KID", "chr2", "loss")].inheritance == "inherited from the father"
    assert o[("KID", "chr21", "gain")].inheritance.startswith("new")
    assert o[("DAD", "chr10", "gain")].inheritance == "passed to the child"
    assert o[("MOM", "chr8", "gain")].inheritance == "not passed to the child"
    # a deletion and a uniparental disomy break Mendel at the informative sites; a trisomy and a mosaic do not
    assert o[("KID", "chr18", "loss")].mie_rate > 0.1 and o[("KID", "chr7", "LOH")].mie_rate > 0.1
    assert o[("KID", "chr21", "gain")].mie_rate < 0.01 and o[("KID", "chr12", "gain")].mie_rate < 0.01
    assert run["summary"]["mie_rate_genome"] < 0.01


def test_sex_chromosomes_and_sex_check(run, mock):
    s = run["summary"]
    assert s["child_x_copies"] == 1 and s["father_x_copies"] == 1 and s["mother_x_copies"] == 2
    assert s["child_x_check"] == "agrees" and s["mother_x_check"] == "agrees"
    assert (s["child_y_copies"], s["father_y_copies"], s["mother_y_copies"]) == (1, 1, 0) == tuple(mock["truth_data"]["y_copies"][k] for k in ("KID", "DAD", "MOM"))
    assert [s["%s_sex_karyotype" % r] for r in ("child", "father", "mother")] == ["XY", "XY", "XX"]
    assert all(s["%s_sex_check" % r] == "agrees" for r in ("child", "father", "mother"))
    assert abs(s["child_x_copies_raw"] - 1) < 0.1 and abs(s["mother_x_copies_raw"] - 2) < 0.1 and abs(s["father_y_copies_raw"] - 1) < 0.12 and s["mother_y_copies_raw"] == 0
    assert abs(s["y_father_son_log2"]) < 0.1 and s["y_father_son_sites"] > 500
    assert not [e for e in run["events"] if e.chrom in ("chrX", "chrY")]
    # the karyotype string: the constitutional events in the main line, each mosaic its own line
    import re
    k = s["child_karyotype"]
    main = k.split("/")[-1]
    assert k.startswith("mos ") and main.startswith("47,XY,") and main.endswith(",+21mat(MI)"), k            # one constitutional whole-chromosome gain
    # cytogenetic bands from the site-resolution boundaries: the inherited 10q gain's refined start lands in a centromeric band (p11.1 or q11.1)
    assert re.search(r"del\(2\)\(q11\.2q12\.[12]\)pat,upd\(7\)mat\(iso\),dup\(10\)\((p11\.1|q11\.1|q11\.2\d)q26\.3\)pat,upd\(15\)mat\(hetero\),del\(18\)\(q21\.[12]q23\)pat", main), main
    assert "/48,XY,+12pat(MII/mit)[0.30]/" in k and re.search(r"47,XY,loh\(6\)\(p25\.3p11\.1\)mat\[0\.[34]\d\]", k), k
    assert re.fullmatch(r"46,XY,del\(2\)\(q11\.2q12\.[12]\),dup\(10\)\((p11\.1|q11\.1|q11\.2\d)q26\.3\)", s["father_karyotype"]), s["father_karyotype"]
    (e2,) = [e for e in run["events"] if e.sample == "KID" and e.chrom == "chr2"]
    assert e2.bands.startswith("q11.2") and all(e.bands for e in run["events"])
    assert re.fullmatch(r"mos 47,XX,\+8\[0\.1\d\]/46,XX", s["mother_karyotype"]), s["mother_karyotype"]


def test_external_events_matched(run):
    ext = run["external"]
    assert len(ext) == 8 and all(x.inheritance.startswith("matched") for x in ext)
    assert all(e.external for e in run["events"] if e.type not in ("LOH", "UPD"))      # a depth tool sees neither
    assert all(not e.external for e in run["events"] if e.type == "LOH")


def test_outputs_written(run):
    import os
    out = run["out"]
    for f in ("events.tsv", "bins.tsv", "summary.tsv", "summary.json", "index.html", "external.tsv", "events.bed", "crossovers.tsv"):
        assert os.path.getsize(os.path.join(out, f)) > 0
    bed = open(os.path.join(out, "events.bed")).read().splitlines()
    assert bed[0].startswith("track name=") and len(bed) == 1 + len(run["events"]) and all(len(l.split("\t")) == 6 for l in bed[1:])
    assert any("KID|child|gain|p13.2q22.3|f=1.00|extra_copy_maternal" in l or ("KID|child|gain|" in l and "|f=1.00|extra_copy_maternal" in l) for l in bed[1:] if l.startswith("chr21"))
    import pysam
    vf = pysam.VariantFile(os.path.join(out, "events.vcf"))
    recs = list(vf)
    assert list(vf.header.samples) == ["KID", "DAD", "MOM"] and len(recs) == len(run["events"])
    r21 = [r for r in recs if r.chrom == "chr21"][0]
    assert r21.alts == ("<DUP>",) and r21.info["SVTYPE"] == "DUP" and r21.info["ORIGIN"] == "mat" and r21.info["STAGE"] == "MI" and abs(r21.info["CF"] - 1.0) < 0.1
    assert r21.samples["KID"]["GT"] == (0, 1) and abs(r21.samples["KID"]["CN"] - 3.0) < 0.2 and r21.samples["DAD"]["GT"] == (None, None)
    r7 = [r for r in recs if r.chrom == "chr7"][0]
    assert r7.alts == ("<CNLOH>",) and r7.samples["KID"]["GT"] == (1, 1) and r7.info["INHERITANCE"] == "new"
    xo = open(os.path.join(out, "crossovers.tsv")).read().splitlines()
    assert xo[0].split("\t")[:4] == ["trio", "sample", "role", "chrom"] and len(xo) == 1          # no crossover in the default simulation
    figs = run["figures"]
    names = {f["name"] for f in figs}
    assert "genome" in names and "chrom_chr21" in names and "chrom_chr7" in names
    for f in figs:
        for ext in ("png", "svg", "pdf"):
            assert os.path.exists(f[ext])
        assert os.path.exists(os.path.join(out, "figures", f["name"] + ".txt"))
        assert os.path.exists(os.path.join(out, "figures", "legends", f["name"] + "_legend.png"))
    html = open(os.path.join(out, "index.html")).read()
    assert "MOCK" not in html or "mock" in html.lower()
    assert "extra copy maternal" in html and "data:image/png;base64" in html
