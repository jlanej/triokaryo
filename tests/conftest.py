import json
import os

import pytest

from triokaryo.mock import write_mock
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio


@pytest.fixture(scope="session")
def mock(tmp_path_factory):
    d = tmp_path_factory.mktemp("mock")
    paths = write_mock(str(d / "mock"), seed=1)
    paths["truth_data"] = json.load(open(paths["truth"]))
    paths["dir"] = str(d)
    return paths


@pytest.fixture(scope="session")
def run(mock):
    trio = read_trios(mock["trios"])[0]
    out = os.path.join(mock["dir"], "KID")
    res = run_trio(mock["vcf"], trio, out, gc_track=mock["gc"], events_path=mock["events"], figures=True, log=lambda s: None)
    res["out"] = out
    res["trio"] = trio
    return res
