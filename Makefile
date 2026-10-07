# Developer workflow. `make venv` once, then `make test`.
VENV := .venv
PY   := $(VENV)/bin/python
PIP  := $(VENV)/bin/pip
IMAGE ?= triokaryo:local

.PHONY: venv test demo docker docker-test clean

venv:                     ## a virtualenv with the package and the test deps
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -e '.[test]'

test: venv                ## the test suite (writes its mock trio under the test's temp dir)
	$(PY) -m pytest -q

demo: venv                ## the mock trio, run end to end into ./mock_out
	$(PY) -m triokaryo.cli mock --out mock_out/mock
	$(PY) -m triokaryo.cli run --vcf mock_out/mock/mock.vcf.gz --pedigree mock_out/mock/mock.trios.tsv --child KID \
	  --gc-track mock_out/mock/gc.tsv --events mock_out/mock/events.external.tsv --out mock_out/KID
	$(PY) -m triokaryo.cli cohort --runs mock_out/KID --events mock_out/mock/events.external.tsv --out mock_out/cohort
	@echo "open mock_out/KID/index.html"

docker:                   ## the container image
	docker build -t $(IMAGE) .

docker-test: docker       ## the mock, inside the image
	docker run --rm -v "$$PWD:/data" -w /data $(IMAGE) mock --out /data/ci_out/mock
	docker run --rm -v "$$PWD:/data" -w /data $(IMAGE) run --vcf /data/ci_out/mock/mock.vcf.gz --pedigree /data/ci_out/mock/mock.trios.tsv \
	  --child KID --gc-track /data/ci_out/mock/gc.tsv --out /data/ci_out/KID --no-figures
	@grep -q "trisomy\|gain" ci_out/KID/events.tsv && echo "image OK"

clean:
	rm -rf $(VENV) .pytest_cache build dist *.egg-info src/*.egg-info ci_out mock_out
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
