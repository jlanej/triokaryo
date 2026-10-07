# triokaryo - container image for HPC (Docker build -> Apptainer/Singularity run).
#
#   docker build -t triokaryo .
#   docker run --rm -v "$PWD:/data" -w /data triokaryo run --vcf trio.vcf.gz --pedigree trios.tsv --child KID --out out/KID
#
# On HPC via Apptainer:
#   apptainer build triokaryo.sif docker://ghcr.io/jlanej/triokaryo:latest
#   apptainer run triokaryo.sif run --vcf trio.vcf.gz --pedigree trios.tsv --child KID --gc-track gc.tsv --out out/KID
ARG PY=3.12-slim
FROM python:${PY} AS runtime

LABEL org.opencontainers.image.title="triokaryo" \
      org.opencontainers.image.description="Large chromosomal events in a parent-offspring trio from its VCF: LRR, BAF and transmission phasing; cell fraction, parent of origin, inheritance" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends bcftools tabix && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/triokaryo
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir . && python -c "import pysam, numpy, matplotlib, triokaryo"

RUN useradd --create-home --uid 1000 runner
USER runner
WORKDIR /data
ENV MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl

ENTRYPOINT ["triokaryo"]
CMD ["--help"]
