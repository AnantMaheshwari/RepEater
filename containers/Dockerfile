# RepEater — single image
FROM condaforge/miniforge3:26.5.3-0@sha256:66f4b7d5dd07883181ef5370f44b43857008a61e8f5b9bcc07f3f3ee027f7e63


ARG GIT_COMMIT=unknown
ARG REP_EATER_VERSION=0.1.0

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]


RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential \
      zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

COPY containers/environment.yml /tmp/environment.yml
RUN mamba env update -n base -f /tmp/environment.yml \
 && mamba clean --all --yes \
 && rm /tmp/environment.yml

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

COPY . /opt/RepEater
WORKDIR /opt/RepEater
RUN make

RUN pip install --no-cache-dir --no-deps -e /opt/RepEater

RUN TESORTER_DB="$(python -c 'import TEsorter, os; print(os.path.join(os.path.dirname(TEsorter.__file__), "database"))')" \
 && echo "TEsorter database dir: ${TESORTER_DB}" \
 && for hmm in "${TESORTER_DB}"/*.hmm; do \
      echo "hmmpress ${hmm##*/}"; \
      hmmpress -f "$hmm"; \
    done \
 && ls -1 "${TESORTER_DB}"/*.h3i | wc -l | xargs echo "pressed databases:"


RUN TEsorter --version \
 && hmmsearch -h | head -2 \
 && hmmpress -h > /dev/null \
 && blastn -version \
 && blastp -version \
 && makeblastdb -version \
 && RepEater --help > /dev/null \
 && test -x /opt/RepEater/external/FASTAN/FasTAN \
 && test -x /opt/RepEater/external/FASTAN/FastLTR \
 && test -x /opt/RepEater/external/alntools/taco \
 && test -x /opt/RepEater/seq/extract_ltrs \
 && test -x /opt/RepEater/seq/aln_info \
 && test -f /opt/RepEater/algorithms/kmer_seed.so \
 && test -f /opt/RepEater/algorithms/wfa_align.so \
 && python -c 'import numpy, scipy, matplotlib, networkx, Bio; print("python deps OK")'

RUN chmod -R a+rX /opt/RepEater /opt/conda

ENV PYTHONPATH=/opt/RepEater \
    KMP_DUPLICATE_LIB_OK=TRUE \
    REP_EATER_IN_CONTAINER=1 \
    REP_EATER_VERSION=${REP_EATER_VERSION} \
    REP_EATER_GIT_COMMIT=${GIT_COMMIT}

LABEL org.opencontainers.image.title="RepEater" \
      org.opencontainers.image.description="De novo LTR retrotransposon discovery, clustering and classification" \
      org.opencontainers.image.source="https://github.com/AnantMaheshwari/RepEater" \
      org.opencontainers.image.version="${REP_EATER_VERSION}" \
      org.opencontainers.image.revision="${GIT_COMMIT}" \
      io.rep_eater.tesorter="1.5.1" \
      io.rep_eater.hmmer="3.4" \
      io.rep_eater.rmblast="2.17.0"

WORKDIR /data
ENTRYPOINT ["RepEater"]
CMD ["--help"]
