FROM mambaorg/micromamba:2.3.3@sha256:e583f8e151fd1ebbc78b05cf4feb37ddbb518892dbabf6c0c03ea8ad613bc519

USER root
COPY docs/environment-linux-64.lock /tmp/absolve-linux-64.lock
RUN micromamba create --yes --no-rc --prefix /opt/absolve-env --file /tmp/absolve-linux-64.lock \
    && micromamba clean --yes --all

ENV PATH="/opt/absolve-env/bin:${PATH}" \
    PYTHONNOUSERSITE=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1
WORKDIR /opt/absolve
COPY *.py ./
COPY runtime/ runtime/
COPY workers/ workers/
COPY panel_report/ panel_report/
COPY examples/ examples/
LABEL org.opencontainers.image.source="https://github.com/JasonLinjc/ABSOLVE-seq-agent"
USER $MAMBA_USER
ENTRYPOINT []
CMD ["python", "/opt/absolve/paired.py", "--help"]
