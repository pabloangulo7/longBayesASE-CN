FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Point rstan at the conda toolchain, then prove at build time that the model
# actually compiles. A broken image fails here instead of in every ASE task.
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ASE_model.stan /tmp/ASE_model.stan
RUN mkdir -p ~/.R \
    && printf 'CXX14=%s\nCXX14FLAGS=-O2 -fPIC\n' "$(ls /opt/conda/bin/*-linux-gnu-g++ | head -n1)" > ~/.R/Makevars \
    && Rscript -e 'q(status = inherits(try(rstan::stan_model(file="/tmp/ASE_model.stan")), "try-error"))'

LABEL org.opencontainers.image.source="https://github.com/pabloangulo7/longBayesASE-CN" \
      org.opencontainers.image.description="longBayesASE-CN reproducible runtime"
