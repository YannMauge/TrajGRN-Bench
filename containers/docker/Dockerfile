FROM mambaorg/micromamba:1.5.8

ARG ENV_FILE
ARG ENV_NAME=method_env

COPY --chown=$MAMBA_USER:$MAMBA_USER ${ENV_FILE} /tmp/environment.yml

RUN sed -i '/^name:/d' /tmp/environment.yml \
    && micromamba create -y -n ${ENV_NAME} -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV MAMBA_DEFAULT_ENV=${ENV_NAME}
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH

WORKDIR /work
ENTRYPOINT []
