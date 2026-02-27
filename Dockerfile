# syntax=docker/dockerfile:1.7

FROM mambaorg/micromamba:1.5.10 AS builder

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG ENV_NAME=cryptic-ip

WORKDIR /workspace

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt setup.py README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER cryptic_ip ./cryptic_ip

RUN micromamba install -y -n ${ENV_NAME} -c conda-forge \
    python=3.10 \
    pip \
    fpocket \
    freesasa \
    apbs \
    pdb2pqr \
    openmm \
    && micromamba clean --all --yes

RUN pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -e .

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN pip install --no-cache-dir -e .


FROM mambaorg/micromamba:1.5.10 AS runtime

ARG ENV_NAME=cryptic-ip
ENV ENV_NAME=${ENV_NAME}
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:${PATH}
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /workspace

COPY --from=builder /opt/conda/envs/${ENV_NAME} /opt/conda/envs/${ENV_NAME}
COPY --from=builder /workspace /workspace

ENTRYPOINT ["cryptic-ip"]
CMD ["--help"]
