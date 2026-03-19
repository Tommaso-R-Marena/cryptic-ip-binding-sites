FROM ubuntu:22.04

# Avoid tzdata interactive prompt during package installation
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    fpocket \
    apbs \
    pdb2pqr \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Establish secure non-root user execution context
RUN useradd -m pipelineuser

WORKDIR /app

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt pydantic

COPY . .

# Transfer ownership and swap users
RUN chown -R pipelineuser:pipelineuser /app
USER pipelineuser

ENTRYPOINT ["python3", "-m", "pipeline"]
CMD ["--help"]
