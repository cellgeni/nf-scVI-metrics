FROM python:3.11

# use bash shell
SHELL ["/bin/bash", "-c"]

# non-interactive mode
# make apt-get have zero interaction while installing or upgrading
ENV DEBIAN_FRONTEND=noninteractive
ENV VENV_PATH="/env"
ENV PATH="${VENV_PATH}/bin:$PATH"

# Update the package list and install necessary packages
RUN apt-get update -y \
    && apt-get install -y\
    curl \
    wget \
    && apt-get clean

RUN python -m venv "${VENV_PATH}" && \
    . "${VENV_PATH}/bin/activate" && \
    pip install -U pip setuptools wheel && \
    pip install jupyterlab \
    bbknn \
    scib-metrics==0.5.3 \
    scvi-tools==1.3.0


# Copy Dockerfile to the container
COPY Dockerfile /docker/
RUN chmod -R 755 /docker
