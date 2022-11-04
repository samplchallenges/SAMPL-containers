FROM continuumio/miniconda3

WORKDIR /opt/app/

COPY . ./

RUN update-ca-certificates

RUN conda update conda --all && \
    conda install python=3.8
#    conda install -c conda-forge mdtraj && \
#    conda install -c conda-forge rdkit && \
#    conda install -c openbabel openbabel && \
#    conda clean --all --yes

ENV PATH="/root/.local/bin:$PATH"

RUN /opt/app/dependencies/mgl/install.sh

ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

RUN /opt/app/dependencies/adv/bin/vina --help
