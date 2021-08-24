FROM continuumio/miniconda3:4.9.2-alpine


WORKDIR /opt/app/

COPY . ./

RUN conda env update -f environment.yml && \
    conda clean --all --yes

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .

RUN /opt/app/dependencies/mgl/install.sh

ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

RUN /opt/app/dependencies/adv/bin/vina --help



