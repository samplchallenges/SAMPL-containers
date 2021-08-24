FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY . ./

RUN conda env update -f environment.yml && \
    conda clean --all --yes

RUN pip install .

#ENTRYPOINT ["conda", "run", "-n", "oepy37", "oedock"]
