FROM continuumio/miniconda3

WORKDIR /opt/app/
# set the work directory

COPY setup.py main.py ./
# copy necessary files to /opt/app

RUN conda update conda && \
    conda install click && \
    conda clean --all --yes   
# install our command line parser, Click

RUN pip install .
# install setup.py

ENV PATH="/root/.local/bin:$PATH"

RUN ls -l /opt/app

ENTRYPOINT ["run-virtual-screen"]
