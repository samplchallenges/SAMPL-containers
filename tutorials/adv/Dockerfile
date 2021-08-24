FROM osatom/adv-rdkit-base:latest

WORKDIR /opt/app/

ENV PATH="/root/.local/bin:$PATH"

COPY setup.py autodock.py main.py ./

RUN pip install .

ENTRYPOINT ["run-autodock"]

